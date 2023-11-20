/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTBCommon.h"
#include "SKAtom.h"
#include "SKPair.h"
#include "Sparrow/Resources/Dftb/ParameterSets.h"
#include "boost/filesystem.hpp"
#include <Core/Exceptions.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/Scf/MethodExceptions.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace Scine {
namespace Sparrow {

namespace dftb {

struct IncompleteParametersException : public std::exception {
  IncompleteParametersException(int Z1, int Z2) {
    error_ = ("No parameter pair found for the element pair " + Utils::ElementInfo::symbol(Utils::ElementInfo::element(Z1)) +
              " and " + Utils::ElementInfo::symbol(Utils::ElementInfo::element(Z2)));
  }
  const char* what() const noexcept final {
    return error_.c_str();
  }

 private:
  std::string error_;
};

DFTBCommon::DFTBCommon(const Utils::ElementTypeCollection& elements, int& nEl, int& charge,
                       AtomicParameterContainer& atomicPar, DiatomicParameterContainer& diatomicPar)
  : atomTypePresent(nElements_, false),
    atomParameters(atomicPar),
    pairParameters(diatomicPar),
    nElectrons_(nEl),
    molecularCharge_(charge),
    elementTypes_(elements) {
}

DFTBCommon::~DFTBCommon() = default;

void DFTBCommon::reinitializeParameters() {
  atomParameters = AtomicParameterContainer(110);
  pairParameters.clear();
}

bool DFTBCommon::unrestrictedCalculationPossible() const {
  if (dftbType_ == 0) {
    return false;
  }

  return spinPolarizedPossible;
}

void DFTBCommon::setMethodDetails(const std::string& path, unsigned dftbType) {
  path_ = path;
  dftbType_ = dftbType;
}

void DFTBCommon::initialize(const std::string& path, unsigned dftbType) {
  setMethodDetails(path, dftbType);
  initialize(elementTypes_);
}

void DFTBCommon::initialize(const Utils::ElementTypeCollection& elementTypes) {
  elementTypes_ = elementTypes;
  reinitializeParameters();
  nAOs = 0;
  nAtoms = static_cast<unsigned>(elementTypes.size());

  nInitialElectrons_ = 0;
  noInteractionEnergy = 0.0;
  atomTypePresent = std::vector<bool>(nElements_, false);

  // Look what atoms are present
  aoIndexes_ = Utils::AtomsOrbitalsIndexes(nAtoms);
  for (unsigned i = 0; i < nAtoms; i++) {
    const unsigned elementZ = Utils::ElementInfo::Z(elementTypes[i]);
    if (!atomTypePresent[elementZ]) {
      atomTypePresent[elementZ] = true;

      if (!atomParameters[elementZ]) {
        atomParameters[elementZ] = std::make_unique<SKAtom>(elementTypes[i]);
      }
    }
    auto nOrbitalsForAtom = static_cast<unsigned>(atomParameters[elementZ]->getnAOs());
    aoIndexes_.addAtom(nOrbitalsForAtom);
    nAOs += nOrbitalsForAtom;
  }

  // Find the unique atomic numbers we need
  std::vector<int> Zs;
  Zs.reserve(elementTypes.size());
  std::transform(std::begin(elementTypes), std::end(elementTypes), std::back_inserter(Zs),
                 [](Utils::ElementType e) { return Utils::ElementInfo::Z(e); });
  std::sort(std::begin(Zs), std::end(Zs));
  Zs.erase(std::unique(std::begin(Zs), std::end(Zs)), std::end(Zs));

  ParameterSet parameters;
  if (boost::filesystem::exists(path_) && boost::filesystem::is_directory(path_)) {
    parameters = ParameterSet::collect(path_, Zs);
  }
  else {
    parameters = embeddedParameters(path_, Zs).value_or_eval([this]() -> ParameterSet {
      throw std::runtime_error("No embedded parameters named '" + path_ + "'");
      return {};
    });
  }

  if (nAtoms > 0) {
    // Create SKPairs
    for (auto& atomPair : parameters.pairData) {
      const auto& pairZ = atomPair.first;
      pairParameters.emplace(std::piecewise_construct, std::forward_as_tuple(pairZ),
                             std::forward_as_tuple(atomParameters[pairZ.first].get(),
                                                   atomParameters[pairZ.second].get(), std::move(atomPair.second)));
    }
  }

  // Check completeness of the parameters w.r.t. the elements
  for (int Z1 : Zs) {
    for (int Z2 : Zs) {
      if (pairParameters.find(std::make_pair(Z1, Z2)) == pairParameters.end()) {
        throw IncompleteParametersException(Z1, Z2);
      }
    }
  }

  for (int Z1 : Zs) {
    for (int Z2 : Zs) {
      auto& p1 = pairParameters.at(std::make_pair(Z1, Z2));
      auto& p2 = pairParameters.at(std::make_pair(Z2, Z1));
      p1.complete(&p2);
      p2.complete(&p1);
      p1.precalculateGammaTerms();
      p2.precalculateGammaTerms();
    }
  }

  // Get number of electrons
  coreCharges_.clear();
  for (auto e : elementTypes) {
    double coreCharge = atomParameters[Utils::ElementInfo::Z(e)]->getOccupation();
    coreCharges_.push_back(coreCharge);
    nInitialElectrons_ += static_cast<unsigned>(coreCharge);
    noInteractionEnergy += atomParameters[Utils::ElementInfo::Z(e)]->getEnergy();
  }
  nElectrons_ = nInitialElectrons_ - molecularCharge_;

  if (dftbType_ > 0 && parameters.spin) {
    for (auto& spinPair : parameters.spin->map) {
      if (atomParameters[spinPair.first]) {
        atomParameters[spinPair.first]->setSpinConstants(std::move(spinPair.second));
      }
    }

    spinPolarizedPossible =
        std::all_of(std::begin(elementTypes), std::end(elementTypes), [&](const Utils::ElementType element) -> bool {
          const unsigned Z = Utils::ElementInfo::Z(element);
          return atomParameters[Z]->hasSpinConstants();
        });
  }

  if (dftbType_ == 3) {
    if (!parameters.hubbard) {
      throw Core::InitializationException(
          "DFTB3 is not possible because there are no Hubbard derivatives in this parameter set.\n");
    }
    for (const auto& hubbardPair : parameters.hubbard->map) {
      if (atomParameters[hubbardPair.first]) {
        atomParameters[hubbardPair.first]->setHubbardDerivative(hubbardPair.second);
      }
    }

    DFTB3Possible = std::all_of(std::begin(elementTypes), std::end(elementTypes), [&](const Utils::ElementType element) -> bool {
      const unsigned Z = Utils::ElementInfo::Z(element);
      return atomParameters[Z]->hasHubbardDerivative();
    });
    if (!DFTB3Possible) {
      throw Core::InitializationException(
          "DFTB3 is not possible because there are no Hubbard derivatives for some element.\n");
    }
  }
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
