/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTBCommon.h"
#include "SKAtom.h"
#include "SKPair.h"
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Math/AtomicSecondDerivativeCollection.h>
#include <Utils/MethodEssentials/util/AtomsOrbitalsIndexes.h>
#include <Utils/MethodEssentials/util/MethodExceptions.h>
#include <cmath>
#include <fstream>
#include <string>

namespace Scine {
namespace Sparrow {

namespace dftb {

DFTBCommon::DFTBCommon(const Utils::ElementTypeCollection& elements, int& nEl, int& charge,
                       AtomicParameterContainer& atomicPar, DiatomicParameterContainer& diatomicPar)
  : nElements_(110),
    atomTypePresent(nElements_, false),
    atomParameters(atomicPar),
    pairParameters(diatomicPar),
    nElectrons_(nEl),
    molecularCharge_(charge),
    elementTypes_(elements) {
}

DFTBCommon::~DFTBCommon() = default;

void DFTBCommon::reinitializeParameters() {
  atomParameters = AtomicParameterContainer(110);
  pairParameters = DiatomicParameterContainer();
}

bool DFTBCommon::unrestrictedCalculationPossible() const {
  if (dftbType_ == 0)
    return false;
  return spinPolarizedPossible;
}

void DFTBCommon::setMethodDetails(const std::string& path, unsigned dftbType) {
  path_ = Utils::NativeFilenames::addTrailingSeparator(path);
  dftbType_ = dftbType;
}

void DFTBCommon::initialize(const std::string& path, unsigned dftbType) {
  setMethodDetails(path, dftbType);
  initialize(elementTypes_);
}

void DFTBCommon::initialize(const Utils::ElementTypeCollection& elementTypes) {
  reinitializeParameters();
  nAOs = 0;
  nAtoms = static_cast<unsigned int>(elementTypes.size());
  listAtomTypes.clear();

  nInitialElectrons_ = 0;
  noInteractionEnergy = 0.0;
  atomTypePresent = std::vector<bool>(nElements_, false);

  // Look what atoms are present
  aoIndexes_ = Utils::AtomsOrbitalsIndexes(nAtoms);
  for (unsigned int i = 0; i < nAtoms; i++) {
    if (!atomTypePresent[static_cast<int>(elementTypes[i])]) {
      atomTypePresent[static_cast<int>(elementTypes[i])] = true;
      if (!atomParameters[static_cast<int>(elementTypes[i])])
        atomParameters[static_cast<int>(elementTypes[i])] = std::make_unique<SKAtom>(elementTypes[i]);
      std::string elementString = Utils::ElementInfo::symbol(elementTypes[i]);
      listAtomTypes.push_back(elementString);
    }
    auto nOrbitalsForAtom = static_cast<unsigned int>(atomParameters[static_cast<int>(elementTypes[i])]->getnAOs());
    aoIndexes_.addAtom(nOrbitalsForAtom);
    nAOs += nOrbitalsForAtom;
  }

  // Add parameters for atom pairs
  for (const auto& atom1 : listAtomTypes) {
    for (const auto& atom2 : listAtomTypes) {
      auto val1 = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atom1));
      auto val2 = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atom2));
      if (!pairParameters[val1][val2])
        pairParameters[val1][val2] =
            std::make_unique<SKPair>(atom1, atom2, atomParameters[val1].get(), atomParameters[val2].get(), path_);
    }
  }

  // Eliminate redundancy and complete them
  for (const auto& atom1 : listAtomTypes) {
    for (const auto& atom2 : listAtomTypes) {
      auto val1 = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atom1));
      auto val2 = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atom2));
      auto p1 = pairParameters[val1][val2].get();
      auto p2 = pairParameters[val2][val1].get();
      if (val1 > val2)
        std::swap(p1, p2);
      p1->complete(p2);
      p1->precalculateGammaTerms();
      p2->precalculateGammaTerms();
    }
  }

  // Get number of electrons
  coreCharges_.clear();
  for (auto e : elementTypes) {
    double coreCharge = atomParameters[static_cast<int>(e)]->getOccupation();
    coreCharges_.push_back(coreCharge);
    nInitialElectrons_ += static_cast<unsigned>(coreCharge);
    noInteractionEnergy += atomParameters[static_cast<int>(e)]->getEnergy();
  }
  nElectrons_ = nInitialElectrons_ - molecularCharge_;

  if (dftbType_ > 0)
    readSpinConstants(path_);
  if (dftbType_ == 3)
    readHubbardDerivatives(path_);
}

void DFTBCommon::readSpinConstants(const std::string& path) {
  // Get spin constants if available
  std::string spinFilename = path + "spin.dat";
  std::ifstream fSpin(spinFilename.c_str());
  fSpin.imbue(std::locale("C"));
  // Output warning if it doesn't exist
  if (!fSpin)
    ; // std::cout << "Warning: no spin constants file available at " << spinFilename << "." << std::endl;
  else {
    // std::cout << "Reading spin constants from " << spinFilename << ".\n";
    std::string atom;
    getline(fSpin, atom); // Skip first line
    double sc[3][3];
    fSpin >> atom;
    while (!fSpin.eof()) {
      fSpin >> sc[0][0] >> sc[0][1] >> sc[1][0] >> sc[1][1] >> sc[0][2] >> sc[1][2] >> sc[2][2] >> sc[2][0] >> sc[2][1];
      auto valatom = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atom));
      if (atomTypePresent[valatom])
        atomParameters[valatom]->setSpinConstants(sc);
      fSpin >> atom;
    }
  }

  // Check if SDFTB is possible
  spinPolarizedPossible = true;
  for (const auto& atomStr : listAtomTypes) {
    auto vali = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atomStr));
    if (!(atomParameters[vali]->hasSpinConstants())) {
      spinPolarizedPossible = false;
      return;
    }
  }
}

void DFTBCommon::readHubbardDerivatives(const std::string& path) {
  // Get Hubbard derivatives if available
  std::string HubbardFilename = path + "hubbard.dat";
  std::ifstream fHubbard(HubbardFilename.c_str());
  fHubbard.imbue(std::locale("C"));
  // Output warning if it doesn't exist
  if (!fHubbard)
    throw Utils::Methods::ParameterFileCannotBeOpenedException(HubbardFilename);

  std::string atom;
  fHubbard >> atom;
  while (!fHubbard.eof()) {
    double hubbard;
    fHubbard >> hubbard;
    auto valatom = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atom));
    if (atomTypePresent[valatom])
      atomParameters[valatom]->setHubbardDerivative(hubbard);
    fHubbard >> atom;
  }

  // Check if DFTB3 is allowed
  DFTB3Possible = true;
  for (const auto& atomStr : listAtomTypes) {
    auto vali = static_cast<int>(Utils::ElementInfo::elementTypeForSymbol(atomStr));
    if (!(atomParameters[vali]->hasHubbardDerivative())) {
      DFTB3Possible = false;
      // std::cout << "DFTB3 is not possible because there are no Hubbard derivatives\n";
      return;
    }
  }
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
