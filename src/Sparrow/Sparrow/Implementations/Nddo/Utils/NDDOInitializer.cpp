/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NDDOInitializer.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParameterProcessor.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/Typenames.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

namespace nddo {

void NDDOInitializer::applyRawParameters(const Utils::ElementTypeCollection& elements) {
  elementParameters_.clear();
  elementPairParameters_.clear();
  oneCenterIntegrals_.clear();

  RawParameterProcessor processor(rawParameters_, basisFunctions_);

  // Fill atomic parameters
  for (auto e : elements) {
    if (!elementParameters_.isSet(e)) {
      if (rawParameters_.atomic.count(Utils::ElementInfo::Z(e)) == 0) {
        throw Utils::Methods::ParametersDoNotExistForElementException(e);
      }

      auto par = processor.processAtomicParameters(e);
      elementParameters_.set(e, std::move(par.first));
      oneCenterIntegrals_.set(e, std::move(par.second));
    }
  }

  // Fill pairwise parameters
  if (hasDiatomicParameters_) {
    for (int i = 0; i < 110; i++) {
      for (int j = 0; j <= i; j++) {
        auto e1 = Utils::ElementInfo::element(i);
        auto e2 = Utils::ElementInfo::element(j);
        if (elementParameters_.isSet(e1) && elementParameters_.isSet(e2)) {
          if (rawParameters_.diatomic.count(Parameters::key(e1, e2)) == 0) {
            throw Utils::Methods::ParametersDoNotExistForElementPairException(e1, e2);
          }

          elementPairParameters_.set(e1, e2, processor.runtimeDiatomicParameters(e1, e2));
        }
      }
    }
  }
}

void NDDOInitializer::readParameters(const std::string& parameterPath) {
  rawParameters_ = Parameters::read(parameterPath);
}

void NDDOInitializer::saveParameters(const std::string& fileName) {
  rawParameters_.write(fileName);
}

void NDDOInitializer::initialize(const Utils::ElementTypeCollection& elements) {
  applyRawParameters(elements);
  nElectronsForUnchargedSpecies_ = 0;
  coreCharges_.clear();
  aoIndexes_ = Utils::AtomsOrbitalsIndexes(elements.size());
  for (auto element : elements) {
    unsigned int nOrbitalsForAtom = elementParameters_.get(element).nAOs();
    aoIndexes_.addAtom(nOrbitalsForAtom);
    double coreCharge = elementParameters_.get(element).coreCharge();
    nElectronsForUnchargedSpecies_ += static_cast<unsigned int>(std::lround(coreCharge));
    coreCharges_.push_back(coreCharge);
  }
}

Parameters& NDDOInitializer::getRawParameters() {
  return rawParameters_;
}

const Parameters& NDDOInitializer::getRawParameters() const {
  return rawParameters_;
}

Utils::AtomsOrbitalsIndexes NDDOInitializer::getAtomsOrbitalsIndexes() const {
  return aoIndexes_;
}

unsigned NDDOInitializer::getNumberElectronsForUnchargedSpecies() const {
  return nElectronsForUnchargedSpecies_;
}

std::vector<double> NDDOInitializer::getCoreCharges() const {
  return coreCharges_;
}

const ElementParameters& NDDOInitializer::getElementParameters() {
  return elementParameters_;
}

const ElementPairParameters& NDDOInitializer::getElementPairParameters() {
  return elementPairParameters_;
}

const OneCenterIntegralContainer& NDDOInitializer::getOneCenterIntegrals() {
  return oneCenterIntegrals_;
}

bool NDDOInitializer::unrestrictedCalculationPossible() const {
  return true;
}

BasisFunctions NDDOInitializer::getBasisFunctions() const {
  return basisFunctions_;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
