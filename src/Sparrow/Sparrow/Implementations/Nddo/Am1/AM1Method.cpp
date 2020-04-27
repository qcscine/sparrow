/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AM1Method.h"
#include "AM1RepulsionEnergy.h"
#include <Sparrow/Implementations/Nddo/Utils/FockMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/OverlapMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDODensityGuess.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOElectronicEnergyCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOInitializer.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <Utils/Scf/MethodExceptions.h>
#include <utility>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {

AM1Method::AM1Method() : ScfMethod(true, Utils::derivOrder::two, true) {
  am1Settings_ = std::make_unique<NDDOInitializer>(BasisFunctions::sp, false);
  overlapCalculator_ =
      std::make_unique<OverlapMatrix>(elementTypes_, positions_, aoIndexes_, am1Settings_->getElementParameters());
  am1Fock_ = std::make_shared<FockMatrix>(elementTypes_, positions_, densityMatrix_,
                                          am1Settings_->getOneCenterIntegrals(), am1Settings_->getElementParameters(),
                                          aoIndexes_, *overlapCalculator_, unrestrictedCalculationRunning_);
  rep_ = std::make_unique<AM1RepulsionEnergy>(elementTypes_, positions_, am1Settings_->getElementParameters());
  densityMatrixGuess_ = std::make_unique<NDDODensityGuess>(elementTypes_, am1Settings_->getElementParameters(),
                                                           *overlapCalculator_, nElectrons_, nAOs_);

  electronicPart_ = am1Fock_;
  initializer_ = am1Settings_;
}

AM1Method::~AM1Method() = default;

void AM1Method::setStructure(const Utils::AtomCollection& atoms, std::string parameterPath) {
  readParameters(std::move(parameterPath));
  setAtomCollection(atoms);
  initialize();
}

void AM1Method::readParameters(const std::string& parameterPath) {
  am1Settings_->readParameters(parameterPath);
}

void AM1Method::saveParameters(const std::string& fileName) {
  am1Settings_->saveParameters(fileName);
}

RawParametersContainer& AM1Method::getRawParameters() {
  return am1Settings_->getRawParameters();
}

const RawParametersContainer& AM1Method::getRawParameters() const {
  return am1Settings_->getRawParameters();
}

const OneElectronMatrix& AM1Method::getOneElectronMatrix() const {
  return am1Fock_->getOneElectronMatrix();
}

const TwoElectronMatrix& AM1Method::getTwoElectronMatrix() const {
  return am1Fock_->getTwoElectronMatrix();
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
