/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MNDOMethod.h"
#include "MNDORepulsionEnergy.h"
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

MNDOMethod::MNDOMethod() : ScfMethod(true, Utils::derivOrder::two, true) {
  mndoSettings_ = std::make_unique<NDDOInitializer>(BasisFunctions::sp, false);
  overlapCalculator_ =
      std::make_unique<OverlapMatrix>(elementTypes_, positions_, aoIndexes_, mndoSettings_->getElementParameters());
  mndoFock_ = std::make_shared<FockMatrix>(elementTypes_, positions_, densityMatrix_,
                                           mndoSettings_->getOneCenterIntegrals(), mndoSettings_->getElementParameters(),
                                           aoIndexes_, *overlapCalculator_, unrestrictedCalculationRunning_);
  rep_ = std::make_unique<MNDORepulsionEnergy>(elementTypes_, positions_, mndoSettings_->getElementParameters());
  densityMatrixGuess_ = std::make_unique<NDDODensityGuess>(elementTypes_, mndoSettings_->getElementParameters(),
                                                           *overlapCalculator_, nElectrons_, nAOs_);

  electronicPart_ = mndoFock_;
  initializer_ = mndoSettings_;
}

MNDOMethod::~MNDOMethod() = default;

void MNDOMethod::setStructure(const Utils::AtomCollection& atoms, std::string parameterPath) {
  readParameters(std::move(parameterPath));
  setAtomCollection(atoms);
  initialize();
}

void MNDOMethod::readParameters(const std::string& parameterPath) {
  mndoSettings_->readParameters(parameterPath);
}

void MNDOMethod::saveParameters(const std::string& fileName) {
  mndoSettings_->saveParameters(fileName);
}

RawParametersContainer& MNDOMethod::getRawParameters() {
  return mndoSettings_->getRawParameters();
}

const RawParametersContainer& MNDOMethod::getRawParameters() const {
  return mndoSettings_->getRawParameters();
}

const nddo::OneElectronMatrix& MNDOMethod::getOneElectronMatrix() const {
  return mndoFock_->getOneElectronMatrix();
}

const nddo::TwoElectronMatrix& MNDOMethod::getTwoElectronMatrix() const {
  return mndoFock_->getTwoElectronMatrix();
}
} // namespace nddo
} // namespace Sparrow
} // namespace Scine
