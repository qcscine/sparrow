/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PM6Method.h"
#include "PM6RepulsionEnergy.h"
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

PM6Method::PM6Method() : ScfMethod(true, Utils::derivOrder::two, true) {
  pm6Settings_ = std::make_unique<NDDOInitializer>(BasisFunctions::spd, true);
  overlapCalculator_ =
      std::make_unique<OverlapMatrix>(elementTypes_, positions_, aoIndexes_, pm6Settings_->getElementParameters());
  pm6Fock_ = std::make_shared<FockMatrix>(elementTypes_, positions_, densityMatrix_,
                                          pm6Settings_->getOneCenterIntegrals(), pm6Settings_->getElementParameters(),
                                          aoIndexes_, *overlapCalculator_, unrestrictedCalculationRunning_);
  rep_ = std::make_unique<PM6RepulsionEnergy>(elementTypes_, positions_, pm6Settings_->getElementParameters(),
                                              pm6Settings_->getElementPairParameters());
  densityMatrixGuess_ = std::make_unique<NDDODensityGuess>(elementTypes_, pm6Settings_->getElementParameters(),
                                                           *overlapCalculator_, nElectrons_, nAOs_);

  electronicPart_ = pm6Fock_;
  initializer_ = pm6Settings_;
}

PM6Method::~PM6Method() = default;

void PM6Method::setStructure(const Utils::AtomCollection& atoms, std::string parameterPath) {
  readParameters(std::move(parameterPath));
  setAtomCollection(atoms);
  initialize();
}

void PM6Method::readParameters(const std::string& parameterPath) {
  pm6Settings_->readParameters(parameterPath);
}

void PM6Method::saveParameters(const std::string& fileName) {
  pm6Settings_->saveParameters(fileName);
}

RawParametersContainer& PM6Method::getRawParameters() {
  return pm6Settings_->getRawParameters();
}

const RawParametersContainer& PM6Method::getRawParameters() const {
  return pm6Settings_->getRawParameters();
}

const nddo::OneElectronMatrix& PM6Method::getOneElectronMatrix() const {
  return pm6Fock_->getOneElectronMatrix();
}

const nddo::TwoElectronMatrix& PM6Method::getTwoElectronMatrix() const {
  return pm6Fock_->getTwoElectronMatrix();
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
