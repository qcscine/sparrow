/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Dftb/Dftb0/DFTB0.h>
#include <Sparrow/Implementations/Dftb/Dftb3/DFTB3.h>
#include <Sparrow/Implementations/GenericMethodWrapper.h>
#include <Sparrow/StatesHandling/SparrowState.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Sparrow {

class ADFTBStatesHandlerTest : public Test {
 public:
  std::shared_ptr<Core::Calculator> interfaceDFTB0_;
  std::shared_ptr<Core::Calculator> interfaceDFTB3_;
  dftb::DFTB0 underlyingDFTB0_;
  dftb::DFTB3 underlyingDFTB3_;

  Utils::AtomCollection ethanol_;

  Utils::StatesHandler statesHandlerDftb0, statesHandlerDftb3;

  Core::Log log;

  void SetUp() override {
    log = Core::Log::silent();
#pragma omp critical
    {
      auto& manager = Core::ModuleManager::getInstance();
      interfaceDFTB0_ = manager.get<Core::Calculator>("DFTB0");
      interfaceDFTB0_->setLog(log);
      interfaceDFTB3_ = manager.get<Core::Calculator>("DFTB3");
      interfaceDFTB3_->setLog(log);
      interfaceDFTB0_->setRequiredProperties(Utils::Property::Energy);
      interfaceDFTB3_->setRequiredProperties(Utils::Property::Energy);

      std::stringstream ss("9\n\n"
                           "H      1.9655905060   -0.0263662325    1.0690084915\n"
                           "C      1.3088788172   -0.0403821764    0.1943189946\n"
                           "H      1.5790293586    0.8034866305   -0.4554748131\n"
                           "H      1.5186511399   -0.9518066799   -0.3824432806\n"
                           "C     -0.1561112248    0.0249676675    0.5877379610\n"
                           "H     -0.4682794700   -0.8500294693    1.1854276282\n"
                           "H     -0.4063173598    0.9562730342    1.1264955766\n"
                           "O     -0.8772416674    0.0083263307   -0.6652828084\n"
                           "H     -1.8356000997    0.0539308952   -0.5014877498\n");
      ethanol_ = Utils::XyzStreamHandler::read(ss);

      interfaceDFTB0_->setStructure(ethanol_);
      interfaceDFTB3_->setStructure(ethanol_);
      underlyingDFTB0_.setAtomCollection(ethanol_);
      underlyingDFTB3_.setAtomCollection(ethanol_);
      underlyingDFTB0_.initializeFromParameterPath("3ob-3-1");
      underlyingDFTB3_.initializeFromParameterPath("3ob-3-1");

      statesHandlerDftb0 = Utils::StatesHandler(interfaceDFTB0_);
      statesHandlerDftb3 = Utils::StatesHandler(interfaceDFTB3_);
    }
  }
};

TEST_F(ADFTBStatesHandlerTest, CanSaveStateWithDFTB0) {
  underlyingDFTB0_.calculate(Utils::Derivative::None, log);
  auto results = interfaceDFTB0_->calculate("");

  statesHandlerDftb0.store();
  auto savedState = std::dynamic_pointer_cast<SparrowState>(statesHandlerDftb0.popNewestState());
  const auto& stateDensity = savedState->getDensityMatrix().restrictedMatrix();
  const auto& actualDensity = underlyingDFTB0_.getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < actualDensity.rows(); ++i) {
    for (int j = i; j < actualDensity.cols(); ++j) {
      ASSERT_THAT(stateDensity.col(j)(i), DoubleNear(actualDensity.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ADFTBStatesHandlerTest, LoadsStateCorrectlyWithDFTB0) {
  underlyingDFTB0_.calculate(Utils::Derivative::None, log);
  interfaceDFTB0_->calculate("");

  Utils::DensityMatrix densityMatrix = underlyingDFTB0_.getDensityMatrix();

  auto state = std::make_shared<SparrowState>(densityMatrix);

  statesHandlerDftb0.load(state);
  statesHandlerDftb0.store();

  auto loadedDensityMatrix =
      std::dynamic_pointer_cast<SparrowState>(statesHandlerDftb0.popNewestState())->getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < loadedDensityMatrix.rows(); ++i) {
    for (int j = i; j < loadedDensityMatrix.cols(); ++j) {
      ASSERT_THAT(loadedDensityMatrix.col(j)(i), DoubleNear(densityMatrix.restrictedMatrix().col(j)(i), 1e-7));
    }
  }
}

TEST_F(ADFTBStatesHandlerTest, CanSaveStateWithDFTB3) {
  underlyingDFTB3_.calculate(Utils::Derivative::None, log);
  auto results = interfaceDFTB3_->calculate("");

  statesHandlerDftb3.store();
  auto savedState = std::dynamic_pointer_cast<SparrowState>(statesHandlerDftb3.popNewestState());
  const Eigen::MatrixXd& stateDensity = savedState->getDensityMatrix().restrictedMatrix();
  const Eigen::MatrixXd& actualDensity = underlyingDFTB3_.getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < actualDensity.rows(); ++i) {
    for (int j = i; j < actualDensity.cols(); ++j) {
      ASSERT_THAT(stateDensity.col(j)(i), DoubleNear(actualDensity.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ADFTBStatesHandlerTest, LoadsStateCorrectlyWithDFTB3) {
  underlyingDFTB3_.calculate(Utils::Derivative::None, log);

  auto densityMatrix = underlyingDFTB3_.getDensityMatrix();

  auto state = std::make_shared<SparrowState>(densityMatrix);

  statesHandlerDftb3.load(state);
  statesHandlerDftb3.store();

  Eigen::MatrixXd loadedDensityMatrix =
      std::dynamic_pointer_cast<SparrowState>(statesHandlerDftb3.popNewestState())->getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < loadedDensityMatrix.rows(); ++i) {
    for (int j = i; j < loadedDensityMatrix.cols(); ++j) {
      ASSERT_THAT(loadedDensityMatrix.col(j)(i), DoubleNear(densityMatrix.restrictedMatrix().col(j)(i), 1e-7));
    }
  }
}

TEST_F(ADFTBStatesHandlerTest, ReinitializesStateCorrectlyWithDFTB3) {
  auto& handler = statesHandlerDftb3;
  handler.store();
  interfaceDFTB3_->calculate("");
  auto oldState = handler.getState(0);
  const Eigen::MatrixXd& initialDensityMatrix =
      std::dynamic_pointer_cast<SparrowState>(oldState)->getDensityMatrix().restrictedMatrix();

  auto reinitializedState = std::make_shared<SparrowState>(
      std::dynamic_pointer_cast<GenericMethodWrapper>(interfaceDFTB3_)->getDensityMatrixGuess());
  handler.load(reinitializedState);
  handler.store();
  const Eigen::MatrixXd& endDensityMatrix =
      std::dynamic_pointer_cast<SparrowState>(handler.getState(1))->getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < initialDensityMatrix.rows(); ++i) {
    for (int j = i; j < initialDensityMatrix.cols(); ++j) {
      ASSERT_THAT(endDensityMatrix.col(j)(i), DoubleNear(initialDensityMatrix.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ADFTBStatesHandlerTest, GetsCurrentStateCorrectly) {
  auto& handler = statesHandlerDftb3;
  interfaceDFTB3_->calculate("");
  handler.store();

  const Eigen::MatrixXd& endDensityMatrix =
      std::dynamic_pointer_cast<SparrowState>(handler.getState(0))->getDensityMatrix().restrictedMatrix();

  const auto currentState = interfaceDFTB3_->getState();
  const Eigen::MatrixXd& currentDM =
      std::dynamic_pointer_cast<SparrowState>(currentState)->getDensityMatrix().restrictedMatrix();
  for (int i = 0; i < endDensityMatrix.rows(); ++i) {
    for (int j = i; j < endDensityMatrix.cols(); ++j) {
      ASSERT_THAT(currentDM.col(j)(i), DoubleNear(endDensityMatrix.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ADFTBStatesHandlerTest, CanUseStatesToCloneCalculator) {
  interfaceDFTB0_->calculate("");
  auto state = interfaceDFTB0_->getState();

  auto clonedCalculator = interfaceDFTB0_->clone();
  clonedCalculator->loadState(state);
}
} // namespace Sparrow
} // namespace Scine
