/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "parameters_location.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/GenericMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Pm6/PM6Method.h>
#include <Sparrow/Implementations/Nddo/Utils/OneElectronMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/TwoElectronMatrix.h>
#include <Sparrow/StatesHandling/SparrowState.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XYZStreamHandler.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <boost/dll/runtime_symbol_info.hpp>

using namespace testing;

namespace Scine {
namespace Sparrow {

class ANDDOStatesHandlerTest : public Test {
 public:
  std::shared_ptr<Core::Calculator> interfaceMethod_;
  nddo::PM6Method underlyingMethod_;

  Utils::AtomCollection ethanol_;

 public:
  void SetUp() override {
    auto& manager = Core::ModuleManager::getInstance();
    // Load the Sparrow module if it is not already loaded.
    try {
      auto programPath = boost::dll::program_location();
      auto libPath = programPath.parent_path() / "sparrow";
      manager.load(libPath);
    }
    catch (std::exception& e) {
    }

    interfaceMethod_ = manager.get<Core::Calculator>("PM6");
    auto& settings = interfaceMethod_->settings();
    settings.modifyString(Utils::SettingsNames::parameterRootDirectory, "");
    settings.modifyString(Utils::SettingsNames::parameterFile, parameters_pm6);
    interfaceMethod_->setRequiredProperties(Utils::Property::Energy);

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
    ethanol_ = Utils::XYZStreamHandler::read(ss);
    underlyingMethod_.setStructure(ethanol_, parameters_pm6);
    interfaceMethod_->setStructure(ethanol_);
  };
};

TEST_F(ANDDOStatesHandlerTest, CanSaveMinimalState) {
  underlyingMethod_.calculate(Utils::derivativeType::none);
  auto results = interfaceMethod_->calculate("");

  interfaceMethod_->statesHandler().store(Utils::StateSize::minimal);
  auto savedState = interfaceMethod_->statesHandler().popNewestState();
  const auto& stateDensity = savedState->getMatrixState("Density Matrix");
  const auto& actualDensity = underlyingMethod_.getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < actualDensity.rows(); ++i) {
    for (int j = i; j < actualDensity.cols(); ++j) {
      ASSERT_THAT(stateDensity.col(j)(i), DoubleNear(actualDensity.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ANDDOStatesHandlerTest, CanSaveRegularlySizedState) {
  underlyingMethod_.calculate(Utils::derivativeType::none);
  auto results = interfaceMethod_->calculate("");

  interfaceMethod_->statesHandler().store(Utils::StateSize::regular);
  auto savedState = interfaceMethod_->statesHandler().popNewestState();
  auto& stateFock = savedState->getMatrixState("Fock Matrix");
  auto& actualFock = underlyingMethod_.getFockMatrix().restrictedMatrix();
  auto& stateDensity = savedState->getMatrixState("Density Matrix");
  auto& actualDensity = underlyingMethod_.getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < actualDensity.rows(); ++i) {
    for (int j = i; j < actualDensity.cols(); ++j) {
      ASSERT_THAT(stateDensity.col(j)(i), DoubleNear(actualDensity.col(j)(i), 1e-7));
    }
  }

  for (int i = 0; i < actualFock.rows(); ++i) {
    for (int j = i; j < actualFock.cols(); ++j) {
      ASSERT_THAT(stateFock.col(j)(i), DoubleNear(actualFock.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ANDDOStatesHandlerTest, CanSaveExtensiveSizedState) {
  underlyingMethod_.calculate(Utils::derivativeType::none);
  auto results = interfaceMethod_->calculate("");

  interfaceMethod_->statesHandler().store(Utils::StateSize::extensive);
  auto savedState = interfaceMethod_->statesHandler().popNewestState();
  auto& stateFock = savedState->getMatrixState("Fock Matrix");
  auto& actualFock = underlyingMethod_.getFockMatrix().restrictedMatrix();
  auto stateDensity = savedState->getMatrixState("Density Matrix");
  auto actualDensity = underlyingMethod_.getDensityMatrix().restrictedMatrix();

  for (int i = 0; i < actualDensity.rows(); ++i) {
    for (int j = i; j < actualDensity.cols(); ++j) {
      ASSERT_THAT(stateDensity.col(j)(i), DoubleNear(actualDensity.col(j)(i), 1e-7));
    }
  }
  for (int i = 0; i < actualFock.rows(); ++i) {
    for (int j = i; j < actualFock.cols(); ++j) {
      ASSERT_THAT(stateFock.col(j)(i), DoubleNear(actualFock.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ANDDOStatesHandlerTest, LoadsStateCorrectly) {
  underlyingMethod_.calculate(Utils::derivativeType::none);

  auto densityMatrix = underlyingMethod_.getDensityMatrix();

  SparrowState state(Utils::StateSize::minimal, dynamic_cast<GenericMethodWrapper&>(*interfaceMethod_));
  state.generateDensityMatrixState(densityMatrix, false);

  interfaceMethod_->statesHandler().load(std::make_shared<SparrowState>(state));
  interfaceMethod_->statesHandler().store(Utils::StateSize::minimal);

  auto loadedDensityMatrix = interfaceMethod_->statesHandler().popNewestState()->getMatrixState("Density Matrix");

  for (int i = 0; i < loadedDensityMatrix.rows(); ++i) {
    for (int j = i; j < loadedDensityMatrix.cols(); ++j) {
      ASSERT_THAT(loadedDensityMatrix.col(j)(i), DoubleNear(densityMatrix.restrictedMatrix().col(j)(i), 1e-7));
    }
  }
}

TEST_F(ANDDOStatesHandlerTest, ReinitializesStateCorrectly) {
  auto& handler = interfaceMethod_->statesHandler();
  handler.store(Utils::StateSize::minimal);
  interfaceMethod_->calculate("");
  handler.store(Utils::StateSize::minimal);
  auto oldState = handler.getState(0);
  const auto& initialDensityMatrix = oldState->getMatrixState("Density Matrix");

  handler.getState(1)->initialize();
  const auto& endDensityMatrix = handler.getState(1)->getMatrixState("Density Matrix");

  for (int i = 0; i < initialDensityMatrix.rows(); ++i) {
    for (int j = i; j < initialDensityMatrix.cols(); ++j) {
      ASSERT_THAT(endDensityMatrix.col(j)(i), DoubleNear(initialDensityMatrix.col(j)(i), 1e-7));
    }
  }
}

TEST_F(ANDDOStatesHandlerTest, GetsCurrentStateCorrectly) {
  auto& handler = interfaceMethod_->statesHandler();
  interfaceMethod_->calculate("");
  handler.store(Utils::StateSize::minimal);

  const auto& endDensityMatrix = handler.getState(0)->getMatrixState("Density Matrix");

  const auto currentState = handler.getCurrentState(Utils::StateSize::minimal);
  const auto& currentDM = currentState->getMatrixState("Density Matrix");
  for (int i = 0; i < endDensityMatrix.rows(); ++i) {
    for (int j = i; j < endDensityMatrix.cols(); ++j) {
      ASSERT_THAT(currentDM.col(j)(i), DoubleNear(endDensityMatrix.col(j)(i), 1e-7));
    }
  }
}
} // namespace Sparrow
} // namespace Scine
