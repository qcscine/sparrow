/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/CalculatorWithReference.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Exceptions.h>
#include <Sparrow/Implementations/Nddo/TimeDependent/LinearResponse/CISLinearResponseTimeDependentCalculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Constants.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <omp.h>
#include <chrono>

using namespace testing;

namespace Scine {
namespace Sparrow {

class ASlowCalculation : public Test {
 public:
  std::shared_ptr<Core::Calculator> methodWrapper_;
  std::shared_ptr<Core::CalculatorWithReference> excitedStatesCalculator_;
  std::shared_ptr<Core::CalculatorWithReference> polymorphicCIS;

 protected:
  void SetUp() final {
    std::stringstream ssOxadiazole("31\n\n"
                                   "C  -4.72638293   1.42028071   0.21488213\n"
                                   "C  -6.07154217   1.24300816  -0.11971503\n"
                                   "C  -6.54204856  -0.01457229  -0.47681026\n"
                                   "C  -5.66996128  -1.09762221  -0.49789970\n"
                                   "C  -4.33798633  -0.92650673  -0.14379704\n"
                                   "C  -3.84431524   0.32898246   0.22432885\n"
                                   "H  -6.75814904   2.09751616  -0.10053013\n"
                                   "H  -7.59688216  -0.15140007  -0.73854141\n"
                                   "H  -6.03399299  -2.09011559  -0.78430435\n"
                                   "H  -3.66654734  -1.79464397  -0.14191096\n"
                                   "S  -2.17490031   0.51958250   0.78273141\n"
                                   "C  -4.29262738   2.77600580   0.56133840\n"
                                   "O  -3.58479588   3.54391605  -0.34398808\n"
                                   "C  -3.36442701   4.75261104   0.25680566\n"
                                   "N  -3.91940199   4.72906507   1.46127889\n"
                                   "N  -4.49627788   3.51093434   1.65171700\n"
                                   "H  -1.33769534  -1.54441748  -0.22809850\n"
                                   "H  -2.80707987   5.52123460  -0.28003351\n"
                                   "C  -1.21113637  -0.46501565  -0.43480790\n"
                                   "C   0.22398843  -0.08543102  -0.34496922\n"
                                   "C   1.08285196  -0.79516742   0.49815207\n"
                                   "C   2.42997738  -0.45910455   0.56249805\n"
                                   "C   2.92814395   0.58402121  -0.21112038\n"
                                   "C   2.07614325   1.29295343  -1.05154120\n"
                                   "C   0.72807757   0.96108082  -1.12144039\n"
                                   "H   0.69278006  -1.61504710   1.11234037\n"
                                   "H   3.09835877  -1.01845433   1.22526917\n"
                                   "H   3.98969630   0.84730854  -0.15852854\n"
                                   "H   2.46832221   2.11470261  -1.65984566\n"
                                   "H   0.05861219   1.52405414  -1.78174762\n"
                                   "H  -1.59963993  -0.29678580  -1.45598699");
    auto structure2 = Utils::XyzStreamHandler::read(ssOxadiazole);

    polymorphicCIS = std::make_shared<CISLinearResponseTimeDependentCalculator>();
    polymorphicCIS->setLog(Core::Log::silent());

    auto& manager = Core::ModuleManager::getInstance();
    methodWrapper_ = manager.get<Core::Calculator>("PM3");
    methodWrapper_->settings().modifyDouble(Utils::SettingsNames::densityRmsdCriterion, 1e-9);
    methodWrapper_->setLog(Core::Log::silent());
    methodWrapper_->setStructure(structure2);

    excitedStatesCalculator_ = manager.get<Core::CalculatorWithReference>("CIS-NDDO");
    excitedStatesCalculator_->setLog(Core::Log::silent());
  };
};

TEST_F(ASlowCalculation, BinaryTestRHF) {
  methodWrapper_->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  excitedStatesCalculator_->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 20);
  excitedStatesCalculator_->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 60);
  excitedStatesCalculator_->setReferenceCalculator(methodWrapper_);
  excitedStatesCalculator_->referenceCalculation();

  auto results = excitedStatesCalculator_->calculate();
  Eigen::VectorXd singletEnergies(20);
  singletEnergies << 3.134, 3.502, 3.798, 3.958, 3.968, 4.105, 4.112, 4.346, 4.465, 4.737, 4.976, 5.055, 5.095, 5.208,
      5.235, 5.317, 5.415, 5.498, 5.549, 5.614;
  auto eigenvals = results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues;
  for (int i = 0; i < eigenvals.size(); ++i) {
    EXPECT_NEAR(eigenvals(i) * Utils::Constants::ev_per_hartree, singletEnergies(i), 2e-3);
  }
}
TEST_F(ASlowCalculation, BinaryTestUHF) {
  methodWrapper_->settings().modifyString(
      Utils::SettingsNames::spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  methodWrapper_->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  excitedStatesCalculator_->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 20);
  excitedStatesCalculator_->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 60);
  excitedStatesCalculator_->setReferenceCalculator(methodWrapper_);
  excitedStatesCalculator_->referenceCalculation();

  Eigen::VectorXd unrestrictedEnergies(20);
  unrestrictedEnergies << 2.037, 2.171, 2.255, 2.557, 2.750, 3.134, 3.254, 3.272, 3.291, 3.393, 3.398, 3.502, 3.648,
      3.695, 3.798, 3.822, 3.958, 3.965, 3.968, 4.047;

  auto results = excitedStatesCalculator_->calculate();
  auto eigenvals = results.get<Utils::Property::ExcitedStates>().unrestricted->eigenStates.eigenValues;
  for (int i = 0; i < eigenvals.size(); ++i) {
    EXPECT_NEAR(eigenvals(i) * Utils::Constants::ev_per_hartree, unrestrictedEnergies(i), 2e-3);
  }
}
TEST_F(ASlowCalculation, CanPerformCalculationThroughInterfaceSinglet) {
  polymorphicCIS->setReferenceCalculator(methodWrapper_);
  polymorphicCIS->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  polymorphicCIS->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  polymorphicCIS->applySettings();
  polymorphicCIS->referenceCalculation();
  auto eigenPairs = polymorphicCIS->calculate().get<Utils::Property::ExcitedStates>();
  Eigen::VectorXd singletEnergies(10);
  singletEnergies << 3.134, 3.502, 3.798, 3.958, 3.968, 4.105, 4.112, 4.346, 4.465, 4.737;

  for (int i = 0; i < eigenPairs.singlet->eigenStates.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.singlet->eigenStates.eigenValues(i) * Utils::Constants::ev_per_hartree,
                DoubleNear(singletEnergies(i), 2e-3));
  }
}

TEST_F(ASlowCalculation, CanPerformCalculationThroughInterfaceTriplet) {
  polymorphicCIS->setReferenceCalculator(methodWrapper_);
  polymorphicCIS->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  polymorphicCIS->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 20);
  polymorphicCIS->settings().modifyString(Utils::SettingsNames::spinBlock, "triplet");
  polymorphicCIS->applySettings();
  polymorphicCIS->referenceCalculation();
  auto eigenPairs = polymorphicCIS->calculate().get<Utils::Property::ExcitedStates>();
  Eigen::VectorXd tripletEnergies(10);
  tripletEnergies << 2.037, 2.171, 2.255, 2.557, 2.750, 3.254, 3.272, 3.291, 3.393, 3.398;

  for (int i = 0; i < eigenPairs.triplet->eigenStates.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.triplet->eigenStates.eigenValues(i) * Utils::Constants::ev_per_hartree,
                DoubleNear(tripletEnergies(i), 5e-3));
  }
}

TEST_F(ASlowCalculation, CanPerformUHFCalculationThroughModuleInterfaceWithMultipleCores) {
  auto& manager = Core::ModuleManager::getInstance();
  methodWrapper_->settings().modifyString(
      Utils::SettingsNames::spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  excitedStatesCalculator_ = manager.get<Core::CalculatorWithReference>("CIS-NDDO");
  excitedStatesCalculator_->setLog(Core::Log::silent());
  excitedStatesCalculator_->setReferenceCalculator(methodWrapper_);
  excitedStatesCalculator_->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 20);
  excitedStatesCalculator_->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 20);
  excitedStatesCalculator_->applySettings();
  excitedStatesCalculator_->referenceCalculation();

  omp_set_num_threads(1);
  auto eigenPairs = excitedStatesCalculator_->calculate().get<Utils::Property::ExcitedStates>();
  Eigen::VectorXd singletEnergies(20);
  singletEnergies << 2.037, 2.171, 2.255, 2.557, 2.750, 3.134, 3.254, 3.272, 3.291, 3.393, 3.398, 3.502, 3.648, 3.695,
      3.798, 3.822, 3.958, 3.965, 3.968, 4.047;
  for (int i = 0; i < eigenPairs.unrestricted->eigenStates.eigenValues.size(); ++i) {
    EXPECT_THAT(eigenPairs.unrestricted->eigenStates.eigenValues(i) * Utils::Constants::ev_per_hartree,
                DoubleNear(singletEnergies(i), 5e-3));
  }

  omp_set_num_threads(4);
  auto eigenPairs4 = excitedStatesCalculator_->calculate().get<Utils::Property::ExcitedStates>();
  for (int i = 0; i < eigenPairs4.unrestricted->eigenStates.eigenValues.size(); ++i) {
    EXPECT_THAT(eigenPairs4.unrestricted->eigenStates.eigenValues(i) * Utils::Constants::ev_per_hartree,
                DoubleNear(singletEnergies(i), 5e-3));
  }
}
} // namespace Sparrow
} // namespace Scine
