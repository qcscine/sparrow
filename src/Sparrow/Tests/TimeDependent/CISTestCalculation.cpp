/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Sparrow/Implementations/Exceptions.h>
#include <Sparrow/Implementations/Nddo/Am1/Wrapper/AM1TypeMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h>
#include <Sparrow/Implementations/Nddo/TimeDependent/LinearResponse/CISLinearResponseTimeDependentCalculator.h>
#include <Sparrow/Implementations/TimeDependent/DiagonalPreconditionerEvaluator.h>
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

/**
 * Reference values for this test are calculated with the CIS module of the ORCA 4.1.0
 * software.
 * Example input file:
 *
 *  !PM3 RHF PRINTMOS
 *  %cis
 *  NRoots 20
 *  triplets true
 *  end
 *
 *  *xyz 0 1
 *  C     0.000000     0.000000     0.000000
 *  O     0.000000     0.000000     1.212200
 *  H     0.937197     0.000000    -0.584262
 *  H    -0.937197     0.000000    -0.584262
 *  *
 *
 */
class ACISTestCalculation : public Test {
 public:
  std::shared_ptr<Core::Calculator> calculator;
  std::shared_ptr<Core::Calculator> calculatorH;
  std::shared_ptr<Core::Calculator> calculatorHF;
  std::shared_ptr<Core::Calculator> calculatorCH3;
  std::shared_ptr<Core::Calculator> calculatorPor;
  std::shared_ptr<Core::Calculator> calculatorOxa;
  std::shared_ptr<Core::Calculator> methodWrapper_;
  std::shared_ptr<Core::CalculatorWithReference> polymorphicCIS;
  std::shared_ptr<Core::CalculatorWithReference> interfaceCIS;
  CISLinearResponseTimeDependentCalculator CISCalculator;

 protected:
  void SetUp() final {
    std::stringstream ss("4\n\n"
                         "C     0.000000     0.000000     0.000000\n"
                         "O     0.000000     0.000000     1.212200\n"
                         "H     0.937197     0.000000    -0.584262\n"
                         "H    -0.937197     0.000000    -0.584262");
    std::stringstream ssHH("2\n\n"
                           "H     0.11827444     0.000000    0\n"
                           "H     0.88172556     0.000000    0");
    std::stringstream ssHF("2\n\n"
                           "H    0.01720534      0.00000     0.00000000\n"
                           "F    0.98279466      0.00000     0.00000000");
    std::stringstream ssPorphyrine("36\n\n"
                                   "N  -0.3782    1.9570   -0.0014\n"
                                   "N   1.9571    0.3783    0.0000\n"
                                   "N  -1.9571   -0.3782    0.0007\n"
                                   "N   0.3782   -1.9570    0.0017\n"
                                   "C   0.4931    2.9223   -0.0010\n"
                                   "C   2.9224   -0.4930    0.0009\n"
                                   "C  -2.9223    0.4931    0.0003\n"
                                   "C  -0.4931   -2.9224    0.0002\n"
                                   "C  -1.6076    2.5726   -0.0008\n"
                                   "C   2.5726    1.6075    0.0001\n"
                                   "C  -2.5725   -1.6075    0.0007\n"
                                   "C   1.6075   -2.5726    0.0002\n"
                                   "C   1.9692    2.8039   -0.0005\n"
                                   "C   2.8039   -1.9692    0.0003\n"
                                   "C  -2.8040    1.9692   -0.0001\n"
                                   "C  -1.9693   -2.8039    0.0004\n"
                                   "C  -0.0844    4.2732   -0.0003\n"
                                   "C   4.2734    0.0844    0.0012\n"
                                   "C  -4.2733   -0.0843    0.0005\n"
                                   "C   0.0844   -4.2734   -0.0023\n"
                                   "C  -1.3904    3.9982   -0.0001\n"
                                   "C   3.9982    1.3903    0.0007\n"
                                   "C  -3.9982   -1.3903    0.0007\n"
                                   "C   1.3904   -3.9983   -0.0022\n"
                                   "H   2.5488    3.7221   -0.0005\n"
                                   "H   3.7220   -2.5487   -0.0012\n"
                                   "H  -3.7221    2.5488    0.0007\n"
                                   "H  -2.5489   -3.7221    0.0002\n"
                                   "H   0.4202    5.2213    0.0003\n"
                                   "H   5.2214   -0.4202    0.0017\n"
                                   "H  -5.2213    0.4203    0.0005\n"
                                   "H  -0.4202   -5.2214   -0.0039\n"
                                   "H  -2.1793    4.7348    0.0005\n"
                                   "H   4.7349    2.1792    0.0007\n"
                                   "H  -4.7349   -2.1793    0.0009\n"
                                   "H   2.1793   -4.7349   -0.0038");
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
    auto structure = Utils::XyzStreamHandler::read(ss);
    calculator = std::make_shared<PM3MethodWrapper>();
    calculator->setLog(Core::Log::silent());
    calculator->setStructure(structure);
    calculatorH = calculator->clone();
    calculatorH->setStructure(Utils::XyzStreamHandler::read(ssHH));
    calculatorHF = calculator->clone();
    calculatorHF->setStructure(Utils::XyzStreamHandler::read(ssHF));
    calculatorPor = calculator->clone();
    calculatorPor->setStructure(Utils::XyzStreamHandler::read(ssPorphyrine));
    calculatorPor->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
    calculatorOxa = calculator->clone();
    auto structure2 = Utils::XyzStreamHandler::read(ssOxadiazole);
    calculatorOxa->setStructure(structure2);
    calculatorOxa->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
    polymorphicCIS = std::make_shared<CISLinearResponseTimeDependentCalculator>();

    CISCalculator.setLog(Core::Log::silent());
    polymorphicCIS->setLog(Core::Log::silent());
  };
};

TEST_F(ACISTestCalculation, CanAssignCalculator) {
  CISCalculator.setReferenceCalculator(calculator);
}

TEST_F(ACISTestCalculation, ThrowsIfMethodIsNotNDDO) {
  std::shared_ptr<Core::Calculator> dftbCalculator = std::make_shared<DFTB2MethodWrapper>();
  ASSERT_THROW(CISCalculator.setReferenceCalculator(dftbCalculator), InvalidCalculatorTypeForCIS);
}

TEST_F(ACISTestCalculation, CanPerformReferenceCalculation) {
  CISCalculator.setReferenceCalculator(calculator->clone());
  CISCalculator.referenceCalculation();
  calculator->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 1);
  calculator->setRequiredProperties(Utils::Property::Energy);
  calculator->calculate("CIS reference calculation.");
  ASSERT_THAT(calculator->results().get<Utils::Property::Energy>(),
              DoubleNear(CISCalculator.getReferenceCalculator().results().get<Utils::Property::Energy>(), 1e-6));
}

TEST_F(ACISTestCalculation, ThrowsIfNoReferenceCalculationPerformed) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "singlet");
  CISCalculator.setReferenceCalculator(calculator->clone());
  ASSERT_THROW(CISCalculator.calculate(), InvalidReferenceCalculationException);
}

TEST_F(ACISTestCalculation, RHFCISPreconditionerEvaluatedCorrectly) {
  Utils::SingleParticleEnergies energies;
  auto mos = Utils::MolecularOrbitals::createFromRestrictedCoefficients(Eigen::MatrixXd::Random(5, 5));
  Utils::LcaoUtils::ElectronicOccupation occupation;
  occupation.fillLowestRestrictedOrbitalsWithElectrons(4);
  Eigen::VectorXd arbitraryEnergies(5);
  arbitraryEnergies << -3.2, -1.0, 2.4, 5.7, 6.2;
  energies.setRestricted(arbitraryEnergies);
  Eigen::VectorXd energyDifferences(6);
  energyDifferences << 5.6, 8.9, 9.4, 3.4, 6.7, 7.2;

  Eigen::VectorXd diff1 = energyDifferences.array() - 1;
  Eigen::VectorXd diff2 = energyDifferences.array() - 2;
  Eigen::VectorXd invDiff1 = diff1.cwiseInverse();
  Eigen::VectorXd invDiff2 = diff2.cwiseInverse();

  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDifferenceVector;
  TimeDependentUtils::generateEnergyDifferenceVector(energies, occupation, energyDifferenceVector);

  DiagonalPreconditionerEvaluator eval(energyDifferenceVector);
  auto result1 = eval.evaluate(Eigen::VectorXd::Ones(energyDifferenceVector.restricted.size()), 1);
  auto result2 = eval.evaluate(Eigen::VectorXd::Ones(energyDifferenceVector.restricted.size()), 2);

  for (int i = 0; i < result1.size(); ++i) {
    EXPECT_THAT(-result1(i), DoubleNear(invDiff1(i), 1e-6));
    EXPECT_THAT(-result2(i), DoubleNear(invDiff2(i), 1e-6));
  }

  Eigen::VectorXd orderedEnDiffs(6);
  orderedEnDiffs << 3.4, 5.6, 6.7, 7.2, 8.9, 9.4;
  diff1 = orderedEnDiffs.array() - 1;
  diff2 = orderedEnDiffs.array() - 2;
  invDiff1 = diff1.cwiseInverse();
  invDiff2 = diff2.cwiseInverse();

  DiagonalPreconditionerEvaluator evalOrdered(energyDifferenceVector, OrderTag{});
  result1 = evalOrdered.evaluate(Eigen::VectorXd::Ones(energyDifferenceVector.restricted.size()), 1);
  result2 = evalOrdered.evaluate(Eigen::VectorXd::Ones(energyDifferenceVector.restricted.size()), 2);

  for (int i = 0; i < result1.size(); ++i) {
    EXPECT_THAT(-result1(i), DoubleNear(invDiff1(i), 1e-6));
    EXPECT_THAT(-result2(i), DoubleNear(invDiff2(i), 1e-6));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnSingletBlockForHF) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, Utils::SettingsNames::SpinBlocks::singlet);

  CISCalculator.setReferenceCalculator(calculatorHF->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 4);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 8);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().singlet->eigenStates;
  Eigen::VectorXd singletEnergies(4);
  singletEnergies << 8.178, 8.178, 8.741, 17.759;

  for (int i = 0; i < eigenPairs.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 5e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnSingletBlockForHFWithUHF) {
  CISCalculator.setReferenceCalculator(calculatorHF->clone());
  CISCalculator.getReferenceCalculator().settings().modifyString(
      Utils::SettingsNames::spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 4);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 8);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().unrestricted->eigenStates;
  Eigen::VectorXd singletEnergies(4);
  singletEnergies << 7.832, 7.832, 8.178, 8.178;

  for (int i = 0; i < eigenPairs.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 5e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnSingletBlockForHydrogen) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "singlet");
  CISCalculator.setReferenceCalculator(calculatorH->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 1);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 1);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().singlet->eigenStates;
  Eigen::VectorXd singletEnergies(1);
  singletEnergies << 10.061;

  for (int i = 0; i < eigenPairs.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 5e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnTripletBlockForHydrogen) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "triplet");
  CISCalculator.setReferenceCalculator(calculatorH->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 1);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 1);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().triplet->eigenStates;
  Eigen::VectorXd tripletEnergies(1);
  // From Orca PM3
  tripletEnergies << 6.908;

  for (int i = 0; i < eigenPairs.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(tripletEnergies(i), 1e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnSingletBlock) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "singlet");
  CISCalculator.setReferenceCalculator(calculator->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().singlet->eigenStates;
  Eigen::VectorXd singletEnergies(10);
  singletEnergies << 2.716, 5.476, 6.506, 7.986, 8.526, 9.111, 9.114, 9.606, 10.719, 11.041;

  for (int i = 0; i < eigenPairs.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 2e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnTripletBlock) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "triplet");
  CISCalculator.setReferenceCalculator(calculator->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().triplet->eigenStates;
  Eigen::VectorXd tripletEnergies(10);
  tripletEnergies << 2.341, 4.316, 4.868, 5.706, 7.094, 7.908, 8.533, 9.037, 10.271, 10.369;

  for (int i = 0; i < eigenPairs.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(tripletEnergies(i), 5e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnSingletBlockForBiggerSystem) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "singlet");
  CISCalculator.setReferenceCalculator(calculatorPor->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().singlet->eigenStates;
  Eigen::VectorXd singletEnergies(10);
  singletEnergies << 2.178, 2.346, 2.537, 2.537, 2.769, 3.183, 3.183, 3.216, 3.503, 3.868;

  for (int i = 0; i < 3; ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 2e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnSingletBlockForOxadiazoleRHF) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "singlet");
  CISCalculator.setReferenceCalculator(calculatorOxa->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 20);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().singlet->eigenStates;
  Eigen::VectorXd singletEnergies(20);
  singletEnergies << 3.134, 3.502, 3.798, 3.958, 3.968, 4.105, 4.112, 4.346, 4.465, 4.737, 4.976, 5.055, 5.095, 5.208,
      5.235, 5.317, 5.415, 5.498, 5.549, 5.614;
  for (int i = 0; i < 3; ++i) {
    EXPECT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 2e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnTripletBlockForOxadiazoleRHF) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "triplet");
  CISCalculator.setReferenceCalculator(calculatorOxa->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 20);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().triplet->eigenStates;
  Eigen::VectorXd singletEnergies(20);
  singletEnergies << 2.037, 2.171, 2.255, 2.557, 2.750, 3.254, 3.272, 3.291, 3.393, 3.398, 3.648, 3.695, 3.822, 3.965,
      4.047, 4.337, 4.583, 4.607, 4.639, 4.788;
  for (int i = 0; i < 3; ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(singletEnergies(i), 2e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationForOxadiazoleUHF) {
  CISCalculator.setReferenceCalculator(calculatorOxa->clone());
  CISCalculator.getReferenceCalculator().settings().modifyString(
      Utils::SettingsNames::spinMode, Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 20);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().unrestricted->eigenStates;
  // put 20 just to know all the data
  Eigen::VectorXd unrestrictedEnergies(20);
  unrestrictedEnergies << 2.037, 2.171, 2.255, 2.557, 2.750, 3.134, 3.254, 3.272, 3.291, 3.393, 3.398, 3.502, 3.648,
      3.695, 3.798, 3.822, 3.958, 3.965, 3.968, 4.047;
  for (int i = 0; i < 3; ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(unrestrictedEnergies(i), 1.5e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformDenseCalculationOnTripletBlockForBiggerSystem) {
  CISCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, "triplet");
  CISCalculator.setReferenceCalculator(calculatorPor->clone());
  CISCalculator.referenceCalculation();
  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 20);
  auto eigenPairs = CISCalculator.calculate().get<Utils::Property::ExcitedStates>().triplet->eigenStates;
  Eigen::VectorXd tripletEnergies(20);
  tripletEnergies << 1.163, 1.561, 1.561, 1.849, 2.026, 2.053, 2.076, 2.076, 2.178, 2.178, 2.219, 2.340, 3.068, 3.318,
      3.318, 3.497, 4.196, 4.218, 4.218, 4.297;

  for (int i = 0; i < 3; ++i) {
    ASSERT_THAT(eigenPairs.eigenValues(i) * Utils::Constants::ev_per_hartree, DoubleNear(tripletEnergies(i), 2e-3));
  }
}

TEST_F(ACISTestCalculation, CanAskForCISThroughInterface) {
  auto& manager = Core::ModuleManager::getInstance();

  interfaceCIS = manager.get<Core::CalculatorWithReference>("CIS-NDDO");
}

TEST_F(ACISTestCalculation, CanPerformCalculationThroughModuleInterface) {
  auto& manager = Core::ModuleManager::getInstance();

  interfaceCIS = manager.get<Core::CalculatorWithReference>("CIS-NDDO");
  interfaceCIS->setLog(Core::Log::silent());
  interfaceCIS->setReferenceCalculator(calculator->clone());
  interfaceCIS->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  interfaceCIS->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  interfaceCIS->settings().modifyString(Utils::SettingsNames::spinBlock, "singlet");
  interfaceCIS->applySettings();
  interfaceCIS->referenceCalculation();
  auto eigenPairs = interfaceCIS->calculate().get<Utils::Property::ExcitedStates>();

  Eigen::VectorXd singletEnergies(10);
  singletEnergies << 2.716, 5.476, 6.506, 7.986, 8.526, 9.111, 9.114, 9.606, 10.719, 11.041;

  for (int i = 0; i < eigenPairs.singlet->eigenStates.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.singlet->eigenStates.eigenValues(i) * Utils::Constants::ev_per_hartree,
                DoubleNear(singletEnergies(i), 5e-3));
  }
}

TEST_F(ACISTestCalculation, CanPerformUHFCalculationThroughModuleInterface) {
  auto& manager = Core::ModuleManager::getInstance();

  calculator->settings().modifyString(Utils::SettingsNames::spinMode,
                                      Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  interfaceCIS = manager.get<Core::CalculatorWithReference>("CIS-NDDO");
  interfaceCIS->setLog(Core::Log::silent());
  interfaceCIS->setReferenceCalculator(calculator->clone());
  interfaceCIS->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  interfaceCIS->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  interfaceCIS->applySettings();
  interfaceCIS->referenceCalculation();
  omp_set_num_threads(1);
  auto eigenPairs = interfaceCIS->calculate().get<Utils::Property::ExcitedStates>();
  omp_set_num_threads(4);
  Eigen::VectorXd singletEnergies(10);
  singletEnergies << 2.341, 2.716, 4.316, 4.868, 5.476, 5.706, 6.506, 7.094, 7.908, 7.986;

  for (int i = 0; i < eigenPairs.unrestricted->eigenStates.eigenValues.size(); ++i) {
    ASSERT_THAT(eigenPairs.unrestricted->eigenStates.eigenValues(i) * Utils::Constants::ev_per_hartree,
                DoubleNear(singletEnergies(i), 5e-3));
  }
}
TEST_F(ACISTestCalculation, GeneratesCorrectEnergyOrderMap) {
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDifferenceVector;
  energyDifferenceVector.restricted = Eigen::VectorXd(4);
  energyDifferenceVector.restricted << 3.5, 2.2, 6.4, 1.0; // Order should be 3, 1, 0, 2
  auto enMap = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);
  ASSERT_EQ(enMap.size(), 4);
  EXPECT_EQ(enMap[0], 3);
  EXPECT_EQ(enMap[1], 1);
  EXPECT_EQ(enMap[2], 0);
  EXPECT_EQ(enMap[3], 2);

  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> unrestrictedEnergyDifferenceVector;
  unrestrictedEnergyDifferenceVector.alpha = Eigen::VectorXd(4);
  unrestrictedEnergyDifferenceVector.alpha << 3.5, 2.2, 6.4, 1.0; // Order should be 3, 1, 0, 2
  unrestrictedEnergyDifferenceVector.beta = Eigen::VectorXd(4);
  unrestrictedEnergyDifferenceVector.beta << 3.2, 7.4, 1.2, -1.2; // Order should be 3, 2, 0, 1
  auto unEnMap = TimeDependentUtils::generateEnergyOrderMap(unrestrictedEnergyDifferenceVector);
  ASSERT_EQ(unEnMap.size(), 8);
  EXPECT_EQ(unEnMap[0], 7);
  EXPECT_EQ(unEnMap[1], 3);
  EXPECT_EQ(unEnMap[2], 6);
  EXPECT_EQ(unEnMap[3], 1);
  EXPECT_EQ(unEnMap[4], 4);
  EXPECT_EQ(unEnMap[5], 0);
  EXPECT_EQ(unEnMap[6], 2);
  EXPECT_EQ(unEnMap[7], 5);
}

TEST_F(ACISTestCalculation, CorrectlyTransformsToOrders) {
  std::vector<int> orderMap{3, 1, 0, 2};
  Eigen::VectorXd vecToTransform(4);
  vecToTransform << 0, 1, 2, 3;
  std::vector<int> stdVecToTransform(4);
  stdVecToTransform[0] = 0;
  stdVecToTransform[1] = 1;
  stdVecToTransform[2] = 2;
  stdVecToTransform[3] = 3;

  Eigen::VectorXd transformedVector;
  TimeDependentUtils::transformOrder(vecToTransform, transformedVector, orderMap, TimeDependentUtils::Direction::To);
  ASSERT_EQ(transformedVector.size(), 4);
  EXPECT_EQ(transformedVector(0), 3);
  EXPECT_EQ(transformedVector(1), 1);
  EXPECT_EQ(transformedVector(2), 0);
  EXPECT_EQ(transformedVector(3), 2);
  Eigen::VectorXd backTransformedVec;
  TimeDependentUtils::transformOrder(transformedVector, backTransformedVec, orderMap, TimeDependentUtils::Direction::From);
  ASSERT_EQ(backTransformedVec.size(), 4);
  EXPECT_EQ(backTransformedVec(0), vecToTransform(0));
  EXPECT_EQ(backTransformedVec(1), vecToTransform(1));
  EXPECT_EQ(backTransformedVec(2), vecToTransform(2));
  EXPECT_EQ(backTransformedVec(3), vecToTransform(3));

  std::vector<int> transformedStdVector;
  TimeDependentUtils::transformOrder(stdVecToTransform, transformedStdVector, orderMap, TimeDependentUtils::Direction::To);
  ASSERT_EQ(transformedStdVector.size(), 4);
  EXPECT_EQ(transformedStdVector[0], 3);
  EXPECT_EQ(transformedStdVector[1], 1);
  EXPECT_EQ(transformedStdVector[2], 0);
  EXPECT_EQ(transformedStdVector[3], 2);
  std::vector<int> backTransformedStdVec;
  TimeDependentUtils::transformOrder(transformedStdVector, backTransformedStdVec, orderMap, TimeDependentUtils::Direction::From);
  ASSERT_EQ(backTransformedStdVec.size(), 4);
  EXPECT_EQ(backTransformedStdVec[0], stdVecToTransform[0]);
  EXPECT_EQ(backTransformedStdVec[1], stdVecToTransform[1]);
  EXPECT_EQ(backTransformedStdVec[2], stdVecToTransform[2]);
  EXPECT_EQ(backTransformedStdVec[3], stdVecToTransform[3]);

  Eigen::MatrixXd matToTransform(4, 2);
  matToTransform << 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5;
  Eigen::MatrixXd transformedMatrix;
  TimeDependentUtils::transformOrder(matToTransform, transformedMatrix, orderMap, TimeDependentUtils::Direction::To);
  ASSERT_EQ(transformedMatrix.rows(), 4);
  ASSERT_EQ(transformedMatrix.cols(), 2);
  EXPECT_EQ(transformedMatrix(0, 0), 3);
  EXPECT_EQ(transformedMatrix(0, 1), 3.5);
  EXPECT_EQ(transformedMatrix(1, 0), 1);
  EXPECT_EQ(transformedMatrix(1, 1), 1.5);
  EXPECT_EQ(transformedMatrix(2, 0), 0);
  EXPECT_EQ(transformedMatrix(2, 1), 0.5);
  EXPECT_EQ(transformedMatrix(3, 0), 2);
  EXPECT_EQ(transformedMatrix(3, 1), 2.5);
  Eigen::MatrixXd backTransformedMat;
  TimeDependentUtils::transformOrder(transformedMatrix, backTransformedMat, orderMap, TimeDependentUtils::Direction::From);
  ASSERT_EQ(backTransformedMat.rows(), 4);
  ASSERT_EQ(backTransformedMat.cols(), 2);
  EXPECT_EQ(backTransformedMat.row(0), matToTransform.row(0));
  EXPECT_EQ(backTransformedMat.row(1), matToTransform.row(1));
  EXPECT_EQ(backTransformedMat.row(2), matToTransform.row(2));
  EXPECT_EQ(backTransformedMat.row(3), matToTransform.row(3));

  matToTransform = Eigen::MatrixXd::Random(4, 20);
  TimeDependentUtils::transformOrder(matToTransform.leftCols(2), transformedMatrix, orderMap,
                                     TimeDependentUtils::Direction::To);
  ASSERT_EQ(transformedMatrix.rows(), 4);
  ASSERT_EQ(transformedMatrix.cols(), 2);
  EXPECT_EQ(transformedMatrix.row(0), matToTransform.leftCols(2).row(3));
  EXPECT_EQ(transformedMatrix.row(1), matToTransform.leftCols(2).row(1));
  EXPECT_EQ(transformedMatrix.row(2), matToTransform.leftCols(2).row(0));
  EXPECT_EQ(transformedMatrix.row(3), matToTransform.leftCols(2).row(2));
  TimeDependentUtils::transformOrder(transformedMatrix, backTransformedMat, orderMap, TimeDependentUtils::Direction::From);
  ASSERT_EQ(backTransformedMat.rows(), 4);
  ASSERT_EQ(backTransformedMat.cols(), 2);
  EXPECT_EQ(backTransformedMat.row(0), matToTransform.leftCols(2).row(0));
  EXPECT_EQ(backTransformedMat.row(1), matToTransform.leftCols(2).row(1));
  EXPECT_EQ(backTransformedMat.row(2), matToTransform.leftCols(2).row(2));
  EXPECT_EQ(backTransformedMat.row(3), matToTransform.leftCols(2).row(3));
}

TEST_F(ACISTestCalculation, CanCalculateHydrogenFluoride) {
  calculatorHF->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-9);
  calculatorHF->settings().modifyString(Utils::SettingsNames::spinMode,
                                        Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  CISCalculator.setReferenceCalculator(calculatorHF);
  CISCalculator.referenceCalculation();

  CISCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 4);
  CISCalculator.settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 4);
  auto results = CISCalculator.calculate();
  std::vector<std::string> notOrderedLabels{"0a -> 4a", "1a -> 4a", "2a -> 4a", "3a -> 4a",
                                            "0b -> 4b", "1b -> 4b", "2b -> 4b", "3b -> 4b"};
  auto labels = results.get<Utils::Property::ExcitedStates>().transitionLabels;
  std::vector<int> order(8);
  for (unsigned int i = 0; i < notOrderedLabels.size(); ++i) {
    auto it = std::find(notOrderedLabels.begin(), notOrderedLabels.end(), labels[i]);
    order[i] = std::distance(notOrderedLabels.begin(), it);
  }

  int nAlphaExcitations = 4;
  auto eigenPairs = results.get<Utils::Property::ExcitedStates>().unrestricted->eigenStates;
  Eigen::MatrixXd reorderedEigenVectors;
  TimeDependentUtils::transformOrder(eigenPairs.eigenVectors, reorderedEigenVectors, order, TimeDependentUtils::Direction::From);
  EXPECT_THAT(reorderedEigenVectors(0, 2), DoubleNear(reorderedEigenVectors(0 + nAlphaExcitations, 2), 1e-8));
  EXPECT_THAT(reorderedEigenVectors(1, 2), DoubleNear(reorderedEigenVectors(1 + nAlphaExcitations, 2), 1e-8));
  EXPECT_THAT(reorderedEigenVectors(2, 2), DoubleNear(reorderedEigenVectors(2 + nAlphaExcitations, 2), 1e-8));
  EXPECT_THAT(reorderedEigenVectors(3, 2), DoubleNear(reorderedEigenVectors(3 + nAlphaExcitations, 2), 1e-8));
}

} // namespace Sparrow
} // namespace Scine
