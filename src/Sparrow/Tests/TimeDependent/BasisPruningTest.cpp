/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/BasisPruner.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBCalculator.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBData.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/TransitionChargesCalculator.h>
#include <Sparrow/Implementations/TimeDependent/DiagonalPreconditionerEvaluator.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Math/IterativeDiagonalizer/DiagonalizerSettings.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <chrono>

using namespace testing;

namespace Scine {
namespace Sparrow {

class ABasisPruningTest : public Test {
 public:
  std::shared_ptr<Core::Calculator> calculator;
  std::shared_ptr<Core::Calculator> calculatorHF;
  std::shared_ptr<Core::Calculator> calculatorOxa;
  TDDFTBCalculator tddftbCalculator;

 private:
  void SetUp() final {
    std::stringstream ss("4\n\n"
                         "C     0.000000     0.000000     0.000000\n"
                         "O     0.000000     0.000000     1.212200\n"
                         "H     0.937197     0.000000    -0.584262\n"
                         "H    -0.937197     0.000000    -0.584262");
    std::stringstream ssHF("2\n\n"
                           "H    0.01720534      0.00000     0.00000000\n"
                           "F    0.98279466      0.00000     0.00000000");
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
    calculator = std::make_shared<DFTB2MethodWrapper>();
    calculator->setLog(Core::Log::silent());
    calculator->settings().modifyString("method_parameters", "3ob-3-1");
    calculator->setStructure(structure);
    calculatorOxa = calculator->clone();
    calculatorHF = calculator->clone();
    auto structureHF = Utils::XyzStreamHandler::read(ssHF);
    calculatorHF->setStructure(structureHF);
    calculatorHF->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
    auto structureOxa = Utils::XyzStreamHandler::read(ssOxadiazole);
    calculatorOxa->setStructure(structureOxa);
    calculatorOxa->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-8);
    tddftbCalculator.setLog(Core::Log::silent());
  }
};

TEST_F(ABasisPruningTest, ResultsWithBasisPrunerSameWithInfiniteEnergyThresholdRestricted) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculator);
  dftbMethod->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  dftbMethod->calculate("TDDFTB reference calculation.");

  auto tddftbData = dftbMethod->getTDDFTBData();
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDifferenceVector;
  TimeDependentUtils::generateEnergyDifferenceVector(tddftbData.MOEnergies, tddftbData.occupation, energyDifferenceVector);
  TransitionChargesCalculator tcCalculator(tddftbData.molecularOrbitals, tddftbData.overlapMatrix, tddftbData.AOInfo);
  auto excitations = TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(tddftbData.molecularOrbitals,
                                                                                           tddftbData.occupation);
  Eigen::MatrixXd transitionCharges =
      tcCalculator.calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(tddftbData.occupation);

  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);

  OrderedInput<Utils::Reference::Restricted> input{energyDifferenceVector, excitations, transitionCharges, orderMap};

  BasisPruner<Utils::Reference::Restricted> basisPruner(input, tddftbData.gammaMatrix, tddftbData.spinConstants);

  auto prunedData = basisPruner.prune(EnergyThreshold{std::numeric_limits<double>::max()}, PerturbativeThreshold{0.0},
                                      Utils::SpinTransition::Singlet);

  Eigen::VectorXd orderedEnergies;
  TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(energyDifferenceVector), orderedEnergies, orderMap,
                                     TimeDependentUtils::Direction::To);

  Eigen::MatrixXd orderedTransitionCharges;
  TimeDependentUtils::transformOrder(transitionCharges, orderedTransitionCharges, orderMap, TimeDependentUtils::Direction::To);

  std::vector<Utils::Excitation> orderedExcitations;
  TimeDependentUtils::transformOrder(
      TimeDependentUtils::flatten(TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(
          tddftbData.molecularOrbitals, tddftbData.occupation)),
      orderedExcitations, orderMap, TimeDependentUtils::Direction::To);

  ASSERT_EQ(prunedData.excitations().size(), orderedExcitations.size());
  for (int i = 0; i < orderedEnergies.size(); ++i) {
    EXPECT_DOUBLE_EQ(prunedData.energyDifferences()(i), orderedEnergies(i));
    EXPECT_EQ(prunedData.excitations()[i].occ, orderedExcitations[i].occ);
    EXPECT_EQ(prunedData.excitations()[i].vir, orderedExcitations[i].vir);
    for (int atom = 0; atom < 4; ++atom)
      EXPECT_DOUBLE_EQ(prunedData.transitionCharges()(i, atom), orderedTransitionCharges(i, atom));
  }

  ASSERT_EQ(tddftbData.gammaMatrix->rows(), 4);
  ASSERT_EQ(tddftbData.gammaMatrix->cols(), 4);
  ASSERT_EQ(tddftbData.spinConstants->size(), 4);
  for (int atom_i = 0; atom_i < 4; ++atom_i) {
    for (int atom_j = 0; atom_j < 4; ++atom_j) {
      EXPECT_DOUBLE_EQ((*tddftbData.gammaMatrix)(atom_i, atom_j), (*tddftbData.gammaMatrix)(atom_i, atom_j));
    }
  }

  auto prunedDataTriplet = basisPruner.prune(EnergyThreshold{std::numeric_limits<double>::max()},
                                             PerturbativeThreshold{0.0}, Utils::SpinTransition::Triplet);
  ASSERT_EQ(prunedDataTriplet.excitations().size(), orderedExcitations.size());
  for (int i = 0; i < orderedEnergies.size(); ++i) {
    EXPECT_DOUBLE_EQ(prunedDataTriplet.energyDifferences()(i), orderedEnergies(i));
    EXPECT_EQ(prunedDataTriplet.excitations()[i].occ, orderedExcitations[i].occ);
    EXPECT_EQ(prunedDataTriplet.excitations()[i].vir, orderedExcitations[i].vir);
    for (int atom = 0; atom < 4; ++atom)
      EXPECT_DOUBLE_EQ(prunedDataTriplet.transitionCharges()(i, atom), orderedTransitionCharges(i, atom));
  }

  ASSERT_EQ(tddftbData.gammaMatrix->rows(), 4);
  ASSERT_EQ(tddftbData.gammaMatrix->cols(), 4);
  ASSERT_TRUE(tddftbData.spinConstants);
  for (int atom_i = 0; atom_i < 4; ++atom_i) {
    for (int atom_j = 0; atom_j < 4; ++atom_j) {
      EXPECT_DOUBLE_EQ((*tddftbData.gammaMatrix)(atom_i, atom_j), (*tddftbData.gammaMatrix)(atom_i, atom_j));
    }
  }
}

TEST_F(ABasisPruningTest, ResultsWithBasisPrunerSameWithInfiniteEnergyThresholdUnrestricted) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculator);
  dftbMethod->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  dftbMethod->settings().modifyString(Utils::SettingsNames::spinMode,
                                      Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  dftbMethod->calculate("TDDFTB reference calculation.");

  auto tddftbData = dftbMethod->getTDDFTBData();
  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> energyDifferenceVector;
  TimeDependentUtils::generateEnergyDifferenceVector(tddftbData.MOEnergies, tddftbData.occupation, energyDifferenceVector);
  TransitionChargesCalculator tcCalculator(tddftbData.molecularOrbitals, tddftbData.overlapMatrix, tddftbData.AOInfo);

  Eigen::MatrixXd transitionCharges =
      tcCalculator.calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(tddftbData.occupation);

  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);

  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<Utils::Excitation>> excitations =
      TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(dftbMethod->getTDDFTBData().molecularOrbitals,
                                                                              dftbMethod->getTDDFTBData().occupation);

  OrderedInput<Utils::Reference::Unrestricted> input{energyDifferenceVector, excitations, transitionCharges, orderMap};
  BasisPruner<Utils::Reference::Unrestricted> basisPruner(input, tddftbData.gammaMatrix, tddftbData.spinConstants);

  auto prunedData = basisPruner.prune(EnergyThreshold{std::numeric_limits<double>::max()}, PerturbativeThreshold{0.0});

  Eigen::VectorXd orderedEnergies;
  TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(energyDifferenceVector), orderedEnergies, orderMap,
                                     TimeDependentUtils::Direction::To);

  Eigen::MatrixXd orderedTransitionCharges;
  TimeDependentUtils::transformOrder(transitionCharges, orderedTransitionCharges, orderMap, TimeDependentUtils::Direction::To);

  std::vector<Utils::Excitation> orderedExcitations;
  TimeDependentUtils::transformOrder(
      TimeDependentUtils::flatten(TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(
          tddftbData.molecularOrbitals, tddftbData.occupation)),
      orderedExcitations, orderMap, TimeDependentUtils::Direction::To);

  ASSERT_EQ(prunedData.excitations().size(), orderedExcitations.size());
  for (int i = 0; i < orderedEnergies.size(); ++i) {
    EXPECT_DOUBLE_EQ(prunedData.energyDifferences()(i), orderedEnergies(i));
    EXPECT_EQ(prunedData.excitations()[i].occ, orderedExcitations[i].occ);
    EXPECT_EQ(prunedData.excitations()[i].vir, orderedExcitations[i].vir);
    for (int atom = 0; atom < 4; ++atom)
      EXPECT_DOUBLE_EQ(prunedData.transitionCharges()(i, atom), orderedTransitionCharges(i, atom));
  }

  ASSERT_EQ(tddftbData.gammaMatrix->rows(), 4);
  ASSERT_EQ(tddftbData.gammaMatrix->cols(), 4);
  ASSERT_TRUE(tddftbData.spinConstants);
  for (int atom_i = 0; atom_i < 4; ++atom_i) {
    for (int atom_j = 0; atom_j < 4; ++atom_j) {
      EXPECT_DOUBLE_EQ((*tddftbData.gammaMatrix)(atom_i, atom_j), (*tddftbData.gammaMatrix)(atom_i, atom_j));
    }
  }
}

TEST_F(ABasisPruningTest, ResultsWithBasisPrunerSameWithZeroPerturbativeThresholdRestricted) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculator);
  dftbMethod->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  dftbMethod->calculate("TDDFTB reference calculation.");

  auto tddftbData = dftbMethod->getTDDFTBData();
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDifferenceVector;
  TimeDependentUtils::generateEnergyDifferenceVector(tddftbData.MOEnergies, tddftbData.occupation, energyDifferenceVector);
  TransitionChargesCalculator tcCalculator(tddftbData.molecularOrbitals, tddftbData.overlapMatrix, tddftbData.AOInfo);

  Eigen::MatrixXd transitionCharges =
      tcCalculator.calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(tddftbData.occupation);

  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);

  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, std::vector<Utils::Excitation>> excitations =
      TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(dftbMethod->getTDDFTBData().molecularOrbitals,
                                                                            dftbMethod->getTDDFTBData().occupation);

  OrderedInput<Utils::Reference::Restricted> input{energyDifferenceVector, excitations, transitionCharges, orderMap};

  BasisPruner<Utils::Reference::Restricted> basisPruner(input, tddftbData.gammaMatrix, tddftbData.spinConstants);

  auto prunedData = basisPruner.prune(EnergyThreshold{1.0}, PerturbativeThreshold{0.0}, Utils::SpinTransition::Singlet);

  Eigen::VectorXd orderedEnergies;
  TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(energyDifferenceVector), orderedEnergies, orderMap,
                                     TimeDependentUtils::Direction::To);

  Eigen::MatrixXd orderedTransitionCharges;
  TimeDependentUtils::transformOrder(transitionCharges, orderedTransitionCharges, orderMap, TimeDependentUtils::Direction::To);

  std::vector<Utils::Excitation> orderedExcitations;
  TimeDependentUtils::transformOrder(
      TimeDependentUtils::flatten(TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(
          tddftbData.molecularOrbitals, tddftbData.occupation)),
      orderedExcitations, orderMap, TimeDependentUtils::Direction::To);

  ASSERT_EQ(prunedData.excitations().size(), orderedExcitations.size());
  for (int i = 0; i < orderedEnergies.size(); ++i) {
    EXPECT_DOUBLE_EQ(prunedData.energyDifferences()(i), orderedEnergies(i));
    EXPECT_EQ(prunedData.excitations()[i].occ, orderedExcitations[i].occ);
    EXPECT_EQ(prunedData.excitations()[i].vir, orderedExcitations[i].vir);
    for (int atom = 0; atom < 4; ++atom)
      EXPECT_DOUBLE_EQ(prunedData.transitionCharges()(i, atom), orderedTransitionCharges(i, atom));
  }

  ASSERT_EQ(tddftbData.gammaMatrix->rows(), 4);
  ASSERT_EQ(tddftbData.gammaMatrix->cols(), 4);
  ASSERT_EQ(tddftbData.spinConstants->size(), 4);
  for (int atom_i = 0; atom_i < 4; ++atom_i) {
    for (int atom_j = 0; atom_j < 4; ++atom_j) {
      EXPECT_DOUBLE_EQ((*tddftbData.gammaMatrix)(atom_i, atom_j), (*tddftbData.gammaMatrix)(atom_i, atom_j));
    }
  }

  auto prunedDataTriplet = basisPruner.prune(EnergyThreshold{1.0}, PerturbativeThreshold{0.0}, Utils::SpinTransition::Triplet);
  ASSERT_EQ(prunedDataTriplet.excitations().size(), orderedExcitations.size());
  for (int i = 0; i < orderedEnergies.size(); ++i) {
    EXPECT_DOUBLE_EQ(prunedDataTriplet.energyDifferences()(i), orderedEnergies(i));
    EXPECT_EQ(prunedDataTriplet.excitations()[i].occ, orderedExcitations[i].occ);
    EXPECT_EQ(prunedDataTriplet.excitations()[i].vir, orderedExcitations[i].vir);
    for (int atom = 0; atom < 4; ++atom)
      EXPECT_DOUBLE_EQ(prunedDataTriplet.transitionCharges()(i, atom), orderedTransitionCharges(i, atom));
  }

  ASSERT_EQ(tddftbData.gammaMatrix->rows(), 4);
  ASSERT_EQ(tddftbData.gammaMatrix->cols(), 4);
  ASSERT_TRUE(tddftbData.spinConstants);
  for (int atom_i = 0; atom_i < 4; ++atom_i) {
    for (int atom_j = 0; atom_j < 4; ++atom_j) {
      EXPECT_DOUBLE_EQ((*tddftbData.gammaMatrix)(atom_i, atom_j), (*tddftbData.gammaMatrix)(atom_i, atom_j));
    }
  }
}

TEST_F(ABasisPruningTest, ResultsWithBasisPrunerSameWithZeroPerturbativeThresholdUnrestricted) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculator);
  dftbMethod->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  dftbMethod->settings().modifyString(Utils::SettingsNames::spinMode,
                                      Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  dftbMethod->calculate("TDDFTB reference calculation.");

  auto tddftbData = dftbMethod->getTDDFTBData();
  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> energyDifferenceVector;
  TimeDependentUtils::generateEnergyDifferenceVector(tddftbData.MOEnergies, tddftbData.occupation, energyDifferenceVector);
  TransitionChargesCalculator tcCalculator(tddftbData.molecularOrbitals, tddftbData.overlapMatrix, tddftbData.AOInfo);

  Eigen::MatrixXd transitionCharges =
      tcCalculator.calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(tddftbData.occupation);

  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(energyDifferenceVector);

  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<Utils::Excitation>> excitations =
      TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(dftbMethod->getTDDFTBData().molecularOrbitals,
                                                                              dftbMethod->getTDDFTBData().occupation);

  OrderedInput<Utils::Reference::Unrestricted> input{energyDifferenceVector, excitations, transitionCharges, orderMap};

  BasisPruner<Utils::Reference::Unrestricted> basisPruner(input, tddftbData.gammaMatrix, tddftbData.spinConstants);

  auto prunedData = basisPruner.prune(EnergyThreshold{1.0}, PerturbativeThreshold{0.0});

  LinearResponseCalculator::GuessSpecifier fakeGuess;
  fakeGuess.unrestricted =
      Eigen::MatrixXd::Random(energyDifferenceVector.alpha.size() + energyDifferenceVector.beta.size(), 3);
  auto prunedGuess = basisPruner.prune(fakeGuess);

  Eigen::MatrixXd orderedGuess;
  TimeDependentUtils::transformOrder(fakeGuess.unrestricted, orderedGuess, orderMap, TimeDependentUtils::Direction::To);

  Eigen::VectorXd orderedEnergies;
  TimeDependentUtils::transformOrder(TimeDependentUtils::flatten(energyDifferenceVector), orderedEnergies, orderMap,
                                     TimeDependentUtils::Direction::To);

  Eigen::MatrixXd orderedTransitionCharges;
  TimeDependentUtils::transformOrder(transitionCharges, orderedTransitionCharges, orderMap, TimeDependentUtils::Direction::To);

  std::vector<Utils::Excitation> orderedExcitations;
  TimeDependentUtils::transformOrder(
      TimeDependentUtils::flatten(TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(
          tddftbData.molecularOrbitals, tddftbData.occupation)),
      orderedExcitations, orderMap, TimeDependentUtils::Direction::To);

  ASSERT_EQ(prunedData.excitations().size(), orderedExcitations.size());
  for (int i = 0; i < orderedEnergies.size(); ++i) {
    EXPECT_DOUBLE_EQ(prunedData.energyDifferences()(i), orderedEnergies(i));
    EXPECT_DOUBLE_EQ(fakeGuess.unrestricted.col(0)(i), prunedGuess->unrestricted.col(0)(i));
    EXPECT_DOUBLE_EQ(fakeGuess.unrestricted.col(1)(i), prunedGuess->unrestricted.col(1)(i));
    EXPECT_DOUBLE_EQ(fakeGuess.unrestricted.col(2)(i), prunedGuess->unrestricted.col(2)(i));
    EXPECT_EQ(prunedData.excitations()[i].occ, orderedExcitations[i].occ);
    EXPECT_EQ(prunedData.excitations()[i].vir, orderedExcitations[i].vir);
    for (int atom = 0; atom < 4; ++atom)
      EXPECT_DOUBLE_EQ(prunedData.transitionCharges()(i, atom), orderedTransitionCharges(i, atom));
  }

  ASSERT_EQ(tddftbData.gammaMatrix->rows(), 4);
  ASSERT_EQ(tddftbData.gammaMatrix->cols(), 4);
  ASSERT_TRUE(tddftbData.spinConstants);
  for (int atom_i = 0; atom_i < 4; ++atom_i) {
    for (int atom_j = 0; atom_j < 4; ++atom_j) {
      EXPECT_DOUBLE_EQ((*tddftbData.gammaMatrix)(atom_i, atom_j), (*tddftbData.gammaMatrix)(atom_i, atom_j));
    }
  }
}

TEST_F(ABasisPruningTest, StandardThresholdsGiveSimilarResultsRestrictedSinglet) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculatorOxa);
  tddftbCalculator.setReferenceCalculator(dftbMethod);
  tddftbCalculator.referenceCalculation();
  // tddftbCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 40);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::energyThreshold, 0.22);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::perturbativeThreshold, 1e-4);
  tddftbCalculator.settings().modifyString(Utils::SettingsNames::pruneBasis, Utils::SettingsNames::PruningOptions::energy);
  auto result2 = tddftbCalculator.calculate();
  Eigen::VectorXd eigenvalues2 = result2.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues;
  Eigen::VectorXd eigenvalues1(42);
  eigenvalues1 << 0.153478, 0.155182, 0.160037, 0.163228, 0.164069, 0.164803, 0.167295, 0.169175, 0.174264, 0.176455,
      0.177058, 0.178246, 0.178561, 0.180074, 0.182639, 0.185729, 0.186709, 0.191444, 0.194681, 0.195288, 0.196206,
      0.199013, 0.201433, 0.20627, 0.207909, 0.208517, 0.20956, 0.209771, 0.210431, 0.212724, 0.2131, 0.214664,
      0.215186, 0.215928, 0.216672, 0.219224, 0.21991, 0.220202, 0.221662, 0.222714, 0.22511, 0.228115;
  double rmsd = Utils::Constants::ev_per_hartree *
                std::sqrt((eigenvalues2 - eigenvalues1).array().square().sum() / eigenvalues2.size());
  ASSERT_LE(rmsd, 0.05);
}

TEST_F(ABasisPruningTest, StandardThresholdsGiveSimilarResultsRestrictedTriplet) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculatorOxa);
  tddftbCalculator.setReferenceCalculator(dftbMethod);
  tddftbCalculator.referenceCalculation();
  // tddftbCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 40);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::energyThreshold, 0.22);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::perturbativeThreshold, 1e-4);
  tddftbCalculator.settings().modifyString(Utils::SettingsNames::pruneBasis, Utils::SettingsNames::PruningOptions::energy);
  tddftbCalculator.settings().modifyString(Utils::SettingsNames::spinBlock, Utils::SettingsNames::SpinBlocks::triplet);
  auto result2 = tddftbCalculator.calculate();
  Eigen::VectorXd eigenvalues2 = result2.get<Utils::Property::ExcitedStates>().triplet->eigenStates.eigenValues;
  Eigen::VectorXd eigenvalues1(42);
  eigenvalues1 << 0.141363, 0.149523, 0.159815, 0.160017, 0.162503, 0.163034, 0.163483, 0.164589, 0.166753, 0.175354,
      0.177039, 0.177993, 0.178213, 0.179226, 0.17927, 0.179779, 0.181208, 0.183094, 0.184982, 0.185995, 0.187676,
      0.188506, 0.189596, 0.190461, 0.191475, 0.195242, 0.197594, 0.198516, 0.199019, 0.208461, 0.209063, 0.209613,
      0.21009, 0.210378, 0.212282, 0.212743, 0.213878, 0.214877, 0.215494, 0.21657, 0.218427, 0.219568;

  double rmsd = Utils::Constants::ev_per_hartree *
                std::sqrt((eigenvalues2 - eigenvalues1).array().square().sum() / eigenvalues2.size());
  ASSERT_LE(rmsd, 0.05);
}
TEST_F(ABasisPruningTest, StandardThresholdsGiveSimilarResultsUnrestricted) {
  auto dftbMethod = std::dynamic_pointer_cast<DFTBMethodWrapper>(calculatorOxa);
  dftbMethod->settings().modifyString(Utils::SettingsNames::spinMode,
                                      Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  tddftbCalculator.setReferenceCalculator(dftbMethod);
  tddftbCalculator.referenceCalculation();
  // tddftbCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 40);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::energyThreshold, 0.22);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::perturbativeThreshold, 1e-4);
  tddftbCalculator.settings().modifyString(Utils::SettingsNames::pruneBasis, Utils::SettingsNames::PruningOptions::energy);
  auto result2 = tddftbCalculator.calculate();
  Eigen::VectorXd eigenvalues2 = result2.get<Utils::Property::ExcitedStates>().unrestricted->eigenStates.eigenValues;
  Eigen::VectorXd eigenvalues1(84);
  eigenvalues1 << 0.141363, 0.149523, 0.153478, 0.155182, 0.159815, 0.160017, 0.160037, 0.162503, 0.163034, 0.163228,
      0.163483, 0.164069, 0.164589, 0.164803, 0.166753, 0.167295, 0.169175, 0.174264, 0.175354, 0.176455, 0.177039,
      0.177058, 0.177993, 0.178213, 0.178246, 0.178561, 0.179226, 0.17927, 0.179779, 0.180074, 0.181208, 0.182639,
      0.183094, 0.184982, 0.185729, 0.185995, 0.186709, 0.187676, 0.188506, 0.189596, 0.190461, 0.191444, 0.191475,
      0.194681, 0.195242, 0.195288, 0.196206, 0.197594, 0.198516, 0.199013, 0.199019, 0.201433, 0.20627, 0.207909,
      0.208461, 0.208517, 0.209063, 0.20956, 0.209613, 0.209771, 0.21009, 0.210378, 0.210431, 0.212282, 0.212724,
      0.212743, 0.2131, 0.213878, 0.214664, 0.214877, 0.215186, 0.215494, 0.215928, 0.21657, 0.216672, 0.218427,
      0.219224, 0.219568, 0.21991, 0.220202, 0.220224, 0.221662, 0.221944, 0.222714;

  double rmsd = Utils::Constants::ev_per_hartree *
                std::sqrt((eigenvalues2 - eigenvalues1).array().square().sum() / eigenvalues2.size());
  ASSERT_LE(rmsd, 0.05);
}

TEST_F(ABasisPruningTest, BinaryTestPrunedRHF) {
  calculatorOxa->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  tddftbCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  tddftbCalculator.settings().modifyString(Utils::SettingsNames::pruneBasis, "energy");
  tddftbCalculator.settings().modifyString("gep_algo", "standard");
  // tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::energyThreshold, 0.22);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::perturbativeThreshold, 1e-4);
  tddftbCalculator.applySettings();
  calculatorOxa->calculate("");
  tddftbCalculator.setReferenceCalculator(calculatorOxa);
  Eigen::VectorXd eigenvalues1(3);
  eigenvalues1 << 0.153478, 0.155182, 0.160037;
  auto results = tddftbCalculator.calculate();
  Eigen::VectorXd eigenvalues2 = results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues;
  double rmsd = Utils::Constants::ev_per_hartree *
                std::sqrt((eigenvalues2 - eigenvalues1).array().square().sum() / eigenvalues2.size());
  ASSERT_LE(rmsd, 0.05);
}
TEST_F(ABasisPruningTest, BinaryTestPrunedUHF) {
  calculatorOxa->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixMO);
  calculatorOxa->settings().modifyString(Utils::SettingsNames::spinMode,
                                         Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  tddftbCalculator.settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 3);
  tddftbCalculator.settings().modifyString(Utils::SettingsNames::pruneBasis, "energy");
  // tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::energyThreshold, 0.22);
  tddftbCalculator.settings().modifyDouble(Utils::SettingsNames::perturbativeThreshold, 1e-4);
  tddftbCalculator.applySettings();
  calculatorOxa->calculate("");
  tddftbCalculator.setReferenceCalculator(calculatorOxa);

  auto results = tddftbCalculator.calculate();
  Eigen::VectorXd eigenvalues2 = results.get<Utils::Property::ExcitedStates>().unrestricted->eigenStates.eigenValues;
  Eigen::VectorXd eigenvalues1(3);
  eigenvalues1 << 0.141363, 0.149523, 0.153478;
  double rmsd = Utils::Constants::ev_per_hartree *
                std::sqrt((eigenvalues2 - eigenvalues1).array().square().sum() / eigenvalues2.size());
  ASSERT_LE(rmsd, 0.05);
}
} // namespace Sparrow
} // namespace Scine
