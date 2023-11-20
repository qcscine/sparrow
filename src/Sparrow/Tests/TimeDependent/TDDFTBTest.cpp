/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Dftb/Dftb0/Wrapper/DFTB0MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/OrderedInput.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBCalculator.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBData.h>
#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBSigmaVectorEvaluator.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/TransitionChargesCalculator.h>
#include <Sparrow/Implementations/TimeDependent/DiagonalPreconditionerEvaluator.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/Math/IterativeDiagonalizer/DavidsonDiagonalizer.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Sparrow {
class ATDDFTBTestCalculation : public Test {
 public:
  std::unique_ptr<TransitionChargesCalculator> transChargeCalc;
  std::shared_ptr<dftb::DFTB2> calculatorEtene;
  std::shared_ptr<DFTB2MethodWrapper> dftbMethod;
  std::shared_ptr<DFTB2MethodWrapper> dftbPorphMethod;
  std::shared_ptr<DFTB0MethodWrapper> dftb0Method;
  std::shared_ptr<TDDFTBCalculator> tddftbCalculator;
  std::shared_ptr<Core::CalculatorWithReference> interfaceCalculator;

  Core::Log log;

 private:
  void SetUp() final {
    std::stringstream etene{"6\n\n"
                            "C          0.94815        0.05810        0.05008\n"
                            "C          2.28393        0.05810        0.05008\n"
                            "H          0.38826       -0.56287        0.74233\n"
                            "H          0.38826        0.67907       -0.64217\n"
                            "H          2.84382       -0.56287        0.74233\n"
                            "H          2.84382        0.67907       -0.64217\n"};
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
    log = Core::Log::silent();
    auto porphyrine = Utils::XyzStreamHandler::read(ssPorphyrine);
    dftbMethod = std::make_shared<DFTB2MethodWrapper>();
    dftbMethod->setLog(Core::Log::silent());
    dftbMethod->settings().modifyDouble(Utils::SettingsNames::densityRmsdCriterion, 1e-9);
    dftbPorphMethod = std::make_shared<DFTB2MethodWrapper>();
    dftbPorphMethod->setLog(Core::Log::silent());
    dftbPorphMethod->settings().modifyDouble(Utils::SettingsNames::densityRmsdCriterion, 1e-15);
    dftbPorphMethod->setStructure(porphyrine);
    dftb0Method = std::make_shared<DFTB0MethodWrapper>();
    dftb0Method->setLog(Core::Log::silent());
    tddftbCalculator = std::make_shared<TDDFTBCalculator>();
    tddftbCalculator->setLog(Core::Log::silent());

    auto& manager = Core::ModuleManager::getInstance();
    interfaceCalculator = manager.get<Core::CalculatorWithReference>("TD-DFTB");
    interfaceCalculator->setLog(Core::Log::silent());
  }
};

/**
 * @brief Test that mathematical operations are coherent.
 * @test CanCalculateTransitionCharges
 * It is very difficult to test the intermediates, so I can at least test that
 * my math is on point by implementing a naive version of the formulae in
 * A. Dominguez, B. Aradi, T. Frauenheim, V. Lutsker, T.A. Niehaus,
 * Extensions of the Time-Dependent Density Functional Based Tight-Binding Approach,
 * 2013.
 *
 */
TEST_F(ATDDFTBTestCalculation, CanCalculateTransitionCharges) {
  Eigen::MatrixXd refTransQ(16, 5);
  refTransQ << -0.53298242221662728, 0.11417760560823831, 5.0202999256975245E-003, -0.15779757137823897,
      -0.50682540581525803, -9.3961145881681196E-002, 2.2128962148245557E-002, 8.4458671606008895E-002,
      -0.58856449960896340, 0.13798084313095874, -8.1635585946912920E-002, -0.22039922685955038, 0.25239430861461465,
      -4.6798155309363126E-003, 4.2762195371146637E-002, -0.62634353135472653, 7.6804730210780395E-002,
      9.2958122398791879E-002, 2.7838312403571853E-002, 0.29831850550636407, 0.36100857097500522, 0.10806390910889868,
      1.1588514120173796E-002, 1.4035057310773324E-002, 4.2117083302872053E-003, 6.8854455783794705E-002,
      9.3808682355889628E-002, 0.26552435748908076, 0.11350302880698007, 3.3958182542152618E-002,
      1.0733089818728723E-002, 8.2789202878895365E-002, 4.2593725640938810E-003, -8.5082057018763362E-003,
      -2.4012110428906070E-002, -9.2521703735310764E-003, 1.8470782192934326E-002, 5.2149840617228486E-002,
      -6.6003367183572664E-002, 0.13148099626342488, 0.37178494263476258, 1.7505496159755874E-002,
      5.5456709060627267E-002, -3.7697365832783547E-002, -0.11043739395445905, -0.31234844635070702,
      -0.26012241396842062, 0.21750310706019091, 0.44970800499409458, -0.18096977536096984, -2.5289585473866647E-002,
      -0.12146374182251074, 4.8874395070582166E-002, 6.8279911327645522E-003, 3.7543238107786581E-003,
      -1.5076820477378421E-003, -2.0936171034587438E-004, -0.24116430904945088, -0.10478835326909049,
      6.4904550634281438E-002, 4.2158354882298389E-002, 5.8870754320221822E-003, -1.8548910602008706E-003,
      5.5677861677582402E-002, 2.2103144476581394E-003, -1.7657746944184108E-002, 1.6443083573503683E-002,
      -9.8050219320831942E-003, 7.8471657576736104E-002, -7.3080594977210783E-002, 2.8531567104374857E-002,
      -0.22846704313246924, 0.21277721035425962, 1.6823513974941532E-002, 3.7158547799486415E-002,
      -7.2332315431028138E-002, -0.29761829834943415, 0.27718300390746847, 0.20848201983874620, 0.27037335973805776;

  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.7    0.72    0.64\n"
                       "H     -0.635   -0.635    0.635\n"
                       "H     -0.71     0.7287000000   -0.7287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  dftb::DFTB2 method;
  method.setAtomCollection(Utils::XyzStreamHandler::read(ss));
  method.initializeFromParameterPath("mio-1-1");
  method.setConvergenceCriteria({{}, 1e-9});
  method.calculate(Utils::Derivative::None, log);

  Eigen::MatrixXd ovlp = method.getOverlapMatrix();
  Eigen::MatrixXd coeffs = method.getMolecularOrbitals().restrictedMatrix();

  double chargeOnC34 = 0.0;
  for (int aoOnC = 0; aoOnC < 4; ++aoOnC) {
    for (int contractedAO = 0; contractedAO < coeffs.rows(); ++contractedAO) {
      chargeOnC34 += 0.5 * (coeffs(aoOnC, 3) * ovlp(aoOnC, contractedAO) * coeffs(contractedAO, 4) +
                            coeffs(aoOnC, 4) * ovlp(aoOnC, contractedAO) * coeffs(contractedAO, 3));
    }
  }

  transChargeCalc = std::make_unique<TransitionChargesCalculator>(method.getMolecularOrbitals(), method.getOverlapMatrix(),
                                                                  method.getAtomsOrbitalsIndexesHolder());
  ASSERT_THAT(transChargeCalc->calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(
                  method.getElectronicOccupation())(3 * 4 + 0, 0),
              DoubleNear(chargeOnC34, 1e-5));
  method.setUnrestrictedCalculation(true);
  method.calculate(Utils::Derivative::None, log);
  transChargeCalc->fillOverlapProductMatrix();
  Eigen::MatrixXd alphaCoeffs = method.getMolecularOrbitals().alphaMatrix();
  Eigen::MatrixXd betaCoeffs = method.getMolecularOrbitals().betaMatrix();
  double chargeOnC34Alpha = 0.0;
  double chargeOnC34Beta = 0.0;
  for (int aoOnC = 0; aoOnC < 4; ++aoOnC) {
    for (int contractedAO = 0; contractedAO < coeffs.rows(); ++contractedAO) {
      chargeOnC34Alpha += 0.5 * (alphaCoeffs(aoOnC, 3) * ovlp(aoOnC, contractedAO) * alphaCoeffs(contractedAO, 4) +
                                 alphaCoeffs(aoOnC, 4) * ovlp(aoOnC, contractedAO) * alphaCoeffs(contractedAO, 3));
      chargeOnC34Beta += 0.5 * (betaCoeffs(aoOnC, 3) * ovlp(aoOnC, contractedAO) * betaCoeffs(contractedAO, 4) +
                                betaCoeffs(aoOnC, 4) * ovlp(aoOnC, contractedAO) * betaCoeffs(contractedAO, 3));
    }
  }
  ASSERT_THAT(transChargeCalc->calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(
                  method.getElectronicOccupation())(3 * 4 + 0, 0),
              DoubleNear(chargeOnC34Alpha, 1e-5));
  ASSERT_THAT(transChargeCalc->calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(
                  method.getElectronicOccupation())(3 * 4 + 16, 0),
              DoubleNear(chargeOnC34Beta, 1e-5));
}

// Check Gamma Matrix against gamma matrix calculated with DFTB+
TEST_F(ATDDFTBTestCalculation, TDDFTBDataRecordsRightConstants) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.7    0.72    0.64\n"
                       "H     -0.635   -0.635    0.635\n"
                       "H     -0.71     0.7287000000   -0.7287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  DFTB2MethodWrapper calc;
  calc.setLog(log);
  calc.settings().modifyString(Utils::SettingsNames::methodParameters, "mio-1-1");
  calc.setStructure(Utils::XyzStreamHandler::read(ss));
  calc.calculate("");
  auto data = calc.getTDDFTBData();
  // Check spin constants against the ones given in the mio-1-1 parameter file.
  ASSERT_EQ(data.spinConstants->size(), 5);
  EXPECT_DOUBLE_EQ((*data.spinConstants)(0), -0.02265);
  EXPECT_DOUBLE_EQ((*data.spinConstants)(1), -0.07174);
  EXPECT_DOUBLE_EQ((*data.spinConstants)(2), -0.07174);
  EXPECT_DOUBLE_EQ((*data.spinConstants)(3), -0.07174);
  EXPECT_DOUBLE_EQ((*data.spinConstants)(4), -0.07174);

  // Check gamma matrix
  ASSERT_EQ(data.gammaMatrix->cols(), 5);
  ASSERT_EQ(data.gammaMatrix->rows(), 5);
  Eigen::MatrixXd refMatrix(5, 5);
  refMatrix << 0.36470000000000002, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
      0.31363790781240219, 0.41949999999999998, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
      0.32240853157413785, 0.25358898147087322, 0.41949999999999998, 0.0000000000000000, 0.0000000000000000,
      0.30775942178169513, 0.24780029015390459, 0.25100708718481635, 0.41949999999999998, 0.0000000000000000,
      0.32345428437310619, 0.25822129786110948, 0.26457379404845255, 0.25294438901675798, 0.41949999999999998;
  refMatrix = refMatrix.selfadjointView<Eigen::Lower>();

  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      EXPECT_THAT((*data.gammaMatrix)(row, col), DoubleNear(refMatrix(row, col), 1e-4));
    }
  }
  calc.settings().modifyString(Utils::SettingsNames::spinMode,
                               Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  calc.calculate("");
  auto data2 = calc.getTDDFTBData();
  // Check spin constants against the ones given in the mio-1-1 parameter file.
  ASSERT_EQ(data2.spinConstants->size(), 5);
  EXPECT_DOUBLE_EQ((*data2.spinConstants)(0), -0.02265);
  EXPECT_DOUBLE_EQ((*data2.spinConstants)(1), -0.07174);
  EXPECT_DOUBLE_EQ((*data2.spinConstants)(2), -0.07174);
  EXPECT_DOUBLE_EQ((*data2.spinConstants)(3), -0.07174);
  EXPECT_DOUBLE_EQ((*data2.spinConstants)(4), -0.07174);

  // Check gamma matrix
  ASSERT_EQ(data2.gammaMatrix->cols(), 5);
  ASSERT_EQ(data2.gammaMatrix->rows(), 5);
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 5; ++col) {
      EXPECT_THAT((*data2.gammaMatrix)(row, col), DoubleNear(refMatrix(row, col), 1e-4));
    }
  }
}

TEST_F(ATDDFTBTestCalculation, CanPerformRestrictedSingletTDDFTCalculation) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.7    0.72    0.64\n"
                       "H     -0.635   -0.635    0.635\n"
                       "H     -0.71     0.7287000000   -0.7287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  DFTB2MethodWrapper calc;
  calc.setLog(log);
  calc.settings().modifyString(Utils::SettingsNames::methodParameters, "mio-1-1");
  calc.setStructure(Utils::XyzStreamHandler::read(ss));
  calc.calculate("");
  auto data = calc.getTDDFTBData();

  TransitionChargesCalculator tcCalc(data.molecularOrbitals, data.overlapMatrix, data.AOInfo);
  Eigen::MatrixXd transitionCharges =
      tcCalc.calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(data.occupation);

  auto excitations =
      TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(data.molecularOrbitals, data.occupation);
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> orbEnergyDifference;
  TimeDependentUtils::generateEnergyDifferenceVector(data.MOEnergies, data.occupation, orbEnergyDifference);
  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(orbEnergyDifference);

  OrderedInput<Utils::Reference::Restricted> input{orbEnergyDifference, excitations, transitionCharges, orderMap};

  auto sve = std::make_unique<TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>>(
      data.gammaMatrix, data.spinConstants, input, Utils::SpinTransition::Singlet);
  Utils::NonOrthogonalDavidson diag(16, 16);
  diag.setSigmaVectorEvaluator(std::move(sve));
  diag.setPreconditionerEvaluator(std::make_unique<DiagonalPreconditionerEvaluator>(input.energyDifferences()));

  auto result = diag.solve(log);

  Eigen::VectorXd energies = Utils::Constants::ev_per_hartree * result.eigenValues.cwiseSqrt();
  EXPECT_THAT(energies(0), DoubleNear(13.355, 1e-3));
  EXPECT_THAT(energies(1), DoubleNear(13.485, 1e-3));
  EXPECT_THAT(energies(2), DoubleNear(14.582, 1e-3));
  EXPECT_THAT(energies(3), DoubleNear(15.581, 1e-3));
  EXPECT_THAT(energies(4), DoubleNear(15.594, 1e-3));
  EXPECT_THAT(energies(5), DoubleNear(16.894, 1e-3));
  EXPECT_THAT(energies(6), DoubleNear(18.039, 1e-3));
  EXPECT_THAT(energies(7), DoubleNear(18.660, 1e-3));
  EXPECT_THAT(energies(8), DoubleNear(19.323, 1e-3));
  EXPECT_THAT(energies(9), DoubleNear(20.513, 1e-3));
}

TEST_F(ATDDFTBTestCalculation, CanPerformRestrictedTripletTDDFTCalculation) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.7    0.72    0.64\n"
                       "H     -0.635   -0.635    0.635\n"
                       "H     -0.71     0.7287000000   -0.7287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  DFTB2MethodWrapper calc;
  calc.setLog(log);
  calc.settings().modifyString(Utils::SettingsNames::methodParameters, "mio-1-1");
  calc.setStructure(Utils::XyzStreamHandler::read(ss));
  calc.calculate("");
  auto data = calc.getTDDFTBData();

  auto excitations =
      TimeDependentUtils::generateExcitations<Utils::Reference::Restricted>(data.molecularOrbitals, data.occupation);
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> orbEnergyDifference;
  TimeDependentUtils::generateEnergyDifferenceVector(data.MOEnergies, data.occupation, orbEnergyDifference);
  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(orbEnergyDifference);

  TransitionChargesCalculator tcCalc(data.molecularOrbitals, data.overlapMatrix, data.AOInfo);
  Eigen::MatrixXd transitionCharges =
      tcCalc.calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(data.occupation);

  OrderedInput<Utils::Reference::Restricted> input{orbEnergyDifference, excitations, transitionCharges, orderMap};

  auto sve = std::make_unique<TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>>(
      data.gammaMatrix, data.spinConstants, input, Utils::SpinTransition::Triplet);
  Utils::NonOrthogonalDavidson diag(16, 16);
  diag.setSigmaVectorEvaluator(std::move(sve));
  diag.setPreconditionerEvaluator(std::make_unique<DiagonalPreconditionerEvaluator>(input.energyDifferences()));

  auto result = diag.solve(log);
  Eigen::VectorXd energies = Utils::Constants::ev_per_hartree * result.eigenValues.cwiseSqrt();
  EXPECT_THAT(energies(0), DoubleNear(11.381, 1e-3));
  EXPECT_THAT(energies(1), DoubleNear(13.168, 1e-3));
  EXPECT_THAT(energies(2), DoubleNear(13.467, 1e-3));
  EXPECT_THAT(energies(3), DoubleNear(14.001, 1e-3));
  EXPECT_THAT(energies(4), DoubleNear(15.117, 1e-3));
  EXPECT_THAT(energies(5), DoubleNear(15.479, 1e-3));
  EXPECT_THAT(energies(6), DoubleNear(17.111, 1e-3));
  EXPECT_THAT(energies(7), DoubleNear(17.876, 1e-3));
  EXPECT_THAT(energies(8), DoubleNear(18.063, 1e-3));
  EXPECT_THAT(energies(9), DoubleNear(19.585, 1e-3));
}

TEST_F(ATDDFTBTestCalculation, CanPerformUnrestrictedTDDFTCalculation) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.7    0.72    0.64\n"
                       "H     -0.635   -0.635    0.635\n"
                       "H     -0.71     0.7287000000   -0.7287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  DFTB2MethodWrapper calc;
  calc.setLog(log);
  calc.settings().modifyString(Utils::SettingsNames::methodParameters, "mio-1-1");
  calc.settings().modifyString(Utils::SettingsNames::spinMode,
                               Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  calc.setStructure(Utils::XyzStreamHandler::read(ss));
  calc.calculate("");
  auto data = calc.getTDDFTBData();

  auto excitations =
      TimeDependentUtils::generateExcitations<Utils::Reference::Unrestricted>(data.molecularOrbitals, data.occupation);
  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> orbEnergyDifference;
  TimeDependentUtils::generateEnergyDifferenceVector(data.MOEnergies, data.occupation, orbEnergyDifference);
  auto orderMap = TimeDependentUtils::generateEnergyOrderMap(orbEnergyDifference);

  TransitionChargesCalculator tcCalc(data.molecularOrbitals, data.overlapMatrix, data.AOInfo);
  Eigen::MatrixXd transitionCharges =
      tcCalc.calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(data.occupation);

  OrderedInput<Utils::Reference::Unrestricted> input{orbEnergyDifference, excitations, transitionCharges, orderMap};

  auto sve = std::make_unique<TDDFTBSigmaVectorEvaluator<Utils::Reference::Unrestricted>>(data.gammaMatrix,
                                                                                          data.spinConstants, input);
  Utils::NonOrthogonalDavidson diag(32, 32);
  diag.setSigmaVectorEvaluator(std::move(sve));
  diag.setPreconditionerEvaluator(std::make_unique<DiagonalPreconditionerEvaluator>(input.energyDifferences()));

  auto result = diag.solve(log);
  Eigen::VectorXd energies = Utils::Constants::ev_per_hartree * result.eigenValues.cwiseSqrt();
  EXPECT_THAT(energies(0), DoubleNear(11.381, 1e-3));
  EXPECT_THAT(energies(1), DoubleNear(13.168, 1e-3));
  EXPECT_THAT(energies(2), DoubleNear(13.355, 1e-3));
  EXPECT_THAT(energies(3), DoubleNear(13.467, 1e-3));
  EXPECT_THAT(energies(4), DoubleNear(13.485, 1e-3));
  EXPECT_THAT(energies(5), DoubleNear(14.001, 1e-3));
  EXPECT_THAT(energies(6), DoubleNear(14.582, 1e-3));
  EXPECT_THAT(energies(7), DoubleNear(15.117, 1e-3));
  EXPECT_THAT(energies(8), DoubleNear(15.479, 1e-3));
  EXPECT_THAT(energies(9), DoubleNear(15.581, 1e-3));
  EXPECT_THAT(energies(10), DoubleNear(15.594, 1e-3));
  EXPECT_THAT(energies(11), DoubleNear(16.894, 1e-3));
  EXPECT_THAT(energies(12), DoubleNear(17.111, 1e-3));
  EXPECT_THAT(energies(13), DoubleNear(17.876, 1e-3));
  EXPECT_THAT(energies(14), DoubleNear(18.039, 1e-3));
  EXPECT_THAT(energies(15), DoubleNear(18.063, 1e-3));
  EXPECT_THAT(energies(16), DoubleNear(18.660, 1e-3));
  EXPECT_THAT(energies(17), DoubleNear(19.323, 1e-3));
  EXPECT_THAT(energies(18), DoubleNear(19.585, 1e-3));
  EXPECT_THAT(energies(19), DoubleNear(20.513, 1e-3));
}

TEST_F(ATDDFTBTestCalculation, CanCalculateThroughTDDFTBCalculatorWithdOrbitals) {
  std::stringstream ZnOss("2\n\n"
                          "Zn   0.00000    0.00000    0.00000\n"
                          "O    1.20000    0.00000    0.00000");
  auto dftb2Method = std::make_shared<DFTB2MethodWrapper>();
  dftb2Method->setLog(log);
  dftb2Method->settings().modifyString("method_parameters", "3ob-3-1");
  dftb2Method->setStructure(Utils::XyzStreamHandler::read(ZnOss));

  // NOTE: I had to request 15 roots in DFTB+ because it does not find the lowest roots as
  //       easily as Sparrow.
  Eigen::VectorXd referenceEnergies(10);
  referenceEnergies << 4.2645809988E+00, 4.2645809988E+00, 5.0828471919E+00, 9.5066345607E+00, 9.5066345607E+00,
      9.6859336276E+00, 9.6859336276E+00, 9.6859336276E+00, 1.0202049455E+01, 1.0202049455E+01;

  tddftbCalculator->setReferenceCalculator(dftb2Method);

  tddftbCalculator->referenceCalculation();
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  auto results = tddftbCalculator->calculate();
  Eigen::VectorXd actualEnergies =
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues * Utils::Constants::ev_per_hartree;

  for (int i = 0; i < 10; ++i)
    EXPECT_THAT(actualEnergies(i), DoubleNear(referenceEnergies(i), 1e-4));
}

TEST_F(ATDDFTBTestCalculation, CanCalculateEtheneThroughTDDFTBCalculator) {
  std::stringstream etene{"6\n\n"
                          "C          0.94815        0.05810        0.05008\n"
                          "C          2.28393        0.05810        0.05008\n"
                          "H          0.38826       -0.56287        0.74233\n"
                          "H          0.38826        0.67907       -0.64217\n"
                          "H          2.84382       -0.56287        0.74233\n"
                          "H          2.84382        0.67907       -0.64217\n"};
  auto structure = Utils::XyzStreamHandler::read(etene);
  dftbPorphMethod = std::make_shared<DFTB2MethodWrapper>();
  dftbPorphMethod->setLog(log);
  dftbPorphMethod->settings().modifyDouble(Utils::SettingsNames::densityRmsdCriterion, 1e-9);
  dftbPorphMethod->setStructure(structure);

  tddftbCalculator->setReferenceCalculator(dftbPorphMethod);
  tddftbCalculator->referenceCalculation();
  tddftbCalculator->settings().modifyString(Utils::SettingsNames::spinBlock, Utils::SettingsNames::SpinBlocks::singlet);
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  tddftbCalculator->settings().modifyDouble("convergence", 1e-5);
  auto results = tddftbCalculator->calculate();
  Eigen::VectorXd sortedEigenVector =
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenVectors.col(1).array().abs();
  std::sort(sortedEigenVector.data(), sortedEigenVector.data() + sortedEigenVector.size(), std::greater<>());
  Eigen::VectorXd actualEigenvalues = results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues;

  Eigen::VectorXd referenceEigenvalues(10);
  referenceEigenvalues << 6.4372123108708559E-002, 7.9594182113271461E-002, 9.4467270990339705E-002,
      0.11432557148235396, 0.19206918916273796, 0.31418674501003346, 0.36066681266480594, 0.37129421370841414,
      0.42560368108501884, 0.47506834940552584;
  Eigen::VectorXd referenceEigenvector(36);
  referenceEigenvector << 0.99991195, 0.00796583, 0.00780009, 0.00714281, 0.00079503, 0.00033460, 0.00016785, 0.00004490,
      0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
      0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000;
  Eigen::VectorXd referenceOscStrength(10);
  referenceOscStrength << 0.0, 0.27779743, 0.0, 0.0, 0.0, 0.0, 0.0, 0.50585082, 0.12344951, 0.0;

  for (int i = 0; i < 36; ++i)
    EXPECT_THAT(sortedEigenVector(i), DoubleNear(referenceEigenvector(i), 1e-5));

  Eigen::VectorXd oscStrength = Utils::TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(
      results.get<Utils::Property::ExcitedStates>().singlet->transitionDipoles,
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues);

  for (int i = 0; i < 10; ++i) {
    EXPECT_THAT(oscStrength(i), DoubleNear(referenceOscStrength(i), 1e-5));
    EXPECT_THAT(actualEigenvalues(i), DoubleNear(std::sqrt(referenceEigenvalues(i)), 1e-6));
  }
}
TEST_F(ATDDFTBTestCalculation, CanCalculatePorphyrineThroughTDDFTBCalculator) {
  tddftbCalculator->setReferenceCalculator(dftbPorphMethod);
  tddftbCalculator->referenceCalculation();
  tddftbCalculator->settings().modifyString(Utils::SettingsNames::spinBlock, Utils::SettingsNames::SpinBlocks::singlet);
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  auto results = tddftbCalculator->calculate();

  Eigen::VectorXd referenceEnergies(10);
  Eigen::VectorXd referenceOscStrength(10);
  Eigen::VectorXd referenceSPTransitions(56);
  referenceSPTransitions << 0.99998738, 0.00325800, 0.00325648, 0.00116549, 0.00088365, 0.00061298, 0.00061266,
      0.00058739, 0.00058685, 0.00034416, 0.00028714, 0.00028671, 0.00019486, 0.00013887, 0.00013878, 0.00012238,
      0.00012197, 0.00008484, 0.00007650, 0.00006958, 0.00006920, 0.00005983, 0.00005983, 0.00005338, 0.00005338,
      0.00005294, 0.00004627, 0.00004626, 0.00003971, 0.00003966, 0.00003638, 0.00003278, 0.00003260, 0.00003071,
      0.00002585, 0.00002345, 0.00002338, 0.00002072, 0.00002072, 0.00001956, 0.00001556, 0.00001522, 0.00001515,
      0.00001498, 0.00001350, 0.00001279, 0.00001221, 0.00001220, 0.00001209, 0.00001130, 0.00001108, 0.00001097,
      0.00001095, 0.00001019, 0.00001016, 0.00001008;
  referenceEnergies << 1.2151578473E+00, 1.5922799345E+00, 1.9631328412E+00, 1.9635740483E+00, 2.0982073662E+00,
      2.0983963952E+00, 2.1541952875E+00, 2.4066132528E+00, 2.4185088647E+00, 2.6419966041E+00;
  referenceOscStrength << 0.00000000, 0.00000005, 0.00000000, 0.00000000, 0.00469321, 0.00468552, 0.00000000,
      0.00000000, 0.00000000, 0.01421826;

  Eigen::VectorXd actualEigenvalues =
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues * Utils::Constants::ev_per_hartree;
  Eigen::VectorXd actualOscStrength = Utils::TransitionDipoleCalculator::transitionDipoleMomentToOscillatorStrength(
      results.get<Utils::Property::ExcitedStates>().singlet->transitionDipoles,
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues);
  Eigen::VectorXd firstSortedEigenVector =
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenVectors.col(0);
  std::sort(firstSortedEigenVector.data(), firstSortedEigenVector.data() + firstSortedEigenVector.size(),
            [](auto a, auto b) { return std::abs(a) > std::abs(b); });

  for (int i = 0; i < 10; ++i) {
    EXPECT_THAT(actualEigenvalues(i), DoubleNear(referenceEnergies(i), 1e-3));
    EXPECT_THAT(actualOscStrength(i), DoubleNear(referenceOscStrength(i), 1e-4));
  }
  for (int i = 0; i < 56; ++i)
    EXPECT_THAT(std::abs(firstSortedEigenVector(i)), DoubleNear(referenceSPTransitions(i), 1e-5));
}

// Check that TDA is similar to TDDFTB, I did not think of a better way to test this...
TEST_F(ATDDFTBTestCalculation, CanCalculateEteneThroughTDACalculator) {
  std::stringstream etene{"6\n\n"
                          "C          0.94815        0.05810        0.05008\n"
                          "C          2.28393        0.05810        0.05008\n"
                          "H          0.38826       -0.56287        0.74233\n"
                          "H          0.38826        0.67907       -0.64217\n"
                          "H          2.84382       -0.56287        0.74233\n"
                          "H          2.84382        0.67907       -0.64217\n"};
  auto structure = Utils::XyzStreamHandler::read(etene);
  dftbPorphMethod = std::make_shared<DFTB2MethodWrapper>();
  dftbPorphMethod->setLog(log);
  dftbPorphMethod->settings().modifyDouble(Utils::SettingsNames::densityRmsdCriterion, 1e-9);
  dftbPorphMethod->setStructure(structure);

  tddftbCalculator->setReferenceCalculator(dftbPorphMethod);
  tddftbCalculator->referenceCalculation();
  tddftbCalculator->settings().modifyString(Utils::SettingsNames::spinBlock, Utils::SettingsNames::SpinBlocks::singlet);
  // tddftbCalculator->settings().modifyString("gep_algo", "standard");
  tddftbCalculator->settings().modifyString("gep_algo", "simultaneous_diag");
  tddftbCalculator->settings().modifyDouble("convergence", 1e-8);
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 5);
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::initialSubspaceDimension, 10);
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::maxDavidsonIterations, 1000);
  tddftbCalculator->settings().modifyBool("tda", true);

  // Energies squared
  Eigen::VectorXd referenceTDDFTBEigenvalues(10);
  referenceTDDFTBEigenvalues << 6.4372123108708559E-002, 7.9594182113271461E-002, 9.4467270990339705E-002,
      0.11432557148235396, 0.19206918916273796, 0.31418674501003346, 0.36066681266480594, 0.37129421370841414,
      0.42560368108501884, 0.47506834940552584;
  auto results = tddftbCalculator->calculate();

  Eigen::VectorXd actualEigenvalues =
      results.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues * Utils::Constants::ev_per_hartree;

  Eigen::VectorXd errors =
      referenceTDDFTBEigenvalues.head(5).array().sqrt() * Utils::Constants::ev_per_hartree - actualEigenvalues.array();
  double rmsd = errors.norm() / std::sqrt(1.0 * errors.size());
  for (int i = 0; i < 5; ++i) {
    EXPECT_LE(rmsd, 5e-1);
  }
}

// Test that DFTB0 returns the difference in energy
TEST_F(ATDDFTBTestCalculation, CanCalculateThroughTDDFTB0Calculator) {
  std::stringstream ZnOss("2\n\n"
                          "Zn   0.00000    0.00000    0.00000\n"
                          "O    1.20000    0.00000    0.00000");
  auto dftb0Method = std::make_shared<DFTB0MethodWrapper>();
  dftb0Method->setLog(log);
  dftb0Method->settings().modifyString("method_parameters", "3ob-3-1");
  dftb0Method->setStructure(Utils::XyzStreamHandler::read(ZnOss));
  dftb0Method->setRequiredProperties(Utils::Property::OrbitalEnergies | Utils::Property::ElectronicOccupation);
  auto results = dftb0Method->calculate("");

  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> energyDifferences;
  TimeDependentUtils::generateEnergyDifferenceVector(results.get<Utils::Property::OrbitalEnergies>(),
                                                     results.get<Utils::Property::ElectronicOccupation>(), energyDifferences);

  std::sort(energyDifferences.restricted.data(), energyDifferences.restricted.data() + energyDifferences.restricted.size());

  tddftbCalculator->setReferenceCalculator(dftb0Method);
  tddftbCalculator->referenceCalculation();
  tddftbCalculator->settings().modifyInt(Utils::SettingsNames::numberOfEigenstates, 10);
  auto excitedStateResults = tddftbCalculator->calculate();

  Eigen::VectorXd expectedResults = energyDifferences.restricted.head(10) * Utils::Constants::ev_per_hartree;
  Eigen::VectorXd actualResults = excitedStateResults.get<Utils::Property::ExcitedStates>().singlet->eigenStates.eigenValues *
                                  Utils::Constants::ev_per_hartree;

  for (int i = 0; i < 10; ++i)
    EXPECT_THAT(actualResults(i), DoubleNear(expectedResults(i), 1e-6));
}
} // namespace Sparrow
} // namespace Scine
