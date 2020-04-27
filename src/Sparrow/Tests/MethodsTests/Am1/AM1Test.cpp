/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "../parameters_location.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Nddo/Am1/AM1Method.h>
#include <Sparrow/Implementations/Nddo/Am1/Wrapper/AM1TypeMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/OneElectronMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/ConvergenceAccelerators/ConvergenceAcceleratorFactory.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <Eigen/Core>
#include <boost/dll/runtime_symbol_info.hpp>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace Utils;
using Utils::derivativeType;
using namespace Utils::Constants;

class AAM1Calculation : public Test {
 public:
  AM1Method method;
  std::shared_ptr<AM1MethodWrapper> am1MethodWrapper;
  std::shared_ptr<Core::Calculator> polymorphicMethodWrapper;
  ElementParameters elementParameters;

  void SetUp() override {
    am1MethodWrapper = std::make_shared<AM1MethodWrapper>();
    polymorphicMethodWrapper = std::make_shared<AM1MethodWrapper>();

    method.setMaxIterations(10000);
    method.setConvergenceCriteria(1e-8);
    am1MethodWrapper->settings().modifyInt(Utils::SettingsNames::maxIterations, 10000);
    am1MethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
    am1MethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
    am1MethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);
    polymorphicMethodWrapper->settings().modifyInt(Utils::SettingsNames::maxIterations, 10000);
    polymorphicMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
    polymorphicMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
    polymorphicMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);
    Utils::Log::startConsoleLogging(Utils::SettingsNames::LogLevels::none);
  }
};

TEST_F(AAM1Calculation, HasTheCorrectNumberOfAOsAfterInitialization) {
  std::stringstream ssH("4\n\n"
                        "C      0.0000000000    0.0000000000    0.0000000000\n"
                        "P      0.0529177211   -0.3175063264    0.2645886053\n"
                        "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                        "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ssH);
  method.setStructure(as, parameters_am1);
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(10));
}

TEST_F(AAM1Calculation, GetsCorrectOneElectronPartOfFockForH) {
  std::stringstream ssH("1\n\n"
                        "H -1.0 0.2 -0.3\n");
  auto as = Utils::XyzStreamHandler::read(ssH);
  method.setStructure(as, parameters_am1);

  method.calculateDensityIndependentQuantities();
  ASSERT_THAT(method.getOneElectronMatrix().getMatrix()(0, 0) * ev_per_hartree, DoubleNear(-11.396427, 1e-4));
}

TEST_F(AAM1Calculation, GetsCorrectOneElectronPartOfFockForC) {
  std::stringstream ssC("1\n\n"
                        "C 0.0 0.0 0.0\n");
  auto as = Utils::XyzStreamHandler::read(ssC);
  method.setStructure(as, parameters_am1);
  method.setScfMixer(scf_mixer_t::none);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 4, 4> expected;
  expected << -52.028658, 0, 0, 0, 0, -39.614239, 0, 0, 0, 0, -39.614239, 0, 0, 0, 0, -39.614239;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-4));
}

TEST_F(AAM1Calculation, GetsCorrectOneElectronPartOfFockForH2) {
  std::stringstream ssH2("2\n\n"
                         "H -0.529177  0.105835 -0.158753\n"
                         "H -0.105835  0.105835 -0.158753\n");
  auto as = Utils::XyzStreamHandler::read(ssH2);
  method.setStructure(as, parameters_am1);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 2, 2> expected;
  expected << -23.415584, 0, -5.373535, -23.415584;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-4));
}

TEST_F(AAM1Calculation, GetsCorrectOneElectronPartOfFockForCC) {
  std::stringstream ssC2("2\n\n"
                         "C 0.0 0.0 0.0\n"
                         "C 0.423342  0.0529177 -0.0529177\n");
  auto as = Utils::XyzStreamHandler::read(ssC2);
  method.setStructure(as, parameters_am1);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 8, 8> expected;
  expected << -97.981262, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.550343, -81.933451, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -0.818792, -0.151386, -80.741285, 0.0, 0.0, 0.0, 0.0, 0.0, 0.818792, 0.151386, 0.018923, -80.741285, 0.0, 0.0,
      0.0, 0.0, -14.047246, -4.527135, -0.565891, 0.565891, -97.981262, 0.0, 0.0, 0.0, 4.527135, -4.322779, 0.267067,
      -0.267067, 6.550343, -81.933451, 0.0, 0.0, 0.565891, 0.267067, -6.425929, -0.033383, 0.818792, -0.151386,
      -80.741285, 0.0, -0.565891, -0.267067, -0.033383, -6.425929, -0.818792, 0.151386, 0.018923, -80.741285;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-3));
}

TEST_F(AAM1Calculation, GetsCorrectOneElectronPartOfFockForCH2O) {
  std::stringstream ss("4\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "O      0.0529177211   -0.3175063264    0.2645886053\n"
                       "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                       "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 10, 10> expected;
  expected << -151.724973, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.780226, -127.627565, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 9.081853, 0.367070, -128.660502, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.978496, -0.382423, 1.323895,
      -128.250115, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.096410, -1.279126, 7.674755, -6.395629, -171.464122, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.822641, -12.975084, -0.585656, 0.488047, 2.739877, -147.686804, 0.0, 0.0, 0.0, 0.0, -4.935844,
      -0.585656, -9.558758, -2.928279, -7.579829, 0.468558, -148.529641, 0.0, 0.0, 0.0, 4.113204, 0.488047, -2.928279,
      -10.632461, 6.713787, -0.446867, 1.096667, -148.291036, 0.0, 0.0, -8.838985, 3.273323, -0.654665, 0.981997,
      -7.546637, 4.580123, -3.330998, 3.330998, -133.697936, 0.0, -10.417368, 0.871967, -0.871967, 1.307950, -10.283813,
      1.854250, -4.944667, 4.944667, -5.373535, -144.684916;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-4));
}

TEST_F(AAM1Calculation, GetsCorrectEigenvaluesForC) {
  std::stringstream ss("1\n\n"
                       "C     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);
  method.setScfMixer(scf_mixer_t::none);

  method.convergedCalculation();

  auto eigenvalues = method.getSingleParticleEnergies().getRestrictedEnergies();

  std::vector<double> expected = {-19.28866, -8.02424, -0.04424, -0.04424};
  for (unsigned i = 0; i < eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i] * ev_per_hartree, DoubleNear(expected[i], 1e-4));
  }
}

TEST_F(AAM1Calculation, GetsCorrectUnrestrictedEigenvaluesForLithium) {
  std::stringstream ss("1\n\n"
                       "Li -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);

  method.setUnrestrictedCalculation(true);
  method.setSpinMultiplicity(2);
  method.convergedCalculation(derivativeType::none);

  auto alphaEigenvalues = method.getSingleParticleEnergies().getAlphaEnergies();
  auto betaEigenvalues = method.getSingleParticleEnergies().getBetaEnergies();

  // From Mopac
  double mopacEnergy = -4.93844; // eV
  std::vector<double> alphaExpected = {-4.93844, 5.45631, 5.45631, 5.45631};
  std::vector<double> betaExpected = {0.46149, 5.85628, 5.85628, 5.85628};

  ASSERT_THAT(method.getEnergy() * ev_per_hartree, DoubleNear(mopacEnergy, 1e-5));
  for (unsigned i = 0; i < alphaEigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalues " + std::to_string(i) + ":");
    EXPECT_THAT(alphaEigenvalues[i] * ev_per_hartree, DoubleNear(alphaExpected[i], 1e-5));
    EXPECT_THAT(betaEigenvalues[i] * ev_per_hartree, DoubleNear(betaExpected[i], 1e-5));
  }
}

TEST_F(AAM1Calculation, IsIndependentOfOrderOfAtoms) {
  std::stringstream ss("3\n\n"
                       "O      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.4005871485    0.3100978455    0.0000000000\n"
                       "H     -0.4005871485    0.3100978455    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();
  double energy1 = method.getElectronicEnergy();

  std::stringstream ss2("3\n\n"
                        "H      0.4005871485    0.3100978455    0.0000000000\n"
                        "H     -0.4005871485    0.3100978455    0.0000000000\n"
                        "O      0.0000000000    0.0000000000    0.0000000000\n");
  auto as2 = Utils::XyzStreamHandler::read(ss2);
  method.setStructure(as2, parameters_am1);

  method.convergedCalculation();
  double energy2 = method.getElectronicEnergy();

  ASSERT_THAT(energy1, DoubleNear(energy2, 1e-10));
}

TEST_F(AAM1Calculation, GetsCorrectTotalEnergyForH2) {
  std::stringstream ssH2("2\n\n"
                         "H -0.529177  0.105835 -0.158753\n"
                         "H -0.105835  0.105835 -0.158753\n");
  auto as = Utils::XyzStreamHandler::read(ssH2);
  method.setStructure(as, parameters_am1);
  method.setScfMixer(scf_mixer_t::none);
  method.convergedCalculation();

  auto eigenvalues = method.getSingleParticleEnergies().getRestrictedEnergies();
  std::vector<double> expected = {-16.35554, 6.41069};

  for (unsigned i = 0; i < eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i] * ev_per_hartree, DoubleNear(expected[i], 1e-4));
  }
  ASSERT_THAT(method.getElectronicEnergy() * ev_per_hartree, DoubleNear(-45.14466, 1e-4));
  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(19.14262, 1e-4));
}

TEST_F(AAM1Calculation, GetsCorrectTotalEnergyForH2AtOtherDistance) {
  std::stringstream ss("2\n\n"
                       "H     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();

  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(6.31146, 1e-4));
  ASSERT_THAT(method.getElectronicEnergy() * ev_per_hartree, DoubleNear(-27.48149, 1e-3));
}

TEST_F(AAM1Calculation, GetsCorrectTotalEnergyForC2) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "C     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();

  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(100.68541, 1e-4));
  ASSERT_THAT(method.getElectronicEnergy() * ev_per_hartree, DoubleNear(-341.28101, 1e-3));
}

TEST_F(AAM1Calculation, GetsCorrectRepulsionEnergyForCH) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);

  method.setStructure(as, parameters_am1);
  method.calculate(derivativeType::none);

  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(25.20936, 1e-4));
}

TEST_F(AAM1Calculation, GetsCorrectGradientsForMethane) {
  /* Structure optimized with MOPAC 2016
   *
   *
   * AM1 RHF 1SCF PRECISE GRADIENTS
   *
   *
   * C      0.0000000000    0.0000000000    0.0000000000
   * H      0.6287000000    0.6287000000    0.6287000000
   * H     -0.6287000000   -0.6287000000    0.6287000000
   * H     -0.6287000000    0.6287000000   -0.6287000000
   * H      0.6287000000   -0.6287000000   -0.6287000000
   *
   */

  auto stream = std::stringstream("5\n\n"
                                  "C      0.0000000000    0.0000000000    0.0000000000\n"
                                  "H      0.6287000000    0.6287000000    0.6287000000\n"
                                  "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                                  "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                                  "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  auto structure = Utils::XyzStreamHandler::read(stream);
  method.setStructure(structure, parameters_am1);
  method.convergedCalculation(Utils::derivativeType::first);

  Eigen::RowVector3d cOneGradient(0.000020, 0.000012, 0.000003);
  Eigen::RowVector3d hOneGradient(-12.214666, -12.208263, -12.210581);
  Eigen::RowVector3d hTwoGradient(12.188451, 12.186204, -12.169867);
  Eigen::RowVector3d hThreeGradient(12.186784, -12.169577, 12.188152);
  Eigen::RowVector3d hFourGradient(-12.160589, 12.191624, 12.192292);
  cOneGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hOneGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hTwoGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hThreeGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hFourGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;

  ASSERT_THAT(method.getEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-183.18975, 1e-3));
  ASSERT_TRUE((method.getGradients().row(0) - cOneGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(1) - hOneGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(2) - hTwoGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(3) - hThreeGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(4) - hFourGradient).norm() < 1e-4);
}

TEST_F(AAM1Calculation, GetsZeroGradientForOptimizedMethane) {
  /* Structure optimized with MOPAC 2016
   *
   *
   * AM1 RHF OPT PRECISE GRADIENTS
   *
   *
   * C      0.0000000000    0.0000000000    0.0000000000
   * H      0.6287000000    0.6287000000    0.6287000000
   * H     -0.6287000000   -0.6287000000    0.6287000000
   * H     -0.6287000000    0.6287000000   -0.6287000000
   * H      0.6287000000   -0.6287000000   -0.6287000000
   *
   */
  std::stringstream stream("5\n\n"
                           "C     0.00002085  -0.00000363   -0.00003011\n"
                           "H     0.64131546   0.64266845    0.64155219\n"
                           "H    -0.64127341  -0.64260156    0.64141129\n"
                           "H    -0.64128349   0.64253963   -0.64143849\n"
                           "H     0.64116947  -0.64259398   -0.64142104\n");

  auto structure = Utils::XyzStreamHandler::read(stream);
  method.setStructure(structure, parameters_am1);
  method.convergedCalculation(Utils::derivativeType::first);

  ASSERT_THAT(method.getEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-183.23058, 1e-3));
  ASSERT_TRUE(method.getGradients().row(0).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(1).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(2).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(3).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(4).norm() < 1e-3);
}

TEST_F(AAM1Calculation, GetsCorrectTotalEnergyForEthanol) {
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
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();

  ASSERT_THAT(method.getEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-659.71181, 2e-2));
  ASSERT_THAT(method.getRepulsionEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(1128.16522, 2e-2));
  ASSERT_THAT(method.getElectronicEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-1787.87703, 2e-2));
}

TEST_F(AAM1Calculation, GetsSameTotalEnergyWithAndWithoutWrapper) {
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
  auto as = Utils::XyzStreamHandler::read(ss);
  am1MethodWrapper->setStructure(as);
  auto result = am1MethodWrapper->calculate("");
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();

  ASSERT_THAT(method.getEnergy(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-5));
}

TEST_F(AAM1Calculation, GetsSameTotalEnergyWithPolymorphicMethodAndWithMethod) {
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
  auto as = Utils::XyzStreamHandler::read(ss);
  polymorphicMethodWrapper->setStructure(as);
  auto result = polymorphicMethodWrapper->calculate("");
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();

  ASSERT_THAT(method.getEnergy(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-5));
}

TEST_F(AAM1Calculation, GetsSameTotalEnergyWithDynamicallyLoadedMethod) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyInt(Utils::SettingsNames::maxIterations, 10000);
  dynamicallyLoadedMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);
  auto result = dynamicallyLoadedMethodWrapper->calculate("");
  method.setStructure(as, parameters_am1);
  method.convergedCalculation();
  ASSERT_THAT(result.get<Utils::Property::Energy>() * Utils::Constants::ev_per_hartree, DoubleNear(-659.71181, 2e-2));
  ASSERT_THAT(method.getEnergy(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-5));
}

TEST_F(AAM1Calculation, MethodWrapperCanBeCloned) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_pm6);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);

  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  auto scc = Utils::SettingsNames::selfConsistanceCriterion;
  ASSERT_EQ(cloned->settings().getDouble(scc), dynamicallyLoadedMethodWrapper->settings().getDouble(scc));
}

TEST_F(AAM1Calculation, StructureIsCorrectlyCloned) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);

  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  ASSERT_EQ(cloned->getPositions(), dynamicallyLoadedMethodWrapper->getPositions());
}

TEST_F(AAM1Calculation, ClonedMethodCanCalculate) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);

  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  auto resultCloned = cloned->calculate("");
  auto result = dynamicallyLoadedMethodWrapper->calculate("");
  ASSERT_THAT(resultCloned.get<Utils::Property::Energy>(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-9));
}

TEST_F(AAM1Calculation, ClonedMethodCanCalculateGradients) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);

  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  dynamicallyLoadedMethodWrapper->setRequiredProperties(Utils::Property::Gradients);
  cloned->setRequiredProperties(Utils::Property::Gradients);

  auto resultCloned = cloned->calculate("");
  auto result = dynamicallyLoadedMethodWrapper->calculate("");

  for (int atom = 0; atom < cloned->getPositions().rows(); ++atom) {
    for (int dimension = 0; dimension < 3; ++dimension) {
      ASSERT_THAT(resultCloned.get<Utils::Property::Gradients>().row(atom)(dimension),
                  DoubleNear(result.get<Utils::Property::Gradients>().row(atom)(dimension), 1e-7));
    }
  }
}

TEST_F(AAM1Calculation, ClonedMethodCanCalculateGradientsWithDifferentNumberCores) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, parameters_root);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);

  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  dynamicallyLoadedMethodWrapper->setRequiredProperties(Utils::Property::Gradients);
  cloned->setRequiredProperties(Utils::Property::Gradients);

  auto numThreads = omp_get_num_threads();
  omp_set_num_threads(1);
  auto resultCloned = cloned->calculate("");
  omp_set_num_threads(4);
  auto result = dynamicallyLoadedMethodWrapper->calculate("");
  omp_set_num_threads(numThreads);

  for (int i = 0; i < resultCloned.get<Utils::Property::Gradients>().size(); ++i) {
    ASSERT_THAT(resultCloned.get<Utils::Property::Gradients>()(i),
                DoubleNear(result.get<Utils::Property::Gradients>()(i), 1e-7));
  }
  ASSERT_EQ(resultCloned.get<Utils::Property::SuccessfulCalculation>(), true);
}

TEST_F(AAM1Calculation, ClonedMethodCopiesResultsCorrectly) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("AM1");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_am1);

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
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);

  auto result = dynamicallyLoadedMethodWrapper->calculate();
  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  ASSERT_THAT(cloned->results().get<Utils::Property::Energy>(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-9));
}

TEST_F(AAM1Calculation, AtomCollectionCanBeReturned) {
  std::stringstream ssH("4\n\n"
                        "C      0.0000000000    0.0000000000    0.0000000000\n"
                        "C      0.0529177211   -0.3175063264    0.2645886053\n"
                        "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                        "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto structure = Utils::XyzStreamHandler::read(ssH);
  am1MethodWrapper->setStructure(structure);
  ASSERT_EQ(structure.getPositions(), am1MethodWrapper->getStructure()->getPositions());
  for (int i = 0; i < structure.getElements().size(); ++i)
    ASSERT_EQ(structure.getElements()[i], am1MethodWrapper->getStructure()->getElements()[i]);
}
} // namespace Sparrow
} // namespace Scine
