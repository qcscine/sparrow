/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "../parameters_location.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Nddo/Mndo/MNDOMethod.h>
#include <Sparrow/Implementations/Nddo/Mndo/Wrapper/MNDOMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/OneElectronMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/ConvergenceAccelerators/ConvergenceAcceleratorFactory.h>
#include <gmock/gmock.h>
#include <Eigen/Core>
#include <boost/dll/runtime_symbol_info.hpp>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace Utils::AutomaticDifferentiation;
using namespace Utils::Constants;

class AMNDOCalculation : public Test {
 public:
  MNDOMethod method;
  std::shared_ptr<MNDOMethodWrapper> mndoMethodWrapper;
  std::shared_ptr<Core::Calculator> polymorphicMethodWrapper;
  ElementParameters elementParameters;

  void SetUp() override {
    mndoMethodWrapper = std::make_shared<MNDOMethodWrapper>();
    polymorphicMethodWrapper = std::make_shared<MNDOMethodWrapper>();

    method.setMaxIterations(10000);
    method.setConvergenceCriteria(1e-8);
    mndoMethodWrapper->settings().modifyInt(Utils::SettingsNames::maxIterations, 10000);
    mndoMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
    mndoMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
    mndoMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);
    polymorphicMethodWrapper->settings().modifyInt(Utils::SettingsNames::maxIterations, 10000);
    polymorphicMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
    polymorphicMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
    polymorphicMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);
    Utils::Log::startConsoleLogging(Utils::SettingsNames::LogLevels::none);
  }
};
TEST_F(AMNDOCalculation, HasTheCorrectNumberOfAOsAfterInitialization) {
  std::stringstream ssH("4\n\n"
                        "C      0.0000000000    0.0000000000    0.0000000000\n"
                        "P      0.0529177211   -0.3175063264    0.2645886053\n"
                        "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                        "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ssH);
  method.setStructure(as, parameters_mndo);
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(10));
}

TEST_F(AMNDOCalculation, GetsCorrectOneElectronPartOfFockForH) {
  std::stringstream ssH("1\n\n"
                        "H -1.0 0.2 -0.3\n");
  auto as = Utils::XyzStreamHandler::read(ssH);
  method.setStructure(as, parameters_mndo);

  method.calculateDensityIndependentQuantities();
  ASSERT_THAT(method.getOneElectronMatrix().getMatrix()(0, 0) * ev_per_hartree, DoubleNear(-11.906276, 1e-4));
}

TEST_F(AMNDOCalculation, GetsCorrectOneElectronPartOfFockForC) {
  std::stringstream ssC("1\n\n"
                        "C 0.0 0.0 0.0\n");
  auto as = Utils::XyzStreamHandler::read(ssC);
  method.setStructure(as, parameters_mndo);
  method.setScfMixer(Utils::scf_mixer_t::none);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 4, 4> expected;
  expected << -52.279745, 0, 0, 0, 0, -39.205558, 0, 0, 0, 0, -39.205558, 0, 0, 0, 0, -39.205558;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-4));
}

TEST_F(AMNDOCalculation, GetsCorrectOneElectronPartOfFockForH2) {
  std::stringstream ssH2("2\n\n"
                         "H -0.529177  0.105835 -0.158753\n"
                         "H -0.105835  0.105835 -0.158753\n");
  auto as = Utils::XyzStreamHandler::read(ssH2);
  method.setStructure(as, parameters_mndo);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 2, 2> expected;
  expected << -23.925433, 0, -5.885148, -23.925433;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-4));
}

TEST_F(AMNDOCalculation, GetsCorrectOneElectronPartOfFockForCC) {
  std::stringstream ssC2("2\n\n"
                         "C 0.0 0.0 0.0\n"
                         "C 0.423342  0.0529177 -0.0529177\n");
  auto as = Utils::XyzStreamHandler::read(ssC2);
  method.setStructure(as, parameters_mndo);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 8, 8> expected;
  expected << -98.232349, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.533677, -81.808790, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -0.816709, -0.154931, -80.588709, 0.0, 0.0, 0.0, 0.0, 0.0, 0.816709, 0.154931, 0.019366, -80.588709, 0.0, 0.0,
      0.0, 0.0, -17.012504, -4.999932, -0.624991, 0.624991, -98.232349, 0.0, 0.0, 0.0, 4.999932, -4.103025, 0.299549,
      -0.299549, 6.533677, -81.808790, 0.0, 0.0, 0.624991, 0.299549, -6.461975, -0.037444, 0.816709, -0.154931,
      -80.588709, 0.0, -0.624991, -0.299549, -0.037444, -6.461975, -0.816709, 0.154931, 0.019366, -80.588709;
  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-3));
}

TEST_F(AMNDOCalculation, GetsCorrectOneElectronPartOfFockForCH2O) {
  std::stringstream ss("4\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "O      0.0529177211   -0.3175063264    0.2645886053\n"
                       "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                       "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_mndo);

  method.calculateDensityIndependentQuantities();
  Eigen::Matrix<double, 10, 10> expected;
  expected << -151.976060, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.772741, -127.799476, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 9.079954, 0.379901, -128.887601, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.977843, -0.394638, 1.378475,
      -128.459407, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -20.183154, -1.356416, 8.138497, -6.782081, -173.278431, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.857516, -13.744236, -0.695036, 0.579196, 2.852947, -147.565695, 0.0, 0.0, 0.0, 0.0, -5.145099,
      -0.695036, -9.689862, -3.475178, -7.839927, 0.446619, -148.394089, 0.0, 0.0, 0.0, 4.287582, 0.579196, -3.475178,
      -10.964094, 6.947827, -0.425350, 1.061814, -148.160132, 0.0, 0.0, -10.180955, 3.701731, -0.740346, 1.110519,
      -9.351522, 5.002242, -3.637994, 3.637994, -134.207785, 0.0, -12.119934, 1.029453, -1.029453, 1.544180, -12.681354,
      2.134817, -5.692845, 5.692845, -5.885149, -145.194765;

  ASSERT_TRUE(method.getOneElectronMatrix().getMatrix().isApprox(expected / ev_per_hartree, 1e-4));
}

TEST_F(AMNDOCalculation, GetsCorrectEigenvaluesForC) {
  std::stringstream ss("1\n\n"
                       "C     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_mndo);
  method.setScfMixer(Utils::scf_mixer_t::none);

  method.convergedCalculation();

  auto eigenvalues = method.getSingleParticleEnergies().getRestrictedEnergies();

  std::vector<double> expected = {-19.53974, -7.61556, 0.36444, 0.36444};
  for (unsigned i = 0; i < eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i] * ev_per_hartree, DoubleNear(expected[i], 1e-4));
  }
}

TEST_F(AMNDOCalculation, GetsCorrectUnrestrictedEigenvaluesForBoron) {
  std::stringstream ss("1\n\n"
                       "B     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_mndo);

  method.setUnrestrictedCalculation(true);
  method.setSpinMultiplicity(2);
  method.convergedCalculation(Utils::derivativeType::none);

  auto alphaEigenvalues = method.getSingleParticleEnergies().getAlphaEnergies();
  auto betaEigenvalues = method.getSingleParticleEnergies().getBetaEnergies();

  // From Mopac
  double mopacEnergy = -64.31595; // eV
  std::vector<double> alphaExpected = {-16.20713, -5.81169, 1.54831, 1.54831};
  std::vector<double> betaExpected = {-14.39713, 2.04831, 2.04831, 3.04831};

  ASSERT_THAT(method.getEnergy() * ev_per_hartree, DoubleNear(mopacEnergy, 1e-5));
  for (unsigned i = 0; i < alphaEigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalues " + std::to_string(i) + ":");
    EXPECT_THAT(alphaEigenvalues[i] * ev_per_hartree, DoubleNear(alphaExpected[i], 1e-5));
    EXPECT_THAT(betaEigenvalues[i] * ev_per_hartree, DoubleNear(betaExpected[i], 1e-5));
  }
}

TEST_F(AMNDOCalculation, IsIndependentOfOrderOfAtoms) {
  std::stringstream ss("3\n\n"
                       "O      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.4005871485    0.3100978455    0.0000000000\n"
                       "H     -0.4005871485    0.3100978455    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();
  double energy1 = method.getElectronicEnergy();

  std::stringstream ss2("3\n\n"
                        "H      0.4005871485    0.3100978455    0.0000000000\n"
                        "H     -0.4005871485    0.3100978455    0.0000000000\n"
                        "O      0.0000000000    0.0000000000    0.0000000000\n");
  auto as2 = Utils::XyzStreamHandler::read(ss2);
  method.setStructure(as2, parameters_mndo);

  method.convergedCalculation();
  double energy2 = method.getElectronicEnergy();

  ASSERT_THAT(energy1, DoubleNear(energy2, 1e-10));
}

TEST_F(AMNDOCalculation, GetsCorrectTotalEnergyForH2) {
  std::stringstream ssH2("2\n\n"
                         "H -0.529177  0.105835 -0.158753\n"
                         "H -0.105835  0.105835 -0.158753\n");
  auto as = Utils::XyzStreamHandler::read(ssH2);
  method.setStructure(as, parameters_mndo);
  method.setScfMixer(Utils::scf_mixer_t::none);
  method.convergedCalculation();

  auto eigenvalues = method.getSingleParticleEnergies().getRestrictedEnergies();
  std::vector<double> expected = {-17.37700, 6.41245};

  for (unsigned i = 0; i < eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i] * ev_per_hartree, DoubleNear(expected[i], 1e-4));
  }
  ASSERT_THAT(method.getElectronicEnergy() * ev_per_hartree, DoubleNear(-47.18758, 1e-4));
  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(20.20667, 1e-4));
}

TEST_F(AMNDOCalculation, GetsCorrectTotalEnergyForH2AtOtherDistance) {
  std::stringstream ss("2\n\n"
                       "H     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();

  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(6.35835, 1e-4));
  ASSERT_THAT(method.getElectronicEnergy() * ev_per_hartree, DoubleNear(-28.12803, 1e-3));
}

TEST_F(AMNDOCalculation, GetsCorrectTotalEnergyForC2) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "C     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();

  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(100.49139, 1e-4));
  ASSERT_THAT(method.getElectronicEnergy() * ev_per_hartree, DoubleNear(-340.03345, 1e-3));
}

TEST_F(AMNDOCalculation, GetsCorrectRepulsionEnergyForCH) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);

  method.setStructure(as, parameters_mndo);
  method.calculate(Utils::derivativeType::none);

  ASSERT_THAT(method.getRepulsionEnergy() * ev_per_hartree, DoubleNear(25.27860, 1e-4));
}

TEST_F(AMNDOCalculation, GetsCorrectGradientsForMethane) {
  /* Structure optimized with MOPAC 2016
   *
   *
   * MNDO RHF 1SCF PRECISE GRADIENTS
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
  method.setStructure(structure, parameters_mndo);
  method.convergedCalculation(Utils::derivativeType::first);

  Eigen::RowVector3d cOneGradient(0.000010, 0.000011, 0.000012);
  Eigen::RowVector3d hOneGradient(-8.984101, -8.982663, -8.983207);
  Eigen::RowVector3d hTwoGradient(8.956718, 8.957105, -8.935515);
  Eigen::RowVector3d hThreeGradient(8.957231, -8.933676, 8.958441);
  Eigen::RowVector3d hFourGradient(-8.929858, 8.959223, 8.960270);
  cOneGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hOneGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hTwoGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hThreeGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;
  hFourGradient *= Utils::Constants::angstrom_per_bohr * Utils::Constants::hartree_per_kCalPerMol;

  ASSERT_THAT(method.getEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-185.07194, 1e-3));
  ASSERT_TRUE((method.getGradients().row(0) - cOneGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(1) - hOneGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(2) - hTwoGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(3) - hThreeGradient).norm() < 1e-4);
  ASSERT_TRUE((method.getGradients().row(4) - hFourGradient).norm() < 1e-4);
}

TEST_F(AMNDOCalculation, GetsZeroGradientForOptimizedMethane) {
  /* Structure optimized with MOPAC 2016
   *
   *
   * MNDO RHF OPT PRECISE GRADIENTS
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
                           "C     0.00003086    0.00006509   -0.00021557\n"
                           "H     0.63772679    0.63840547    0.63682059\n"
                           "H    -0.63773658   -0.63845649    0.63671348\n"
                           "H    -0.63765374    0.63821789   -0.63638609\n"
                           "H     0.63755701   -0.63839156   -0.63640383\n");

  auto structure = Utils::XyzStreamHandler::read(stream);
  method.setStructure(structure, parameters_mndo);
  method.convergedCalculation(Utils::derivativeType::first);

  ASSERT_THAT(method.getEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-185.09228, 1e-3));
  ASSERT_TRUE(method.getGradients().row(0).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(1).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(2).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(3).norm() < 1e-3);
  ASSERT_TRUE(method.getGradients().row(4).norm() < 1e-3);
}

TEST_F(AMNDOCalculation, GetsCorrectTotalEnergyForEthanol) {
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
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();

  ASSERT_THAT(method.getEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-663.83891, 7e-2));
  ASSERT_THAT(method.getRepulsionEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(1133.01130, 7e-2));
  ASSERT_THAT(method.getElectronicEnergy() * Utils::Constants::ev_per_hartree, DoubleNear(-1796.85021, 7e-2));
}

TEST_F(AMNDOCalculation, GetsSameTotalEnergyWithAndWithoutWrapper) {
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
  mndoMethodWrapper->setStructure(as);
  auto result = mndoMethodWrapper->calculate("");
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();

  ASSERT_THAT(method.getEnergy(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-5));
}

TEST_F(AMNDOCalculation, GetsSameTotalEnergyWithPolymorphicMethodAndWithMethod) {
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
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();

  ASSERT_THAT(method.getEnergy(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-5));
}

TEST_F(AMNDOCalculation, GetsSameTotalEnergyWithDynamicallyLoadedMethod) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
  dynamicallyLoadedMethodWrapper->settings().modifyInt(Utils::SettingsNames::maxIterations, 10000);
  dynamicallyLoadedMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);
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
  method.setStructure(as, parameters_mndo);
  method.convergedCalculation();

  ASSERT_THAT(result.get<Utils::Property::Energy>() * Utils::Constants::ev_per_hartree, DoubleNear(-663.83891, 7e-2));
  ASSERT_THAT(method.getEnergy(), DoubleNear(result.get<Utils::Property::Energy>(), 1e-5));
  ASSERT_EQ(result.get<Utils::Property::SuccessfulCalculation>(), true);
}

TEST_F(AMNDOCalculation, MethodWrapperCanBeCloned) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
  dynamicallyLoadedMethodWrapper->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-8);
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);

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

TEST_F(AMNDOCalculation, StructureIsCorrectlyCloned) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);

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

TEST_F(AMNDOCalculation, ClonedMethodCanCalculate) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);

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

TEST_F(AMNDOCalculation, ClonedMethodCanCalculateGradients) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);

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

TEST_F(AMNDOCalculation, ClonedMethodCanCalculateGradientsWithDifferentNumberCores) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
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
}

TEST_F(AMNDOCalculation, ClonedMethodCopiesResultsCorrectly) {
  auto& moduleManager = Core::ModuleManager::getInstance();
  auto programPath = boost::dll::program_location();
  auto libPath = programPath.parent_path() / "sparrow";
  try {
    moduleManager.load(libPath);
  }
  catch (const std::runtime_error& e) {
    // Do nothing if module is already loaded.
  }
  auto dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("MNDO");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_mndo);

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

TEST_F(AMNDOCalculation, AtomCollectionCanBeReturned) {
  std::stringstream ssH("4\n\n"
                        "C      0.0000000000    0.0000000000    0.0000000000\n"
                        "C      0.0529177211   -0.3175063264    0.2645886053\n"
                        "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                        "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto structure = Utils::XyzStreamHandler::read(ssH);
  mndoMethodWrapper->setStructure(structure);
  ASSERT_EQ(structure.getPositions(), mndoMethodWrapper->getStructure()->getPositions());
  for (int i = 0; i < structure.getElements().size(); ++i)
    ASSERT_EQ(structure.getElements()[i], mndoMethodWrapper->getStructure()->getElements()[i]);
}
} // namespace Sparrow
} // namespace Scine
