/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMomentCalculator.h>
#include <Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

using namespace testing;

namespace Scine {
namespace Sparrow {

class CommonDipoleCalculation : public Test {
 public:
  std::shared_ptr<PM6MethodWrapper> pm6;
  std::shared_ptr<DFTB2MethodWrapper> dftb2;
  Utils::AtomCollection HF, CO2, Et;
  std::vector<double> ChargesHF;
  std::vector<double> ChargesCO2;
  std::vector<double> ChargesEt;
  std::vector<double> eigenValues_;
  Eigen::MatrixXd eigenVectors_;

  void assignStructure(const Utils::AtomCollection& structure);

 protected:
  void SetUp() override {
    pm6 = std::make_shared<PM6MethodWrapper>();
    pm6->setLog(Core::Log::silent());
    pm6->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);
    pm6->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-9);
    dftb2 = std::make_shared<DFTB2MethodWrapper>();
    dftb2->setLog(Core::Log::silent());
    dftb2->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);
    dftb2->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-9);
    std::stringstream HFss("2\n\n"
                           "H        0.000000000     0.000000000     0.000000000\n"
                           "F        0.965548748     0.000000000     0.000000000\n");
    HF = Utils::XyzStreamHandler::read(HFss);
    std::stringstream CO2ss("3\n\n"
                            "C        0.013554240     0.017673604     0.000000000\n"
                            "O        1.183652472     0.017673604     0.000000000\n"
                            "O       -1.156896800     0.025294937     0.000000000\n");
    CO2 = Utils::XyzStreamHandler::read(CO2ss);

    std::stringstream Ethanolss("9\n\n"
                                "C       -0.025722743    -0.038148774     0.007736703\n"
                                "C        1.494872990    -0.038148774     0.007736703\n"
                                "O        2.017472251     1.299308464     0.007736703\n"
                                "H       -0.418489864     0.516363236     0.871446070\n"
                                "H       -0.435916359     0.434351771    -0.891900626\n"
                                "H       -0.429207593    -1.056100478     0.060014382\n"
                                "H        1.909273199    -0.583326143    -0.859094289\n"
                                "H        1.916182172    -0.455131719     0.942989512\n"
                                "H        1.715886732     1.789449608    -0.779284087\n");
    Et = Utils::XyzStreamHandler::read(Ethanolss);

    // These charges have been obtained with MOPAC (using PM6)
    ChargesHF = {0.270573, -0.270573};
    ChargesCO2 = {0.639403, -0.319550, -0.319854};
    ChargesEt = {-0.537402, 0.052114, -0.555356, 0.172130, 0.157952, 0.162459, 0.104917, 0.138407, 0.304779};
  }
};

void CommonDipoleCalculation::assignStructure(const Utils::AtomCollection& structure) {
  pm6->setStructure(structure);
}

TEST_F(CommonDipoleCalculation, CalculateDipoleCorrectlyHF) {
  CommonDipoleCalculation::assignStructure(HF);
  auto res = pm6->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();
  // Reference value obtained with MOPAC (as implemented in ADF 2016.107)
  ASSERT_THAT(dipole.norm() * 2.541746, DoubleNear(1.436, 0.001));
}

TEST_F(CommonDipoleCalculation, CalculateDipoleCorrectlyCO2) {
  CommonDipoleCalculation::assignStructure(CO2);
  auto res = pm6->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();

  // Reference value obtained with MOPAC (as implemented in ADF 2016.107)
  // Note that the structure is not really an equilibrium structure - hence a nonzero dipole moment
  ASSERT_THAT(dipole.norm() * 2.541746, DoubleNear(0.013, 0.001));
}

TEST_F(CommonDipoleCalculation, CalculateDipoleCorrectlyEthanol) {
  pm6->setStructure(Et);
  auto res = pm6->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();

  // Reference value obtained with MOPAC (as implemented in ADF 2016.107)
  ASSERT_THAT(dipole.norm() * 2.541746, DoubleNear(2.133, 0.001));
}

TEST_F(CommonDipoleCalculation, NDDODipoleIsRotationallyInvariant) {
  pm6->setStructure(Et);
  auto res = pm6->calculate("");
  auto dipole1 = res.get<Utils::Property::Dipole>().norm();

  // 30° rotation on x axis
  std::stringstream EtRotated("9\n\n"
                              "C -0.02572274 -0.03690616 -0.01237421\n"
                              "C 1.49487299 -0.03690616 -0.01237421\n"
                              "O 2.01747225  1.12136579  0.65635441\n"
                              "H -0.41848986  0.01146064  1.01287605\n"
                              "H -0.43591636  0.82210998 -0.55523271\n"
                              "H -0.42920759 -0.94461703 -0.47607626\n"
                              "H 1.90927320 -0.07562811 -1.03566055\n"
                              "H 1.91618217 -0.86565039  0.58908701\n"
                              "H 1.71588673  1.93935086  0.21984499\n");
  auto Et2 = Utils::XyzStreamHandler::read(EtRotated);
  pm6->setStructure(Et2);

  auto res2 = pm6->calculate("");
  auto dipole2 = res2.get<Utils::Property::Dipole>().norm();

  ASSERT_THAT(dipole1, DoubleNear(dipole2, 1e-6));

  // 30° rotation on x axis and translation by 5 Angstrom in z direction
  std::stringstream EtRotatedTranslated("9\n\n"
                                        "C -0.02572274 -0.03690616 4.98762579\n"
                                        "C 1.49487299 -0.03690616 4.98762579\n"
                                        "O 2.01747225  1.12136579  5.65635441\n"
                                        "H -0.41848986  0.01146064  6.01287605\n"
                                        "H -0.43591636  0.82210998 4.44476729\n"
                                        "H -0.42920759 -0.94461703 4.52392374\n"
                                        "H 1.90927320 -0.07562811 3.96433945\n"
                                        "H 1.91618217 -0.86565039  5.58908701\n"
                                        "H 1.71588673  1.93935086  5.21984499\n");
  auto Et3 = Utils::XyzStreamHandler::read(EtRotatedTranslated);
  pm6->setStructure(Et3);

  auto res3 = pm6->calculate("");
  auto dipole3 = res3.get<Utils::Property::Dipole>().norm();

  ASSERT_THAT(dipole1, DoubleNear(dipole3, 1e-6));
}

TEST_F(CommonDipoleCalculation, CysteineDipoleTest) {
  std::stringstream cysteine("14\n\n"
                             "C     -0.2039680165    1.2296177212   -1.7601331753\n"
                             "C     -0.0739379306    2.3710252329   -2.7675709632\n"
                             "C     -0.9367365741    2.0463909311   -4.0072129656\n"
                             "N     -0.4890026632    3.7059180651   -2.2691947201\n"
                             "S      0.9245770748    1.6284396638   -0.3926286873\n"
                             "H      0.2683791134    4.1540294814   -1.7547204163\n"
                             "H      0.5798884925    0.6490262236    0.4712106159\n"
                             "H     -1.2436915763    1.1116634977   -1.4010351046\n"
                             "H      0.9912846668    2.4350033663   -3.1381132275\n"
                             "H      0.0671139651    0.2605154451   -2.2327400800\n"
                             "H     -1.2781819985    3.6324983356   -1.6282366244\n"
                             "O     -1.1677279714    0.9357684341   -4.4147456623\n"
                             "O     -1.4295571152    3.1053969018   -4.7093924764\n"
                             "H     -1.2373572757    4.0001923525   -4.2901399260\n");
  auto cysteineStructure = Utils::XyzStreamHandler::read(cysteine);
  pm6->setStructure(cysteineStructure);

  auto res = pm6->calculate("");
  auto dipole = res.get<Utils::Property::Dipole>().norm() * 2.541746;

  // Reference value obtained with MOPAC (as implemented in mopac 2016)
  ASSERT_THAT(dipole, DoubleNear(4.199, 0.005));
}

TEST_F(CommonDipoleCalculation, DipoleOfOrganicMoleculesWithdOrbitalsCorrectlyCalculated) {
  std::stringstream PSCagess("8\n\n"
                             "S     -0.1879670442   -0.3701558968    1.8628954719\n"
                             "S     -1.0982691763    0.4157385980   -1.2189132042\n"
                             "S     -2.8404586734   -1.4313197107    2.3714571620\n"
                             "S     -3.8507608051   -0.8454252170   -0.5103515145\n"
                             "P     -3.4645106974    0.3030010654    1.2652796617\n"
                             "P     -1.4256258173    1.2371974877    0.7404185869\n"
                             "P     -1.8678459689   -1.5686821956   -0.9176428609\n"
                             "P     -1.3194732153   -2.3026785840    1.5170325276\n");
  auto PSCage = Utils::XyzStreamHandler::read(PSCagess);
  pm6->setStructure(PSCage);

  auto res = pm6->calculate("");
  auto dipole = res.get<Utils::Property::Dipole>().norm() * 2.541746;

  // Reference value obtained with MOPAC (as implemented in mopac 2016)
  ASSERT_THAT(dipole, DoubleNear(0.961, 0.005));
}

TEST_F(CommonDipoleCalculation, DipoleOfTransitionMetalsCorrectlyCalculated) {
  std::stringstream ironComplex1ss("9\n\n"
                                   "O          1.01154        0.40762       -0.35393\n"
                                   "N          2.18761        0.53843       -0.50286\n"
                                   "Fe         2.88109        1.89561       -1.76016\n"
                                   "C          2.86783        0.52302       -3.25074\n"
                                   "O          2.86026       -0.24100       -4.08065\n"
                                   "C          2.89268        3.27041       -0.27155\n"
                                   "O          2.89890        4.03581        0.55706\n"
                                   "N          4.73037        1.37634       -1.29674\n"
                                   "O          5.65279        1.93052       -1.81095\n");
  auto feComplex1 = Utils::XyzStreamHandler::read(ironComplex1ss);
  pm6->setStructure(feComplex1);

  auto res = pm6->calculate("");
  auto dipole = res.get<Utils::Property::Dipole>().norm() * 2.541746;

  // Reference value obtained with MOPAC (as implemented in mopac 2016)
  ASSERT_THAT(dipole, DoubleNear(2.455, 0.005));
}

TEST_F(CommonDipoleCalculation, DFTBCalculatesDipoleCorrectly) {
  dftb2->setStructure(Et);
  dftb2->setRequiredProperties(Utils::Property::Dipole);
  auto dipole = dftb2->calculate("").get<Utils::Property::Dipole>();

  // Reference value obtained with DFTB+
  Eigen::Vector3d referenceDipoleDftbPlus(-0.36023764, -0.21653604, -0.38972594);
  EXPECT_THAT(dipole(0), DoubleNear(referenceDipoleDftbPlus(0), 1e-6));
  EXPECT_THAT(dipole(1), DoubleNear(referenceDipoleDftbPlus(1), 1e-6));
  EXPECT_THAT(dipole(2), DoubleNear(referenceDipoleDftbPlus(2), 1e-6));
}

// Test that if you calculate the dipole with the dipole matrix, you get the same
// as with a Mulliken pop analysis.
TEST_F(CommonDipoleCalculation, DFTBCalculatedDipoleMatrixCorrectly) {
  dftb::DFTB2 dftb2Method;
  Core::Log log = Core::Log::silent();
  auto dipoleMatrixCalc = DFTBDipoleMatrixCalculator<dftb::DFTB2>::create(dftb2Method);
  auto dipoleMomentCalc = DFTBDipoleMomentCalculator<dftb::DFTB2>(dftb2Method);
  dftb2Method.setAtomCollection(Et);
  dftb2Method.initializeFromParameterPath("mio-1-1");
  dftb2Method.convergedCalculation(log, Utils::Derivative::None);
  auto mos = dftb2Method.getMolecularOrbitals().restrictedMatrix();
  int nOcc = dftb2Method.getElectronicOccupation().numberRestrictedElectrons() / 2;

  dipoleMatrixCalc->fillDipoleMatrix({0, 0, 0});
  auto dipoleMatrixMO = dipoleMatrixCalc->getMODipoleMatrix();

  auto dipole = dipoleMomentCalc.calculate();

  Eigen::VectorXd coreCharges(9);
  coreCharges << 4, 4, 6, 1, 1, 1, 1, 1, 1;
  Eigen::RowVector3d dipAtom;
  dipAtom = coreCharges.transpose() * dftb2Method.getPositions();
  EXPECT_THAT(dipAtom.x() - dipoleMatrixMO[0].diagonal().head(nOcc).sum() * 2.0, DoubleNear(dipole.x(), 1e-6));
  EXPECT_THAT(dipAtom.y() - dipoleMatrixMO[1].diagonal().head(nOcc).sum() * 2.0, DoubleNear(dipole.y(), 1e-6));
  EXPECT_THAT(dipAtom.z() - dipoleMatrixMO[2].diagonal().head(nOcc).sum() * 2.0, DoubleNear(dipole.z(), 1e-6));
}

} // namespace Sparrow
} // namespace Scine
