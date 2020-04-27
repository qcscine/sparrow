/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "../MethodsTests/parameters_location.h"
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
    pm6->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);
    pm6->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-9);
    pm6->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, parameters_root);
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

TEST_F(CommonDipoleCalculation, NddoDipoleIsCorrectForBigOrganicMolecule) {
  pm6->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-5);
  std::stringstream BigMoleculeSS("113\n\n"
                                  " O         15.68894000     23.55459000      9.30853000\n"
                                  " O         15.81287000     29.65719000      8.45683000\n"
                                  " O         13.40468000     23.13785000      8.04487000\n"
                                  " O         13.66559000     31.43847000      8.92988000\n"
                                  " O         11.84316000     25.10056000      6.75053000\n"
                                  " O         11.27957000     22.92183000      4.48568000\n"
                                  " O         15.33904000     27.93590000     11.43104000\n"
                                  " O         14.89455000     24.61996000      2.26223000\n"
                                  " O         16.96517000     23.15749000      4.02856000\n"
                                  " O         18.19509000     25.47139000      4.66759000\n"
                                  " O         14.27744000     27.96112000      8.23079000\n"
                                  " O         14.19171000     20.98492000      8.01749000\n"
                                  " O         12.29558000     26.86046000      5.36613000\n"
                                  " O         17.56628000     25.65774000      2.46106000\n"
                                  " N         14.56699000     30.07902000     11.30982000\n"
                                  " C         14.68174000     29.22284000      8.56648000\n"
                                  " C         15.45143000     24.38602000      8.15052000\n"
                                  " C         13.50110000     30.02776000      9.09556000\n"
                                  " C         14.50265000     23.55220000      7.20651000\n"
                                  " C         13.97563000     24.22371000      5.91584000\n"
                                  " C         13.34539000     29.74221000     10.59931000\n"
                                  " C         12.43400000     24.18275000      5.79537000\n"
                                  " C         11.84834000     24.25174000      4.36145000\n"
                                  " C         15.49332000     29.13246000     11.64411000\n"
                                  " C         12.71374000     24.40802000      3.13631000\n"
                                  " C         14.15346000     24.71850000      3.47199000\n"
                                  " C         14.69562000     23.79842000      4.58623000\n"
                                  " C         16.22023000     24.00886000      4.53087000\n"
                                  " C         16.81442000     25.35470000      4.96225000\n"
                                  " C         16.55976000     25.70865000      6.39811000\n"
                                  " C         15.88346000     26.84571000      6.71155000\n"
                                  " C         15.33940000     27.00498000      8.10391000\n"
                                  " C         14.75248000     25.68359000      8.61465000\n"
                                  " C         16.82938000     24.68844000      7.50189000\n"
                                  " C         17.81199000     25.26213000      8.54280000\n"
                                  " C         17.53019000     23.39882000      7.06519000\n"
                                  " C         14.47601000     22.30914000      4.27596000\n"
                                  " C         15.59993000     27.96092000      5.75969000\n"
                                  " C         11.71549000     22.82509000      5.85431000\n"
                                  " C         13.43435000     21.82880000      8.41849000\n"
                                  " C         12.34736000     21.56724000      9.38923000\n"
                                  " C         11.33985000     22.48099000      9.67784000\n"
                                  " C         10.33828000     22.15278000     10.58677000\n"
                                  " C         10.35451000     20.93427000     11.24191000\n"
                                  " C         11.35932000     20.02687000     10.96377000\n"
                                  " C         12.34265000     20.33518000     10.03119000\n"
                                  " C         12.17664000     30.47797000     11.22155000\n"
                                  " C         10.95065000     30.57629000     10.56476000\n"
                                  " C          9.87791000     31.25309000     11.13399000\n"
                                  " C         10.00793000     31.83076000     12.37913000\n"
                                  " C         11.21090000     31.74058000     13.04578000\n"
                                  " C         12.27968000     31.06453000     12.48174000\n"
                                  " C         11.81973000     26.40896000      6.38046000\n"
                                  " C         11.13872000     27.18481000      7.46687000\n"
                                  " C         16.71749000     29.66167000     12.28018000\n"
                                  " C         16.66427000     30.72002000     13.18879000\n"
                                  " C         17.81868000     31.13753000     13.83690000\n"
                                  " C         19.03548000     30.55095000     13.53529000\n"
                                  " C         19.10009000     29.51933000     12.61263000\n"
                                  " C         17.94140000     29.06596000     11.99062000\n"
                                  " C         18.42320000     25.60374000      3.32230000\n"
                                  " C         19.89529000     25.74309000      3.07987000\n"
                                  " H         15.16627000     23.88279000     10.05888000\n"
                                  " H         12.61120000     29.71350000      8.54163000\n"
                                  " H         15.06941000     22.69588000      6.86985000\n"
                                  " H         14.19722000     25.27422000      6.01922000\n"
                                  " H         13.14478000     28.67651000     10.73243000\n"
                                  " H         14.76787000     31.06382000     11.43566000\n"
                                  " H         11.06762000     25.00235000      4.23170000\n"
                                  " H         14.21517000     25.75778000      3.80641000\n"
                                  " H         15.67245000     25.21900000      2.28643000\n"
                                  " H         16.27859000     26.06971000      4.35425000\n"
                                  " H         16.16746000     27.33282000      8.71910000\n"
                                  " H         11.31303000     23.45050000      9.20152000\n"
                                  " H          9.53106000     22.84162000     10.79023000\n"
                                  " H          9.57923000     20.69386000     11.95893000\n"
                                  " H         11.37607000     19.07270000     11.47077000\n"
                                  " H         13.11452000     19.60724000      9.81988000\n"
                                  " H         10.82137000     30.13843000      9.58713000\n"
                                  " H          8.93238000     31.34348000     10.61877000\n"
                                  " H          9.17891000     32.35272000     12.83531000\n"
                                  " H         11.31381000     32.21170000     14.00998000\n"
                                  " H         13.20619000     31.03173000     13.03773000\n"
                                  " H         15.74161000     31.22954000     13.42928000\n"
                                  " H         17.77621000     31.93087000     14.57134000\n"
                                  " H         19.94155000     30.92231000     13.99551000\n"
                                  " H         20.06383000     29.08383000     12.38096000\n"
                                  " H         18.00275000     28.26321000     11.26730000\n"
                                  " H         12.31093000     25.22577000      2.52791000\n"
                                  " H         12.65197000     23.52449000      2.50430000\n"
                                  " H         14.68841000     25.73408000      9.70223000\n"
                                  " H         13.71185000     25.65042000      8.29542000\n"
                                  " H         18.78044000     25.48765000      8.09474000\n"
                                  " H         17.96953000     24.56421000      9.36730000\n"
                                  " H         18.49572000     23.60736000      6.60120000\n"
                                  " H         17.69632000     22.73454000      7.91451000\n"
                                  " H         16.95213000     22.81777000      6.36140000\n"
                                  " H         14.87233000     21.66911000      5.06649000\n"
                                  " H         13.43326000     22.04626000      4.14831000\n"
                                  " H         14.97357000     22.01508000      3.35062000\n"
                                  " H         16.02277000     27.80580000      4.77422000\n"
                                  " H         14.52635000     28.09181000      5.63565000\n"
                                  " H         16.04339000     28.88920000      6.11321000\n"
                                  " H         10.87952000     22.80916000      6.55984000\n"
                                  " H         12.34095000     21.95832000      6.06644000\n"
                                  " H         11.10038000     28.24006000      7.20347000\n"
                                  " H         10.12125000     26.82392000      7.61594000\n"
                                  " H         11.68846000     27.07682000      8.40192000\n"
                                  " H         17.48064000     26.17473000      9.01138000\n"
                                  " H         20.08808000     25.86160000      2.01213000\n"
                                  " H         20.27627000     26.62537000      3.59684000\n"
                                  " H         20.41636000     24.85112000      3.42576000\n"
                                  " H         13.77125000     31.61639000      7.98284000\n");

  auto BigMolecule = Utils::XyzStreamHandler::read(BigMoleculeSS);
  pm6->setStructure(BigMolecule);

  auto res = pm6->calculate("");
  auto dipole = res.get<Utils::Property::Dipole>().norm() * 2.541746;

  // Reference value obtained with MOPAC (as implemented in ADF 2016.107)
  ASSERT_THAT(dipole, DoubleNear(11.693, 0.005));
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
} // namespace Sparrow
} // namespace Scine
