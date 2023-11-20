/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/ConvergenceAccelerators/ConvergenceAcceleratorFactory.h>
#include <gmock/gmock.h>
#include <omp.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using Utils::Derivative;
using Utils::scf_mixer_t;

class ADFTB2Calculation : public Test {
 public:
  dftb::DFTB2 method;
  std::shared_ptr<Core::Calculator> dynamicallyLoadedMethodWrapper;
  std::shared_ptr<DFTB2MethodWrapper> calculator;

  Core::Log log;

  void SetUp() override {
    log = Core::Log::silent();
    calculator = std::make_shared<DFTB2MethodWrapper>();
    calculator->setLog(log);
    auto& moduleManager = Core::ModuleManager::getInstance();
    dynamicallyLoadedMethodWrapper = moduleManager.get<Core::Calculator>("DFTB2");
    dynamicallyLoadedMethodWrapper->setLog(log);
    dynamicallyLoadedMethodWrapper->settings().modifyString(Utils::SettingsNames::methodParameters, "mio-1-1");
    calculator->settings().modifyString(Utils::SettingsNames::methodParameters, "mio-1-1");
  }
};

TEST_F(ADFTB2Calculation, HasTheCorrectNumberOfAtomsAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setMolecularCharge(3);
  method.setAtomCollection(as);
  method.initializeFromParameterPath("mio-1-1");
  ASSERT_THAT(method.getNumberAtoms(), Eq(2));
}

TEST_F(ADFTB2Calculation, HasTheCorrectNumberOfOrbitalsAfterInitialization) {
  std::stringstream ss("2\n\n"
                       "C     0.0000000000    0.0000000000   -0.0000000000\n"
                       "H     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setMolecularCharge(3);
  method.setAtomCollection(as);
  method.initializeFromParameterPath("mio-1-1");
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(5));
}

TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForC) {
  std::stringstream ss("1\n\n"
                       "C     2.0000000000    0.0000000000    0.0000000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setAtomCollection(as);
  method.initializeFromParameterPath("mio-1-1");
  method.calculate(Derivative::First, log);

  // Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(1));
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(4));

  // Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-1.3984936602, 1e-8));

  // Check eigenvalues
  auto eigenvalues = method.getSingleParticleEnergies().getRestrictedEnergies();
  std::vector<double> expected = {-0.50489172, -0.19435511, -0.19435511, -0.19435511};
  for (unsigned i = 0; i < eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i], DoubleNear(expected[i], 1e-10));
  }

  // Check force
  ASSERT_THAT(method.getGradients().row(0).norm(), DoubleEq(0));
}

TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForCH4) {
  std::stringstream ss("5\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.6287000000    0.6287000000    0.6287000000\n"
                       "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                       "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                       "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);

  method.setScfMixer(scf_mixer_t::fock_diis);
  method.setMaxIterations(10000);
  method.setConvergenceCriteria({{}, 1e-8});

  method.setAtomCollection(as);
  method.initializeFromParameterPath("mio-1-1");
  method.calculate(Derivative::First, log);

  // Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(5));
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(8));

  // Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-3.2256725551, 1e-10));

  // Check atomic charges
  std::vector<double> expectedCharges = {-0.30575605, 0.07643901, 0.07643901, 0.07643901};
  for (unsigned long long i = 0; i < expectedCharges.size(); i++) {
    SCOPED_TRACE("... for the charge " + std::to_string(i) + ":");
    EXPECT_THAT(method.getAtomicCharges()[i], DoubleNear(expectedCharges[i], 1e-6));
  }

  // Check eigenvalues
  auto eigenvalues = method.getSingleParticleEnergies().getRestrictedEnergies();
  std::vector<double> expected = {-0.57457347, -0.33599573, -0.33599573, -0.33599573,
                                  0.35138661,  0.35138661,  0.35138661,  0.59576438};

  for (unsigned i = 0; i < eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i], DoubleNear(expected[i], 1e-6));
  }

  // Check force
  Eigen::RowVector3d f0(2.905661822261152E-017, 4.770489558936220E-017, 4.943961906533900E-017);
  Eigen::RowVector3d f1(-4.131279691482748E-005, -4.131279691484135E-005, -4.131279691485523E-005);
  Eigen::RowVector3d f2(4.131279691484829E-005, 4.131279691477197E-005, -4.131279691486217E-005);
  Eigen::RowVector3d f3(4.131279691477197E-005, -4.131279691478584E-005, 4.131279691480666E-005);
  Eigen::RowVector3d f4(-4.131279691484135E-005, 4.131279691482054E-005, 4.131279691486217E-005);

  ASSERT_TRUE((-method.getGradients().row(0) - (f0)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(1) - (f1)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(2) - (f2)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(3) - (f3)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(4) - (f4)).norm() < 1e-5);
}

TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForUnrestrictedCH3) {
  std::stringstream ss("4\n\n"
                       "C      0.0000000000    0.0000000000    0.0000000000\n"
                       "H      0.6287000000    0.6287000000    0.6287000000\n"
                       "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                       "H     -0.6287000000    0.6287000000   -0.6287000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);

  method.setScfMixer(scf_mixer_t::fock_diis);
  method.setMaxIterations(10000);
  method.setConvergenceCriteria({{}, 1e-8});

  method.setUnrestrictedCalculation(true);
  method.setSpinMultiplicity(2);
  method.setAtomCollection(as);
  method.initializeFromParameterPath("mio-1-1");
  method.calculate(Derivative::First, log);

  // Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(4));
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(7));

  // Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-2.7458925732, 1e-8));

  // Check atomic charges
  std::vector<double> expectedCharges = {-0.29151134, 0.09717045, 0.09717045, 0.09717045};
  for (unsigned long long i = 0; i < expectedCharges.size(); i++) {
    SCOPED_TRACE("... for the charge " + std::to_string(i) + ":");
    EXPECT_THAT(method.getAtomicCharges()[i], DoubleNear(expectedCharges[i], 1e-6));
  }

  // Check eigenvalues
  auto alphaEigenvalues = method.getSingleParticleEnergies().getAlphaEnergies();
  auto betaEigenvalues = method.getSingleParticleEnergies().getBetaEnergies();
  std::vector<double> alphaExpected = {-0.58373525, -0.35004087, -0.35004087, -0.23077050,
                                       0.33809319,  0.33809319,  0.51658701};
  std::vector<double> betaExpected = {-0.54684410, -0.32740014, -0.32740014, -0.18248896,
                                      0.36091443,  0.36091443,  0.53408108};

  for (unsigned i = 0; i < alphaEigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(alphaEigenvalues[i], DoubleNear(alphaExpected[i], 1e-6));
    EXPECT_THAT(betaEigenvalues[i], DoubleNear(betaExpected[i], 1e-6));
  }

  // Check force
  Eigen::RowVector3d f0(-1.738173073648139E-002, 1.738173073648139E-002, 1.738173073648149E-002);
  Eigen::RowVector3d f1(1.136823505565547E-002, -3.006747840412945E-003, -3.006747840413014E-003);
  Eigen::RowVector3d f2(3.006747840413014E-003, -1.136823505565546E-002, -3.006747840412952E-003);
  Eigen::RowVector3d f3(3.006747840412931E-003, -3.006747840412986E-003, -1.136823505565552E-002);

  ASSERT_TRUE((-method.getGradients().row(0) - (f0)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(1) - (f1)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(2) - (f2)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(3) - (f3)).norm() < 1e-5);
}

TEST_F(ADFTB2Calculation, MethodWrapperCanBeCloned) {
  dynamicallyLoadedMethodWrapper->settings().modifyInt(Utils::SettingsNames::molecularCharge, 2);

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

  auto mc = Utils::SettingsNames::molecularCharge;
  ASSERT_EQ(cloned->settings().getInt(mc), dynamicallyLoadedMethodWrapper->settings().getInt(mc));
}

TEST_F(ADFTB2Calculation, StructureIsCorrectlyCloned) {
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

TEST_F(ADFTB2Calculation, ClonedMethodCanCalculate) {
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

TEST_F(ADFTB2Calculation, ClonedMethodCanCalculateGradients) {
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
  ASSERT_EQ(resultCloned.get<Utils::Property::SuccessfulCalculation>(), true);
}

TEST_F(ADFTB2Calculation, ClonedMethodCanCalculateBigMoleculeGradient) {
  std::stringstream ss("60\n\n"
                       "C     14.1849078169   -1.8758387788   -3.5550946627\n"
                       "C     13.5788832264   -1.3458097616   -4.8525613198\n"
                       "C     12.2201246777   -0.6724215043   -4.6122599938\n"
                       "C     11.6166853967   -0.1491211064   -5.9220172230\n"
                       "C     10.2494348961    0.5098421841   -5.6891494939\n"
                       "C      9.6394138809    1.0225376340   -6.9581394994\n"
                       "C      8.4418328287    1.6090519467   -7.0548552919\n"
                       "C      7.4893515205    1.8507133434   -5.9271479434\n"
                       "C      7.3853093213    3.3520331769   -5.6036024374\n"
                       "H      8.3805338376    3.7733434390   -5.3324386180\n"
                       "O      6.5747924140    3.3779692042   -4.4040386912\n"
                       "C      6.7087257480    4.1384225509   -6.6940955188\n"
                       "C      7.3926281385    4.8147503515   -7.6250021328\n"
                       "C      6.7508745465    5.5662386384   -8.7045132125\n"
                       "C      7.3460483521    6.5380895687   -9.4117351917\n"
                       "C      8.7471184626    7.0215825541   -9.1836360219\n"
                       "C      9.2259616493    7.9204112461  -10.2854131882\n"
                       "C      9.5651515611    9.2063994777  -10.1471427231\n"
                       "C      9.4963176174    9.9732082613   -8.8610590274\n"
                       "C      9.9033903881   11.4433666463   -9.0243526014\n"
                       "C      9.8124381863   12.1925627082   -7.6866622426\n"
                       "C      9.7666220032   13.6822938370   -7.8383153092\n"
                       "O      9.1388063614   14.0263471971   -9.0254884324\n"
                       "C     10.2017589565   14.6205999547   -6.9803978496\n"
                       "C     10.8580578570   14.3143754852   -5.6793987128\n"
                       "H     15.1527529455   -2.3560760893   -3.7374830840\n"
                       "H     13.5345330729   -2.6190786174   -3.0799063515\n"
                       "H     14.3489848192   -1.0723891428   -2.8281060119\n"
                       "H     13.4642122909   -2.1758580921   -5.5771042776\n"
                       "H     14.2769421740   -0.6269930007   -5.3245847803\n"
                       "H     12.3354077301    0.1579019406   -3.8894450828\n"
                       "H     11.5254684942   -1.3903909635   -4.1352771604\n"
                       "H     11.5123834314   -0.9806650079   -6.6460878717\n"
                       "H     12.3086983872    0.5770265321   -6.3908553877\n"
                       "H     10.3515675699    1.3429928545   -4.9616762311\n"
                       "H      9.5603425964   -0.2172464456   -5.2088033701\n"
                       "H     10.2504963526    0.8878585283   -7.8537047821\n"
                       "H      8.0756046695    1.9610646161   -8.0227553175\n"
                       "H      7.7782223656    1.3128706282   -4.9999221583\n"
                       "H      6.4842658927    1.4526262172   -6.1820450047\n"
                       "H      6.3068846126    4.2899160126   -4.1889131896\n"
                       "H      5.6198586408    4.0981222010   -6.6733333412\n"
                       "H      8.4871508974    4.8262440428   -7.6209294400\n"
                       "H      5.7155517138    5.2879351423   -8.9204333391\n"
                       "H      6.8021161642    7.0550931985  -10.2061767928\n"
                       "H      9.4461128899    6.1581616289   -9.0971102155\n"
                       "H      8.7949869043    7.5271117175   -8.1936996888\n"
                       "H      9.2909208148    7.4310574737  -11.2591840949\n"
                       "H      9.9170272680    9.7750608861  -11.0108607213\n"
                       "H     10.1402890722    9.4856637214   -8.0995802307\n"
                       "H      8.4603605310    9.9194738690   -8.4596544694\n"
                       "H      9.2484335258   11.9354104576   -9.7756011134\n"
                       "H     10.9308547084   11.5204548740   -9.4277574248\n"
                       "H     10.6695620039   11.9109826223   -7.0381512230\n"
                       "H      8.8998258434   11.8769165418   -7.1334176671\n"
                       "H     11.5568320931   13.4697349434   -5.7506551269\n"
                       "H      9.0555975216   15.0032096536   -9.1323867763\n"
                       "H     10.0791176468   15.6781349435   -7.1983001948\n"
                       "H     10.1118372616   14.0635689781   -4.9111453842\n"
                       "H     11.4307176405   15.1746161628   -5.3076887006\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  dynamicallyLoadedMethodWrapper->setStructure(as);
  auto state = dynamicallyLoadedMethodWrapper->getState();

  int nThreads = omp_get_num_threads();
  omp_set_num_threads(1);
  dynamicallyLoadedMethodWrapper->setRequiredProperties(Utils::Property::Gradients);
  auto result = dynamicallyLoadedMethodWrapper->calculate("");

  omp_set_num_threads(2);
  auto cloned = dynamicallyLoadedMethodWrapper->clone();

  cloned->setRequiredProperties(Utils::Property::Gradients);

  cloned->loadState(state);
  auto resultCloned = cloned->calculate("");
  omp_set_num_threads(nThreads);

  ASSERT_EQ(resultCloned.get<Utils::Property::SuccessfulCalculation>(), true);
  for (int atom = 0; atom < cloned->getPositions().rows(); ++atom) {
    for (int dimension = 0; dimension < 3; ++dimension) {
      ASSERT_THAT(resultCloned.get<Utils::Property::Gradients>().row(atom)(dimension),
                  DoubleNear(result.get<Utils::Property::Gradients>().row(atom)(dimension), 1e-10));
    }
  }
}

TEST_F(ADFTB2Calculation, ClonedMethodCopiesResultsCorrectly) {
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

/*
TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForCO) {
  std::stringstream ss("2\n\n"
                         "C      0.0000000000    0.0000000000    0.0000000000\n"
                         "O      0.6287000000    0.6287000000    0.6287000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setAtomCollection(as);
  method.initialize("mio-1-1");
  method.calculate(1);

  //Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(2));
  ASSERT_THAT(method.getnAOs(), Eq(8));

  //Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-5.0281100118, 1e-10));

  //Check eigenvalues
  auto eigenvalues = method.getRestrictedEigenvalues();
  std::vector<double> expected = { -0.96215242, -0.50129687, -0.41598843, -0.41598843, -0.33713518, -0.01580125,
-0.01580125,  2.01165147 };


  for(long long i=0; i<eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i], DoubleNear(expected[i],1e-6));
  }

  //Check force
  Eigen::Vector3d f0 (-2.099482951227899E-002, -2.099482951228265E-002, -2.099482951227699E-002);
  Eigen::Vector3d f1 (2.099482951227899E-002,  2.099482951228265E-002,  2.099482951227699E-002);

  ASSERT_TRUE((method.getGradients()[0]-(-f0)).norm()<1e-4);
  ASSERT_TRUE((method.getGradients()[1]-(-f1)).norm()<1e-4);
}

TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForH2) {
  std::stringstream ss("2\n\n"
                         "H      0.0000000000    0.0000000000    0.0000000000\n"
                         "H      0.6287000000    0.6287000000    0.6287000000\n");
  auto as = Utils::XyzStreamHandler::read(ss);
  method.setAtomCollection(as);
  method.initialize("mio-1-1");
  method.calculate(1);

  //Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(2));
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(2));

  //Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-0.6389548092, 1e-10));

  //Check eigenvalues
  auto eigenvalues = method.getRestrictedEigenvalues();
  std::vector<double> expected = { -0.31947740, -0.04784076 };


  for(long long i=0; i<eigenvalues.size(); i++) {
    SCOPED_TRACE("... for the eigenvalue " + std::to_string(i) + ":");
    EXPECT_THAT(eigenvalues[i], DoubleNear(expected[i],1e-6));
  }

  //Check force
  Eigen::Vector3d f0 (3.390660659605971E-002,  3.390660659605971E-002,  3.390660659605971E-002);
  Eigen::Vector3d f1 (-3.390660659605971E-002, -3.390660659605971E-002, -3.390660659605971E-002);

  ASSERT_TRUE((method.getGradients()[0]-(-f0)).norm()<1e-5);
  ASSERT_TRUE((method.getGradients()[1]-(-f1)).norm()<1e-5);
}
*/

TEST_F(ADFTB2Calculation, AtomCollectionCanBeReturned) {
  std::stringstream ssH("4\n\n"
                        "C      0.0000000000    0.0000000000    0.0000000000\n"
                        "C      0.0529177211   -0.3175063264    0.2645886053\n"
                        "H     -0.5291772107    0.1058354421   -0.1587531632\n"
                        "H     -0.1058354421    0.1058354421   -0.1587531632\n");
  auto structure = Utils::XyzStreamHandler::read(ssH);
  calculator->setStructure(structure);
  ASSERT_EQ(structure.getPositions(), calculator->getStructure()->getPositions());
  for (unsigned int i = 0; i < structure.getElements().size(); ++i) {
    ASSERT_EQ(structure.getElements()[i], calculator->getStructure()->getElements()[i]);
  }
}

TEST_F(ADFTB2Calculation, CanCalculateAgWithDFTB2) {
  std::stringstream ssAg("1\n\n"
                         "Ag     0.000000   0.0000000   0.000000\n");

  calculator->settings().modifyInt("molecular_charge", 1);
  calculator->settings().modifyString("method_parameters", "hyb-0-2");
  calculator->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1.0e-8);
  calculator->setStructure(Utils::XyzStreamHandler::read(ssAg));
  calculator->setRequiredProperties(Utils::Property::Energy);
  calculator->calculate("");

  // From DFTB+
  EXPECT_THAT(calculator->results().get<Utils::Property::Energy>(), DoubleNear(-2.7971092750, 1e-7));
}

TEST_F(ADFTB2Calculation, CanPatchParametersWithDFTB2) {
  std::stringstream ssAg2O("3\n\n"
                           "Ag   -1.000000   0.0000000   0.000000\n"
                           "O     0.000000   0.0000000   0.000000\n"
                           "Ag    1.000000   0.0000000   0.000000");
  calculator->settings().modifyString("method_parameters", "hyb-0-2");
  calculator->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1.0e-8);
  calculator->setStructure(Utils::XyzStreamHandler::read(ssAg2O));
  calculator->setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);
  calculator->calculate("");

  // From DFTB+ (deviations consistent with other tests in this suite)
  EXPECT_THAT(calculator->results().get<Utils::Property::Energy>(), DoubleNear(-2.9939528717, 5e-6));

  // From DFTB+ (deviations consistent with other tests in this suite)
  Utils::GradientCollection referenceForces(3, 3);
  referenceForces << -9.841432449856, 0.000000000000, -0.000000000000, 0.000000000000, -0.000000000000, 0.000000000000,
      9.841432449856, -0.000000000000, 0.000000000000;
  for (int row = 0; row < 3; ++row)
    for (int col = 0; col < 3; ++col)
      EXPECT_THAT(calculator->results().get<Utils::Property::Gradients>()(row, col),
                  DoubleNear(-1.0 * referenceForces(row, col), 1e-5));
}

TEST_F(ADFTB2Calculation, GetsCorrectAtomicHessians) {
  std::stringstream ssH("5\n\n"
                        "C     -4.22875    2.29085   -0.00000\n"
                        "H     -3.42876    2.04248    0.66574\n"
                        "H     -3.90616    3.05518   -0.67574\n"
                        "H     -4.51405    1.42151   -0.55475\n"
                        "H     -5.06606    2.64424    0.56475\n");
  auto structure = Utils::XyzStreamHandler::read(ssH);
  method.setAtomCollection(structure);
  method.initializeFromParameterPath("mio-1-1");
  method.setConvergenceCriteria({1e-8, 1e-8});
  method.calculate(Derivative::SecondAtomic, log);
  auto ah = method.getAtomicSecondDerivatives();
  ASSERT_EQ(ah.size(), 5);

  // Reference data has been generated with original re-implementation in Sparrow
  // TODO: for some reason, these tests fail sometimes
  ASSERT_NEAR(ah.getAtomicHessian(0)(0, 0), 0.55435924542827597, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(0, 1), 5.2069606450091799e-06, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(0, 2), 1.9061722366005939e-05, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(1, 0), 5.2069606450091799e-06, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(1, 1), 0.55436564140093947, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(1, 2), -1.1246696738579331e-05, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(2, 0), 1.9061722366005939e-05, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(2, 1), -1.1246696738579331e-05, 1e-6);
  ASSERT_NEAR(ah.getAtomicHessian(0)(2, 2), 0.55436681809146315, 1e-6);
}

TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForCoComplex) {
  std::stringstream ss("8\n\n"
                       "Co    -0.0947841961    0.0542107769    0.2719105074\n"
                       "O     -2.8865621749   -0.3928440120   -0.3094733154\n"
                       "O      2.5152237470    0.6051735623    1.3747143456\n"
                       "C     -1.7828787564   -0.2201043822   -0.0964346185\n"
                       "O      0.8429878933   -0.4820973727   -2.4182169590\n"
                       "H     -0.5624976624    0.3217268493    1.6135135795\n"
                       "C      1.4939672661    0.3853190465    0.9252282712\n"
                       "C      0.4745438835   -0.2713844680   -1.3612418109\n");
  auto as = Utils::XyzStreamHandler::read(ss);

  method.setScfMixer(scf_mixer_t::fock_diis);
  method.setMaxIterations(10000);
  method.setConvergenceCriteria({{}, 1e-6});

  method.setAtomCollection(as);
  method.initializeFromParameterPath("trans3d-0-1");
  method.calculate(Derivative::First, log);

  // Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(8));
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(34));

  // Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-17.1512552737, 1e-5));

  // Check atomic charges
  std::vector<double> expectedCharges = {0.25059332,  -0.10229016, -0.10228875, 0.13467283,
                                         -0.11066338, -0.25545559, 0.13467202,  0.05075972};
  for (unsigned long long i = 0; i < expectedCharges.size(); i++) {
    SCOPED_TRACE("... for the charge " + std::to_string(i) + ":");
    EXPECT_THAT(method.getAtomicCharges()[i], DoubleNear(expectedCharges[i], 1e-4));
  }

  // Check force
  Eigen::RowVector3d f0(0.000245267388, -0.000139851946, -0.000703921747);
  Eigen::RowVector3d f1(0.000498579755, -0.000028853548, -0.000352593246);
  Eigen::RowVector3d f2(-0.000176331440, -0.000153352472, -0.000562721571);
  Eigen::RowVector3d f3(-0.052892558126, -0.007659772236, -0.007610079972);
  Eigen::RowVector3d f4(0.004003794557, -0.002289258671, -0.011485174294);
  Eigen::RowVector3d f5(-0.009424512040, 0.005390321540, 0.027033448971);
  Eigen::RowVector3d f6(0.047305523100, 0.010852332964, 0.023629559624);
  Eigen::RowVector3d f7(0.010440236806, -0.005971565632, -0.029948517764);

  ASSERT_TRUE((-method.getGradients().row(0) - (f0)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(1) - (f1)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(2) - (f2)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(3) - (f3)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(4) - (f4)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(5) - (f5)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(6) - (f6)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(7) - (f7)).norm() < 1e-5);
}

TEST_F(ADFTB2Calculation, GetsSameResultAsDFTBPlusForBH) {
  std::stringstream ss("2\n\n"
                       "B      0.0 0.0 0.0\n"
                       "H      1.0 0.0 0.0\n");
  auto as = Utils::XyzStreamHandler::read(ss);

  method.setScfMixer(scf_mixer_t::fock_diis);
  method.setMaxIterations(10000);
  method.setConvergenceCriteria({{}, 1e-6});

  method.setAtomCollection(as);
  method.initializeFromParameterPath("borg-0-1");
  method.calculate(Derivative::First, log);

  // Check number of atoms and orbitals
  ASSERT_THAT(method.getNumberAtoms(), Eq(2));
  ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(5));

  // Check energy
  ASSERT_THAT(method.getEnergy(), DoubleNear(-1.2439301670, 1e-5));

  // Check atomic charges
  std::vector<double> expectedCharges = {0.03249036, -0.03249036};

  for (unsigned long long i = 0; i < expectedCharges.size(); i++) {
    SCOPED_TRACE("... for the charge " + std::to_string(i) + ":");
    EXPECT_THAT(method.getAtomicCharges()[i], DoubleNear(expectedCharges[i], 1e-4));
  }

  // Check force
  Eigen::RowVector3d f0(-0.162415248703, 0.0, 0.0);
  Eigen::RowVector3d f1(0.162415248703, 0.0, 0.0);

  ASSERT_TRUE((-method.getGradients().row(0) - (f0)).norm() < 1e-5);
  ASSERT_TRUE((-method.getGradients().row(1) - (f1)).norm() < 1e-5);
}

TEST_F(ADFTB2Calculation, WorksForAllTrans3dCombinations) {
  std::stringstream ss("8\n\n"
                       "X     -0.0947841961    0.0542107769    0.2719105074\n"
                       "N     -2.8865621749   -0.3928440120   -0.3094733154\n"
                       "O      2.5152237470    0.6051735623    1.3747143456\n"
                       "C     -1.7828787564   -0.2201043822   -0.0964346185\n"
                       "O      0.8429878933   -0.4820973727   -2.4182169590\n"
                       "H     -0.5624976624    0.3217268493    1.6135135795\n"
                       "C      1.4939672661    0.3853190465    0.9252282712\n"
                       "C      0.4745438835   -0.2713844680   -1.3612418109\n");
  std::vector<std::string> elements = {"Fe", "Ni", "Sc", "Ti", "Co"};
  for (const auto& ele : elements) {
    auto s = ss.str();
    s.replace(3, 1, ele);
    if (ele == "Sc" || ele == "Co") {
      std::replace(s.begin(), s.end(), 'N', 'O');
    }
    std::stringstream mod(s);
    auto as = Utils::XyzStreamHandler::read(mod);

    method.setScfMixer(scf_mixer_t::fock_diis);
    method.setMaxIterations(10000);
    method.setConvergenceCriteria({{}, 1e-6});

    method.setAtomCollection(as);
    method.initializeFromParameterPath("trans3d-0-1");
    method.calculate(Derivative::First, log);

    // Check number of atoms and orbitals
    ASSERT_THAT(method.getNumberAtoms(), Eq(8));
    ASSERT_THAT(method.getNumberAtomicOrbitals(), Eq(34));

    // Check energy
    ASSERT_LT(method.getEnergy(), 0.0);
  }
}

} // namespace Sparrow
} // namespace Scine
