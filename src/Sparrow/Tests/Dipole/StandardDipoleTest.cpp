/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Sparrow/Implementations/DipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/AtomPairDipole.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/DipoleUtils/NDDODipoleMomentCalculator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/DataStructures/SlaterToGaussian.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <Utils/Settings.h>
#include <Utils/Typenames.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <gmock/gmock.h>
#include <Eigen/Core>
#include <array>
#include <chrono>

using namespace testing;

namespace Scine {
namespace Sparrow {

class SlaterToGaussianDipoleTest : public Test {
 public:
  std::shared_ptr<Core::Calculator> method;

  nddo::PM6Method underlyingMethod;
  std::unique_ptr<NDDODipoleMatrixCalculator<nddo::PM6Method>> matrixCalculator;
  std::unique_ptr<NDDODipoleMomentCalculator<nddo::PM6Method>> dipoleCalculator;

  double const arbitraryExponentA_{1.93834};
  double const arbitraryExponentB_{1.34521};
  Utils::GtoExpansion gtoHs, gtoFs, gtoFp;
  Utils::AtomicGtos AtomicGtosH, AtomicGtosF;

  Eigen::RowVector3d Ri, Rj, Ri2, Rj2;
  Eigen::Vector3d Rij;

  Utils::AtomCollection HF;
  Utils::AtomCollection Ethanol;

  Utils::DipoleMatrix dipoleMatrix_;
  Utils::DipoleMatrix dipoleMatrix2_;

  Eigen::Vector3d evaluationCoordinate;

  Core::Log log;

 private:
  void SetUp() override {
    log = Core::Log::silent();
    method = std::make_shared<PM6MethodWrapper>();
    method->setLog(log);

    underlyingMethod.setConvergenceCriteria({1e-7, 1e-9});

    auto& settings = method->settings();

    settings.modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-7);
    settings.modifyDouble(Utils::SettingsNames::densityRmsdCriterion, 1e-9);
    settings.modifyBool(Utils::SettingsNames::NDDODipoleApproximation, false);
    method->setRequiredProperties(Utils::Property::Energy);

    std::stringstream HFss("2\n\n"
                           "H        0.000000000     0.000000000     0.000000000\n"
                           "F        0.965548748     0.000000000     0.000000000\n");
    HF = Utils::XyzStreamHandler::read(HFss);
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
    Ethanol = Utils::XyzStreamHandler::read(Ethanolss);

    dipoleMatrix_.reset(5);
    dipoleMatrix2_.reset(5);
    gtoHs = Utils::SlaterToGaussian::getGTOExpansion(1, 1, 0, arbitraryExponentA_);
    gtoFs = Utils::SlaterToGaussian::getGTOExpansion(1, 2, 0, arbitraryExponentB_);
    gtoFp = Utils::SlaterToGaussian::getGTOExpansion(1, 2, 1, arbitraryExponentB_);

    AtomicGtosH.s = gtoHs;
    AtomicGtosF.s = gtoFs;
    AtomicGtosF.p = gtoFp;

    Ri << 0.0, 0.0, 0.0;
    Rj << 1.0, 0.0, 0.0;
    Ri2 << -0.5, 0.0, 0.0;
    Rj2 << 0.5, 0.0, 0.0;
    Rij << 1.0, 0.0, 0.0;
    evaluationCoordinate << 0.0, 0.0, 0.0;
  }
};

TEST_F(SlaterToGaussianDipoleTest, TestwithH2) {
  std::stringstream H2ss("2\n\n"
                         "H        0.000000000     0.000000000     0.000000000\n"
                         "H        1.000000000     0.000000000     0.000000000\n");
  auto structure = Utils::XyzStreamHandler::read(H2ss);
  method->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);
  method->setStructure(structure);

  auto res = method->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();

  ASSERT_THAT(dipole.norm(), DoubleNear(0.0, 0.001));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2sSTO_1G) {
  Utils::AtomicGtos AtomicGtosH, AtomicGtosF;
  AtomicGtosH.s = gtoHs;
  AtomicGtosF.s = gtoFs;
  AtomicGtosF.p = gtoFp;

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software
  // multiplying the GTO in each spatial coordinate obtained with the function Integrate[...]
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::DerivativeOrder::Zero>()(0, 1), DoubleNear(0.079595, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2pSTO_1G) {
  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software
  // multiplying the GTO in each spatial coordinate obtained with the function Integrate[...]
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::DerivativeOrder::Zero>()(0, 2), DoubleNear(0.1341609, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2sSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);
  gtoFs = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 0, arbitraryExponentB_);
  gtoFp = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 1, arbitraryExponentB_);

  AtomicGtosH.s = gtoHs;
  AtomicGtosF.s = gtoFs;
  AtomicGtosF.p = gtoFp;

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::DerivativeOrder::Zero>()(0, 1), DoubleNear(0.1288227, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2pSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);
  gtoFs = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 0, arbitraryExponentB_);
  gtoFp = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 1, arbitraryExponentB_);

  AtomicGtosH.s = gtoHs;
  AtomicGtosF.s = gtoFs;
  AtomicGtosF.p = gtoFp;

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::DerivativeOrder::Zero>()(0, 2), DoubleNear(0.1815766, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestyDipoleMatrixBlock1s2pSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);
  gtoFs = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 0, arbitraryExponentB_);
  gtoFp = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 1, arbitraryExponentB_);

  AtomicGtosH.s = gtoHs;
  AtomicGtosF.s = gtoFs;
  AtomicGtosF.p = gtoFp;

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.y().get<Utils::DerivativeOrder::Zero>()(0, 3), DoubleNear(0.3524788, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s1sSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);

  AtomicGtosH.s = gtoHs;

  Rj.x() = 0.5;
  Rij = Eigen::Vector3d::Zero(3);
  auto const firstAOIndexOnH = 0;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnH, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosH, Rj, Rj, Rij, evaluationCoordinate);
  // Diagonal elements of the dipole matrix should be equal to the distance from the evaluation point
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::DerivativeOrder::Zero>()(0, 0), DoubleNear(Rj.x(), 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixHH) {
  std::stringstream HH("2\n\n"
                       "H 0.0 0.0 0.0\n"
                       "H 1.0 0.0 0.0\n");

  auto HHst = Utils::XyzStreamHandler::read(HH);
  method->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixAO);
  method->setStructure(HHst);
  auto result = method->calculate("");
  const auto dipoleMatrix = result.take<Utils::Property::DipoleMatrixAO>();
  // The analytical values were computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix.x().get<Utils::DerivativeOrder::Zero>()(0, 0), DoubleNear(0, 1e-6));
  ASSERT_THAT(dipoleMatrix.x().get<Utils::DerivativeOrder::Zero>()(0, 1), DoubleNear(0.45662245, 1e-5));
  ASSERT_THAT(dipoleMatrix.x().get<Utils::DerivativeOrder::Zero>()(1, 1), DoubleNear(1.8897346, 1e-5));
}

TEST_F(SlaterToGaussianDipoleTest, ObaraSaikaAndAnalyticalAreEqual) {
  method->setStructure(Ethanol);
  method->calculate("");

  underlyingMethod.setStructure(Ethanol);
  underlyingMethod.convergedCalculation(log, Utils::Derivative::None);

  matrixCalculator = NDDODipoleMatrixCalculator<nddo::PM6Method>::create(underlyingMethod);
  dipoleCalculator = NDDODipoleMomentCalculator<nddo::PM6Method>::create(underlyingMethod, *matrixCalculator);

  matrixCalculator->setIntegralMethod(IntegralMethod::ObaraSaika);
  matrixCalculator->fillDipoleMatrix(evaluationCoordinate);
  const auto obaraSaikaDM = matrixCalculator->getAODipoleMatrix();

  matrixCalculator->setIntegralMethod(IntegralMethod::ClosedForm);
  matrixCalculator->fillDipoleMatrix(evaluationCoordinate);
  const auto closedFormDM = matrixCalculator->getAODipoleMatrix();

  for (int dimension = 0; dimension < 3; ++dimension) {
    for (int i = 0; i < obaraSaikaDM.x().get<Utils::DerivativeOrder::Zero>().cols(); ++i) {
      for (int j = 0; j < obaraSaikaDM.x().get<Utils::DerivativeOrder::Zero>().cols(); ++j) {
        ASSERT_NEAR(obaraSaikaDM[dimension](i, j), closedFormDM[dimension](i, j), 5e-6);
      }
    }
  }
}

TEST_F(SlaterToGaussianDipoleTest, NuclearDipoleContribution) {
  // pyscf nuclear dipole = [42.0259269  15.37531606 -1.03465377]
  /*
   * Core charges in pySCF were assigned to the PM6 core charges (C: 4, O: 6, H: 1).
   * The nuclear dipole contribution was calculated separately from the electronic dipole contribution
   */
  underlyingMethod.setStructure(Ethanol);
  underlyingMethod.convergedCalculation(log, Utils::Derivative::None);
  evaluationCoordinate.setZero(3);

  // set dm to 0 and calculate dipole -> nuclear dipole!
  auto const& dm = underlyingMethod.getDensityMatrix();
  auto nEle = dm.numberElectrons();
  auto dim = dm.restrictedMatrix().rows();
  Eigen::MatrixXd newDM = Eigen::MatrixXd::Zero(dim, dim);
  Utils::DensityMatrix emptyDensity;
  emptyDensity.setDensity(std::move(newDM), nEle);
  underlyingMethod.setDensityMatrix(std::move(emptyDensity));

  matrixCalculator = NDDODipoleMatrixCalculator<nddo::PM6Method>::create(underlyingMethod);
  dipoleCalculator = NDDODipoleMomentCalculator<nddo::PM6Method>::create(underlyingMethod, *matrixCalculator);

  dipoleCalculator->useNDDOApproximation(false);
  auto nuclearDipole = dipoleCalculator->calculate();

  ASSERT_THAT(nuclearDipole.x(), DoubleNear(42.0259269, 1e-6));
  ASSERT_THAT(nuclearDipole.y(), DoubleNear(15.37531606, 1e-6));
  ASSERT_THAT(nuclearDipole.z(), DoubleNear(-1.03465377, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, DipoleMatrixEthanolEqualsLibCInt) {
  /*
   * Dipole calculated with LibCInt through pySCF, using the basis of scine, and omitting 1s core orbitals.
   */
  evaluationCoordinate.setZero(3);

  method->setStructure(Ethanol);
  method->setRequiredProperties(Utils::Property::Energy | Utils::Property::DipoleMatrixAO);

  auto results = method->calculate("");

  auto dipoleMatrix = results.take<Utils::Property::DipoleMatrixAO>();

  std::array<std::array<double, 21>, 21> pySCFAOInt = {
      {{{-4.86089394e-02, 7.53562584e-01, -1.26305817e-19, 4.36248881e-20, 2.76487387e-01, -2.58404714e-01,
         3.95494062e-19, -1.96837103e-19, 6.34714563e-03, -3.53257237e-02, -2.96113894e-02, 3.58378269e-20,
         -1.65300944e-01, -1.72549234e-01, -1.69956587e-01, 9.04051444e-02, 8.95789995e-02, 3.01365535e-02}},
       {{7.53562584e-01, -4.86089394e-02, -1.70210969e-18, 6.46638804e-19, 6.26689285e-01, -4.49169782e-01,
         -4.01736305e-19, 5.54291057e-19, 2.56778547e-02, -7.85260918e-02, -8.59181194e-02, -6.74047316e-20,
         4.15444918e-01, 4.23431296e-01, 4.20763782e-01, 2.41168312e-01, 2.39430213e-01, 8.03998098e-02}},
       {{-1.26305817e-19, -1.70210969e-18, -4.86089394e-02, 6.67072375e-37, 7.09559882e-19, 2.28318526e-19,
         2.42307395e-01, -3.88719403e-37, 1.65233713e-02, -5.85754328e-02, -2.36510962e-02, -3.88894736e-20,
         -1.02751624e-01, -9.14936138e-02, 1.94101080e-01, -4.97347900e-02, -3.76785695e-02, 5.45758247e-02}},
       {{4.36248881e-20, 6.46638804e-19, 6.67072375e-37, -4.86089394e-02, -2.04309749e-19, -1.77077517e-19,
         -3.88719403e-37, 2.42307395e-01, 2.55072135e-21, 2.27381825e-20, 2.01172748e-20, 2.61013402e-02,
         -1.60046200e-01, 1.74203123e-01, -9.96820766e-03, -7.90782227e-02, 8.45094228e-02, -2.35020500e-02}},
       {{2.76487387e-01, 6.26689285e-01, 7.09559882e-19, -2.04309749e-19, 2.82490054e+00, 7.53562584e-01,
         7.34024185e-18, -2.53525321e-18, 1.68710002e-01, -1.63500022e-01, -7.17721068e-01, 7.18204812e-19,
         8.99455998e-02, 8.73886717e-02, 8.84281746e-02, 1.32581449e+00, 1.32343363e+00, 2.73390212e-01}},
       {{-2.58404714e-01, -4.49169782e-01, 2.28318526e-19, -1.77077517e-19, 7.53562584e-01, 2.82490054e+00,
         -1.70210969e-18, 6.46638804e-19, 1.61659371e-01, 2.83718718e-01, -6.07155385e-01, -2.66615416e-19,
         -4.58742552e-02, -4.41017017e-02, -4.48220022e-02, 9.06564051e-01, 9.13619131e-01, 1.41382051e-01}},
       {{3.95494062e-19, -4.01736305e-19, 2.42307395e-01, -3.88719403e-37, 7.34024185e-18, -1.70210969e-18,
         2.82490054e+00, -3.87668016e-35, 3.89597697e-01, -2.79412268e-01, -7.77873141e-01, -4.96704414e-19,
         3.26259270e-02, 2.68174223e-02, -5.86290495e-02, -7.45055510e-01, -5.68586850e-01, 4.31996374e-01}},
       {{-1.96837103e-19, 5.54291057e-19, -3.88719403e-37, 2.42307395e-01, -2.53525321e-18, 6.46638804e-19,
         -3.87668016e-35, 2.82490054e+00, 4.91740585e-20, 1.22157230e-19, 4.98258320e-19, 4.43638461e-01,
         5.08182298e-02, -5.10601615e-02, 3.01093914e-03, -1.18463686e+00, 1.27528585e+00, -1.86031095e-01}},
       {{6.34714563e-03, 2.56778547e-02, 1.65233713e-02, 2.55072135e-21, 1.68710002e-01, 1.61659371e-01, 3.89597697e-01,
         4.91740585e-20, 3.81247002e+00, 2.37061798e-01, 1.84654850e-16, -1.51047374e-18, 9.79171040e-03,
         8.67466395e-03, 1.86001827e-03, 4.56106976e-02, 5.57270463e-02, 5.70182494e-01}},
       {{-3.53257237e-02, -7.85260918e-02, -5.85754328e-02, 2.27381825e-20, -1.63500022e-01, 2.83718718e-01,
         -2.79412268e-01, 1.22157230e-19, 2.37061798e-01, 3.81247002e+00, 3.12398005e-17, 3.69594354e-19,
         -2.88790437e-02, -2.55160185e-02, -5.22265776e-03, 2.86524976e-02, 3.42451594e-02, -1.36980893e-01}},
       {{-2.96113894e-02, -8.59181194e-02, -2.36510962e-02, 2.01172748e-20, -7.17721068e-01, -6.07155385e-01,
         -7.77873141e-01, 4.98258320e-19, 1.84654850e-16, 3.12398005e-17, 3.81247002e+00, 0.00000000e+00,
         -1.30427566e-02, -1.26913084e-02, -7.49314084e-03, -2.32740782e-01, -2.66291146e-01, 6.89106286e-01}},
       {{3.58378269e-20, -6.74047316e-20, -3.88894736e-20, 2.61013402e-02, 7.18204812e-19, -2.66615416e-19,
         -4.96704414e-19, 4.43638461e-01, -1.51047374e-18, 3.69594354e-19, 0.00000000e+00, 3.81247002e+00,
         1.43881726e-02, -1.32001693e-02, 1.66308281e-04, -1.07162017e-01, 1.41953852e-01, -1.10649959e+00}},
       {{-1.65300944e-01, 4.15444918e-01, -1.02751624e-01, -1.60046200e-01, 8.99455998e-02, -4.58742552e-02,
         3.26259270e-02, 5.08182298e-02, 9.79171040e-03, -2.88790437e-02, -1.30427566e-02, 1.43881726e-02,
         -7.90831229e-01, -1.31321587e-01, -1.29433544e-01, 2.22812347e-02, 6.35440804e-02, 2.41879386e-02}},
       {{-1.72549234e-01, 4.23431296e-01, -9.14936138e-02, 1.74203123e-01, 8.73886717e-02, -4.41017017e-02,
         2.68174223e-02, -5.10601615e-02, 8.67466395e-03, -2.55160185e-02, -1.26913084e-02, -1.32001693e-02,
         -1.31321587e-01, -8.23762532e-01, -1.32303701e-01, 5.95958955e-02, 2.16755964e-02, 5.28179934e-02}},
       {{-1.69956587e-01, 4.20763782e-01, 1.94101080e-01, -9.96820766e-03, 8.84281746e-02, -4.48220022e-02,
         -5.86290495e-02, 3.01093914e-03, 1.86001827e-03, -5.22265776e-03, -7.49314084e-03, 1.66308281e-04,
         -1.29433544e-01, -1.32303701e-01, -8.11084801e-01, 5.98649707e-02, 5.79995402e-02, 6.53109527e-03}},
       {{9.04051444e-02, 2.41168312e-01, -4.97347900e-02, -7.90782227e-02, 1.32581449e+00, 9.06564051e-01,
         -7.45055510e-01, -1.18463686e+00, 4.56106976e-02, 2.86524976e-02, -2.32740782e-01, -1.07162017e-01,
         2.22812347e-02, 5.95958955e-02, 5.98649707e-02, 3.60800344e+00, 5.50665860e-01, 1.99365731e-01}},
       {{8.95789995e-02, 2.39430213e-01, -3.76785695e-02, 8.45094228e-02, 1.32343363e+00, 9.13619131e-01,
         -5.68586850e-01, 1.27528585e+00, 5.57270463e-02, 3.42451594e-02, -2.66291146e-01, 1.41953852e-01,
         6.35440804e-02, 2.16755964e-02, 5.79995402e-02, 5.50665860e-01, 3.62105951e+00, 8.87843454e-02}},
       {{3.01365535e-02, 8.03998098e-02, 5.45758247e-02, -2.35020500e-02, 2.73390212e-01, 1.41382051e-01,
         4.31996374e-01, -1.86031095e-01, 5.70182494e-01, -1.36980893e-01, 6.89106286e-01, -1.10649959e+00,
         2.41879386e-02, 5.28179934e-02, 6.53109527e-03, 1.99365731e-01, 8.87843454e-02, 3.24255598e+00}}}};

  for (int mu = 0; mu < 18; ++mu) {
    for (int nu = mu; nu < 18; ++nu) {
      ASSERT_THAT(dipoleMatrix.x().get<Utils::DerivativeOrder::Zero>()(mu, nu), DoubleNear(pySCFAOInt[mu][nu], 5e-6));
    }
  }
}

TEST_F(SlaterToGaussianDipoleTest, DiHydrogenDipoleTestStretching) {
  std::stringstream H21ss("2\n\n"
                          "H        0.000000000     0.000000000     0.000000000\n"
                          "H        1.000000000     0.000000000     0.000000000\n");

  method->settings().modifyString(Utils::SettingsNames::spinMode,
                                  Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  method->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 3);
  Utils::AtomCollection H21A = Utils::XyzStreamHandler::read(H21ss);
  Utils::DisplacementCollection Disp(2, 3);
  auto positions = H21A.getPositions();

  Disp.row(0) << 0, 0, 0;
  Disp.row(1) << 1, 0, 0;

  method->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);

  method->setStructure(H21A);
  method->settings().modifyString(Utils::SettingsNames::spinMode,
                                  Utils::SpinModeInterpreter::getStringFromSpinMode(Utils::SpinMode::Unrestricted));
  method->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 3);
  auto res = method->calculate("");
  auto dipole1 = res.get<Utils::Property::Dipole>();

  method->modifyPositions(positions + Disp);
  auto res2 = method->calculate("");
  auto dipole2 = res2.get<Utils::Property::Dipole>();

  Disp.row(1) *= 2;

  method->modifyPositions(positions + Disp);
  res = method->calculate("");
  auto dipole3 = res.get<Utils::Property::Dipole>();

  Disp.row(1) *= 3;
  method->modifyPositions(positions + Disp);
  res = method->calculate("");
  auto dipole4 = res.get<Utils::Property::Dipole>();

  Disp.row(1) *= 5;
  method->modifyPositions(positions + Disp);
  res = method->calculate("");
  auto dipole6 = res.get<Utils::Property::Dipole>();

  Disp.row(1) *= 9;
  method->modifyPositions(positions + Disp);
  res = method->calculate("");
  auto dipole10 = res.get<Utils::Property::Dipole>();

  Disp.row(1) *= 99;
  method->modifyPositions(positions + Disp);
  res = method->calculate("");
  auto dipole100 = res.get<Utils::Property::Dipole>();

  ASSERT_NEAR(dipole1.norm(), 0., 1e-5);
  ASSERT_NEAR(dipole2.norm(), 0., 1e-5);
  ASSERT_NEAR(dipole3.norm(), 0., 1e-5);
  ASSERT_NEAR(dipole4.norm(), 0., 1e-5);
  ASSERT_NEAR(dipole6.norm(), 0., 1e-5);
  ASSERT_NEAR(dipole10.norm(), 0., 1e-5);
  ASSERT_NEAR(dipole100.norm(), 0., 1e-5);
}

TEST_F(SlaterToGaussianDipoleTest, EthanolDipoleTest) {
  method->setStructure(Ethanol);
  method->setRequiredProperties(Utils::Property::Dipole | Utils::Property::Energy);
  auto res = method->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();
  // Value computed with pySCF for the same structure, with the Hartree-Fock Method and the STO-6G basis.
  ASSERT_NEAR(dipole.norm() * 2.541765, 1.54, 6e-2);
}

TEST_F(SlaterToGaussianDipoleTest, HFDipoleTest) {
  method->setStructure(HF);
  method->setRequiredProperties(Utils::Property::Dipole | Utils::Property::Energy);
  auto res = method->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();

  std::stringstream HFRotatedss("2\n\n"
                                "H        0.000000000     0.000000000     0.000000000\n"
                                "F        0.000000000     0.965548748     0.000000000\n");
  Utils::AtomCollection HF2 = Utils::XyzStreamHandler::read(HFRotatedss);
  method->setStructure(HF2);
  auto res2 = method->calculate("");
  auto const& dipole2 = res2.get<Utils::Property::Dipole>();

  ASSERT_THAT(dipole2.norm(), DoubleNear(dipole.norm(), 1e-6));

  std::stringstream HFRotatedTranslatedss("2\n\n"
                                          "H        0.000000000     0.000000000     5.000000000\n"
                                          "F        0.000000000     0.965548748     5.000000000\n");
  Utils::AtomCollection HF3 = Utils::XyzStreamHandler::read(HFRotatedTranslatedss);
  method->setStructure(HF3);
  auto res3 = method->calculate("");
  auto const& dipole3 = res3.get<Utils::Property::Dipole>();

  ASSERT_THAT(dipole3.norm(), DoubleNear(dipole.norm(), 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestCO2Dipole) {
  std::stringstream CO2("3\n\n"
                        "C 0.0 0.0 0.0\n"
                        "O 0.529177210 0.0 0.0\n"
                        "O -0.529177210 0.0 0.0\n");

  auto CO2st = Utils::XyzStreamHandler::read(CO2);
  method->setStructure(CO2st);
  method->setRequiredProperties(Utils::Property::Dipole | Utils::Property::Energy);
  auto res = method->calculate("");
  auto const& dipole = res.get<Utils::Property::Dipole>();

  // dipole should be 0 due to symmetry.
  ASSERT_NEAR(dipole.norm(), 0.0, 1e-6);
}

} // namespace Sparrow
} // namespace Scine
