/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "../MethodsTests/parameters_location.h"
#include <Core/Interfaces/Calculator.h>
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

 private:
  void SetUp() override {
    method = std::make_shared<PM6MethodWrapper>();

    underlyingMethod.setConvergenceCriteria(1e-9);

    auto& settings = method->settings();
    settings.modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-9);
    settings.modifyString(Utils::SettingsNames::parameterRootDirectory, parameters_root);
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

    AtomicGtosH.setS(gtoHs);
    AtomicGtosF.setS(gtoFs);
    AtomicGtosF.setP(gtoFp);

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
  AtomicGtosH.setS(gtoHs);
  AtomicGtosF.setS(gtoFs);
  AtomicGtosF.setP(gtoFp);

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software
  // multiplying the GTO in each spatial coordinate obtained with the function Integrate[...]
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::derivOrder::zero>()(0, 1), DoubleNear(0.079595, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2pSTO_1G) {
  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software
  // multiplying the GTO in each spatial coordinate obtained with the function Integrate[...]
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::derivOrder::zero>()(0, 2), DoubleNear(0.1341609, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2sSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);
  gtoFs = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 0, arbitraryExponentB_);
  gtoFp = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 1, arbitraryExponentB_);

  AtomicGtosH.setS(gtoHs);
  AtomicGtosF.setS(gtoFs);
  AtomicGtosF.setP(gtoFp);

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::derivOrder::zero>()(0, 1), DoubleNear(0.1288227, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s2pSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);
  gtoFs = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 0, arbitraryExponentB_);
  gtoFp = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 1, arbitraryExponentB_);

  AtomicGtosH.setS(gtoHs);
  AtomicGtosF.setS(gtoFs);
  AtomicGtosF.setP(gtoFp);

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::derivOrder::zero>()(0, 2), DoubleNear(0.1815766, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestyDipoleMatrixBlock1s2pSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);
  gtoFs = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 0, arbitraryExponentB_);
  gtoFp = Utils::SlaterToGaussian::getGTOExpansion(2, 2, 1, arbitraryExponentB_);

  AtomicGtosH.setS(gtoHs);
  AtomicGtosF.setS(gtoFs);
  AtomicGtosF.setP(gtoFp);

  auto const firstAOIndexOnH = 0;
  auto const firstAOIndexOnF = 1;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnF, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosF, Ri, Rj, Rij, evaluationCoordinate);
  // The analytical value was computed with the Mathematica software.
  // The integral between GTO pairs were calculated by summing the multiplying the integrals in the three
  // spatial coordinates. The normalized contributions are then summed together
  // (See Molecular Electronic Structure Theory, Trygve Helgaker, Chapter 9)
  ASSERT_THAT(dipoleMatrix_.y().get<Utils::derivOrder::zero>()(0, 3), DoubleNear(0.3524788, 1e-6));
}

TEST_F(SlaterToGaussianDipoleTest, TestxDipoleMatrixBlock1s1sSTO_2G) {
  gtoHs = Utils::SlaterToGaussian::getGTOExpansion(2, 1, 0, arbitraryExponentA_);

  AtomicGtosH.setS(gtoHs);

  Rj.x() = 0.5;
  Rij = Eigen::Vector3d::Zero(3);
  auto const firstAOIndexOnH = 0;

  AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndexOnH, firstAOIndexOnH, IntegralMethod::ObaraSaika,
                                          AtomicGtosH, AtomicGtosH, Rj, Rj, Rij, evaluationCoordinate);
  // Diagonal elements of the dipole matrix should be equal to the distance from the evaluation point
  ASSERT_THAT(dipoleMatrix_.x().get<Utils::derivOrder::zero>()(0, 0), DoubleNear(Rj.x(), 1e-6));
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
  ASSERT_THAT(dipoleMatrix.x().get<Utils::derivOrder::zero>()(0, 0), DoubleNear(0, 1e-6));
  ASSERT_THAT(dipoleMatrix.x().get<Utils::derivOrder::zero>()(0, 1), DoubleNear(0.45662245, 1e-5));
  ASSERT_THAT(dipoleMatrix.x().get<Utils::derivOrder::zero>()(1, 1), DoubleNear(1.8897346, 1e-5));
}

TEST_F(SlaterToGaussianDipoleTest, ObaraSaikaAndAnalyticalAreEqual) {
  method->setStructure(Ethanol);
  method->calculate("");

  underlyingMethod.setStructure(Ethanol, parameters_pm6);
  underlyingMethod.convergedCalculation(Utils::derivativeType::none);

  matrixCalculator = NDDODipoleMatrixCalculator<nddo::PM6Method>::create(underlyingMethod);
  dipoleCalculator = NDDODipoleMomentCalculator<nddo::PM6Method>::create(underlyingMethod, *matrixCalculator);

  matrixCalculator->setIntegralMethod(IntegralMethod::ObaraSaika);
  matrixCalculator->fillDipoleMatrix(evaluationCoordinate);
  const auto obaraSaikaDM = matrixCalculator->getAODipoleMatrix();

  matrixCalculator->setIntegralMethod(IntegralMethod::ClosedForm);
  matrixCalculator->fillDipoleMatrix(evaluationCoordinate);
  const auto closedFormDM = matrixCalculator->getAODipoleMatrix();

  for (int dimension = 0; dimension < 3; ++dimension) {
    for (int i = 0; i < obaraSaikaDM.x().get<Utils::derivOrder::zero>().cols(); ++i) {
      for (int j = 0; j < obaraSaikaDM.x().get<Utils::derivOrder::zero>().cols(); ++j) {
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
  underlyingMethod.setStructure(Ethanol, parameters_pm6);
  underlyingMethod.convergedCalculation(Utils::derivativeType::none);
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
      ASSERT_THAT(dipoleMatrix.x().get<Utils::derivOrder::zero>()(mu, nu), DoubleNear(pySCFAOInt[mu][nu], 5e-6));
    }
  }
}

TEST_F(SlaterToGaussianDipoleTest, DiHydrogenDipoleTestStretching) {
  std::stringstream H21ss("2\n\n"
                          "H        0.000000000     0.000000000     0.000000000\n"
                          "H        1.000000000     0.000000000     0.000000000\n");

  method->settings().modifyBool(Utils::SettingsNames::unrestrictedCalculation, true);
  method->settings().modifyInt(Utils::SettingsNames::spinMultiplicity, 3);
  Utils::AtomCollection H21A = Utils::XyzStreamHandler::read(H21ss);
  Utils::DisplacementCollection Disp(2, 3);
  auto positions = H21A.getPositions();

  Disp.row(0) << 0, 0, 0;
  Disp.row(1) << 1, 0, 0;

  method->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);

  method->setStructure(H21A);
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

TEST_F(SlaterToGaussianDipoleTest, DipoleIsCorrectForBigOrganicMolecule) {
  method->setRequiredProperties(Utils::Property::Energy | Utils::Property::Dipole);
  method->settings().modifyDouble(Utils::SettingsNames::selfConsistanceCriterion, 1e-5);
  std::stringstream BigMoleculeSS("113\n\n"
                                  "C    14.663000    29.225000     8.550000\n"
                                  "O    15.679000    23.530000     9.284000\n"
                                  "C    15.446000    24.383000     8.168000\n"
                                  "N    14.572000    30.082000    11.272000\n"
                                  "O    15.782000    29.646000     8.475000\n"
                                  "C    13.465000    29.993000     9.084000\n"
                                  "O    13.409000    23.131000     8.031000\n"
                                  "C    14.530000    23.550000     7.207000\n"
                                  "O    13.665000    31.396000     8.952000\n"
                                  "C    13.981000    24.239000     5.913000\n"
                                  "O    11.843000    25.089000     6.726000\n"
                                  "C    13.331000    29.725000    10.608000\n"
                                  "O    11.300000    22.870000     4.513000\n"
                                  "C    12.443000    24.165000     5.775000\n"
                                  "O    15.297000    27.944000    11.449000\n"
                                  "C    11.809000    24.239000     4.371000\n"
                                  "O    14.900000    24.579000     2.267000\n"
                                  "C    15.500000    29.146000    11.631000\n"
                                  "O    16.937000    23.198000     4.018000\n"
                                  "C    12.690000    24.440000     3.137000\n"
                                  "O    18.173000    25.452000     4.631000\n"
                                  "C    14.142000    24.761000     3.477000\n"
                                  "O    14.293000    27.968000     8.242000\n"
                                  "C    14.685000    23.777000     4.594000\n"
                                  "O    14.180000    21.013000     8.016000\n"
                                  "C    16.211000    24.019000     4.536000\n"
                                  "O    12.296000    26.819000     5.358000\n"
                                  "C    16.768000    25.373000     4.958000\n"
                                  "O    17.590000    25.630000     2.474000\n"
                                  "C    16.537000    25.716000     6.403000\n"
                                  "C    15.883000    26.853000     6.702000\n"
                                  "C    15.365000    26.987000     8.120000\n"
                                  "C    14.741000    25.692000     8.628000\n"
                                  "C    16.829000    24.694000     7.509000\n"
                                  "C    17.817000    25.283000     8.563000\n"
                                  "C    17.533000    23.387000     7.066000\n"
                                  "C    14.452000    22.283000     4.272000\n"
                                  "C    15.595000    27.978000     5.756000\n"
                                  "C    11.716000    22.817000     5.876000\n"
                                  "C    13.407000    21.834000     8.432000\n"
                                  "C    12.322000    21.565000     9.402000\n"
                                  "C    11.326000    22.468000     9.688000\n"
                                  "C    10.344000    22.126000    10.598000\n"
                                  "C    10.368000    20.934000    11.245000\n"
                                  "C    11.346000    20.041000    10.963000\n"
                                  "C    12.328000    20.345000    10.030000\n"
                                  "C    12.155000    30.467000    11.228000\n"
                                  "C    10.972000    30.594000    10.576000\n"
                                  "C     9.881000    31.266000    11.159000\n"
                                  "C    10.016000    31.826000    12.365000\n"
                                  "C    11.179000    31.735000    12.995000\n"
                                  "C    12.278000    31.046000    12.470000\n"
                                  "C    11.819000    26.403000     6.384000\n"
                                  "C    11.141000    27.185000     7.455000\n"
                                  "C    16.728000    29.689000    12.287000\n"
                                  "C    16.692000    30.723000    13.191000\n"
                                  "C    17.815000    31.137000    13.839000\n"
                                  "C    19.009000    30.546000    13.564000\n"
                                  "C    19.099000    29.524000    12.617000\n"
                                  "C    17.936000    29.083000    11.990000\n"
                                  "C    18.443000    25.606000     3.318000\n"
                                  "C    19.908000    25.754000     3.065000\n"
                                  "H    15.164000    23.911000    10.056000\n"
                                  "H    12.660000    29.703000     8.594000\n"
                                  "H    15.064000    22.834000     6.818000\n"
                                  "H    14.210000    25.180000     6.039000\n"
                                  "H    13.157000    28.794000    10.711000\n"
                                  "H    14.747000    31.020000    11.471000\n"
                                  "H    11.229000    24.988000     4.178000\n"
                                  "H    14.210000    25.658000     3.808000\n"
                                  "H    15.700000    25.275000     2.308000\n"
                                  "H    13.743000    31.792000     9.823000\n"
                                  "H    16.297000    26.065000     4.473000\n"
                                  "H    16.158000    27.238000     8.613000\n"
                                  "H    11.328000    23.337000     9.259000\n"
                                  "H     9.619000    22.738000    10.768000\n"
                                  "H     9.659000    20.728000    11.879000\n"
                                  "H    11.368000    19.196000    11.414000\n"
                                  "H    13.017000    19.699000     9.847000\n"
                                  "H    10.871000    30.182000     9.695000\n"
                                  "H     9.023000    31.355000    10.702000\n"
                                  "H     9.281000    32.288000    12.791000\n"
                                  "H    11.288000    32.169000    13.902000\n"
                                  "H    13.137000    31.044000    12.962000\n"
                                  "H    15.840000    31.211000    13.408000\n"
                                  "H    17.787000    31.834000    14.500000\n"
                                  "H    19.854000    30.900000    13.931000\n"
                                  "H    20.013000    29.129000    12.411000\n"
                                  "H    18.006000    28.339000    11.329000\n"
                                  "H    12.322000    25.156000     2.592000\n"
                                  "H    12.660000    23.648000     2.602000\n"
                                  "H    14.687000    25.706000     9.581000\n"
                                  "H    13.812000    25.658000     8.300000\n"
                                  "H    18.662000    25.467000     8.167000\n"
                                  "H    17.946000    24.653000     9.278000\n"
                                  "H    18.383000    23.576000     6.647000\n"
                                  "H    17.668000    22.810000     7.806000\n"
                                  "H    16.992000    22.906000     6.419000\n"
                                  "H    14.826000    21.733000     4.985000\n"
                                  "H    13.534000    22.092000     4.178000\n"
                                  "H    14.925000    22.044000     3.466000\n"
                                  "H    15.979000    27.765000     4.890000\n"
                                  "H    14.667000    28.076000     5.650000\n"
                                  "H    15.979000    28.770000     6.058000\n"
                                  "H    10.990000    22.810000     6.524000\n"
                                  "H    12.282000    22.068000     6.106000\n"
                                  "H    11.090000    28.124000     7.264000\n"
                                  "H    10.235000    26.879000     7.635000\n"
                                  "H    11.606000    27.094000     8.318000\n"
                                  "H    17.469000    26.089000     8.926000\n"
                                  "H    20.093000    25.850000     2.137000\n"
                                  "H    20.252000    26.520000     3.533000\n"
                                  "H    20.371000    24.964000     3.381000\n");

  auto BigMolecule = Utils::XyzStreamHandler::read(BigMoleculeSS);
  method->setStructure(BigMolecule);

  auto res = method->calculate("");
  const auto dipole = res.get<Utils::Property::Dipole>();
  auto dipoleFromPreciseCalculation = 12.8503; // Debye
  // Dipole must be similar to true one (in Debye-> *2.541746 D/a.u), allow for 10% error.
  // Result of precise calculation in TestScripts
  ASSERT_THAT(dipole.norm() * 2.541746 / dipoleFromPreciseCalculation, DoubleNear(1, 0.1));
}

} // namespace Sparrow
} // namespace Scine
