/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GTOOverlapMatrixBlock.h>
#include <Utils/DataStructures/GtoExpansion.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/DataStructures/SlaterToGaussian.h>
#include <gmock/gmock.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class AGTOOverlapMatrixBlock : public Test {
 public:
  GTOOverlapMatrixBlock<Utils::derivOrder::one> overlapBlock;
  Utils::GtoExpansion gs1, gs2, gp1, gp2, gd1, gd2;
  Eigen::Vector3d arbitraryVector;

  void SetUp() override {
    arbitraryVector = Eigen::Vector3d(1.1, 2.0, 0.23);
    gs1 = Utils::SlaterToGaussian::getGTOExpansion(6, 1, 0);
    gs2 = Utils::SlaterToGaussian::getGTOExpansion(6, 4, 0, 1.42);
    gp1 = Utils::SlaterToGaussian::getGTOExpansion(6, 2, 1, 2.2);
    gp2 = Utils::SlaterToGaussian::getGTOExpansion(6, 3, 1);
    gd1 = Utils::SlaterToGaussian::getGTOExpansion(6, 5, 2, 0.92);
    gd2 = Utils::SlaterToGaussian::getGTOExpansion(6, 4, 2);
  }
};

TEST_F(AGTOOverlapMatrixBlock, ReturnsBlockOfCorrectSize) {
  ASSERT_THAT(overlapBlock.getMatrixBlock(gs1, gs2, arbitraryVector).size(), Eq(1));
  ASSERT_THAT(overlapBlock.getMatrixBlock(gd1, gd2, arbitraryVector).size(), Eq(25));
  ASSERT_THAT(overlapBlock.getMatrixBlock(gp1, gd2, arbitraryVector).rows(), Eq(3));
  ASSERT_THAT(overlapBlock.getMatrixBlock(gp1, gd2, arbitraryVector).cols(), Eq(5));
}

TEST_F(AGTOOverlapMatrixBlock, ReturnsIdentityForOrbitalsOnSameAtom) {
  Utils::MatrixWithDerivatives m;
  m.setOrder(Utils::derivOrder::one);
  m.setDimension(gs1.nAOs(), gs1.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 1) = overlapBlock.getMatrixBlock(gs1, gs1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockSS = m.getMatrixXd();
  m.setDimension(gp1.nAOs(), gp1.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 3, 3) = overlapBlock.getMatrixBlock(gp1, gp1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockPP = m.getMatrixXd();
  m.setDimension(gd1.nAOs(), gd1.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 5, 5) = overlapBlock.getMatrixBlock(gd1, gd1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockDD = m.getMatrixXd();
  m.setDimension(gs1.nAOs(), gp1.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 3) = overlapBlock.getMatrixBlock(gs1, gp1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockSP = m.getMatrixXd();
  m.setDimension(gs1.nAOs(), gd1.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 5) = overlapBlock.getMatrixBlock(gs1, gd1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockSD = m.getMatrixXd();
  m.setDimension(gp1.nAOs(), gd1.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 3, 5) = overlapBlock.getMatrixBlock(gp1, gd1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockPD = m.getMatrixXd();

  ASSERT_TRUE(blockSS.isApprox(Eigen::Matrix<double, 1, 1>::Identity(), 1e-7));
  ASSERT_TRUE(blockPP.isApprox(Eigen::Matrix<double, 3, 3>::Identity(), 1e-7));
  ASSERT_TRUE(blockDD.isApprox(Eigen::Matrix<double, 5, 5>::Identity(), 1e-7));
  ASSERT_TRUE(blockSP.isZero(1e-7));
  ASSERT_TRUE(blockSD.isZero(1e-7));
  ASSERT_TRUE(blockPD.isZero(1e-7));
}

TEST_F(AGTOOverlapMatrixBlock, ReturnsCorrectDerivativesForSOrbitals) {
  Utils::MatrixWithDerivatives m;
  m.setOrder(Utils::derivOrder::one);
  m.setDimension(gs1.nAOs(), gs2.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 1) = overlapBlock.getMatrixBlock(gs1, gs2, arbitraryVector);
  Eigen::MatrixXd S1 = m.getMatrixXd();
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 1) =
      overlapBlock.getMatrixBlock(gs1, gs2, arbitraryVector + Eigen::Vector3d(0.01, 0, 0));
  Eigen::MatrixXd S2 = m.getMatrixXd();
  ASSERT_THAT(((S2 - S1)(0, 0)), DoubleNear(0.01 * m.v1(0, 0).derivatives().x(), 1e-5));
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 1) =
      overlapBlock.getMatrixBlock(gs1, gs2, arbitraryVector + Eigen::Vector3d(0, 0.01, 0));
  S2 = m.getMatrixXd();
  ASSERT_THAT(((S2 - S1)(0, 0)), DoubleNear(0.01 * m.v1(0, 0).derivatives().y(), 1e-5));
  m.get<Utils::derivOrder::one>().block(0, 0, 1, 1) =
      overlapBlock.getMatrixBlock(gs1, gs2, arbitraryVector + Eigen::Vector3d(0, 0, 0.01));
  S2 = m.getMatrixXd();
  ASSERT_THAT(((S2 - S1)(0, 0)), DoubleNear(0.01 * m.v1(0, 0).derivatives().z(), 1e-5));
}

TEST_F(AGTOOverlapMatrixBlock, ReturnsCorrectDerivativesForPOrbitals) {
  Utils::MatrixWithDerivatives m;
  m.setOrder(Utils::derivOrder::one);
  m.setDimension(gp1.nAOs(), gp2.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 3, 3) = overlapBlock.getMatrixBlock(gp1, gp2, arbitraryVector);
  Eigen::MatrixXd S1 = m.getMatrixXd();

  m.get<Utils::derivOrder::one>().block(0, 0, 3, 3) =
      overlapBlock.getMatrixBlock(gp1, gp2, arbitraryVector + Eigen::Vector3d(0.01, 0, 0));
  Eigen::MatrixXd S2 = m.getMatrixXd();
  for (int i = 0; i < gp1.nAOs(); i++)
    for (int j = 0; j < gp2.nAOs(); j++)
      ASSERT_THAT(((S2 - S1)(i, j)), DoubleNear(0.01 * m.v1(i, j).derivatives().x(), 1e-5));

  m.get<Utils::derivOrder::one>().block(0, 0, 3, 3) =
      overlapBlock.getMatrixBlock(gp1, gp2, arbitraryVector + Eigen::Vector3d(0, 0.01, 0));
  S2 = m.getMatrixXd();
  for (int i = 0; i < gp1.nAOs(); i++)
    for (int j = 0; j < gp2.nAOs(); j++)
      ASSERT_THAT(((S2 - S1)(i, j)), DoubleNear(0.01 * m.v1(i, j).derivatives().y(), 1e-5));

  m.get<Utils::derivOrder::one>().block(0, 0, 3, 3) =
      overlapBlock.getMatrixBlock(gp1, gp2, arbitraryVector + Eigen::Vector3d(0, 0, 0.01));
  S2 = m.getMatrixXd();
  for (int i = 0; i < gp1.nAOs(); i++)
    for (int j = 0; j < gp2.nAOs(); j++)
      ASSERT_THAT(((S2 - S1)(i, j)), DoubleNear(0.01 * m.v1(i, j).derivatives().z(), 1e-5));
}

TEST_F(AGTOOverlapMatrixBlock, ReturnsCorrectDerivativesForDOrbitals) {
  Utils::MatrixWithDerivatives m;
  m.setOrder(Utils::derivOrder::one);
  m.setDimension(gd1.nAOs(), gd2.nAOs());
  m.get<Utils::derivOrder::one>().block(0, 0, 5, 5) = overlapBlock.getMatrixBlock(gd1, gd2, arbitraryVector);
  Eigen::MatrixXd S1 = m.getMatrixXd();

  m.get<Utils::derivOrder::one>().block(0, 0, 5, 5) =
      overlapBlock.getMatrixBlock(gd1, gd2, arbitraryVector + Eigen::Vector3d(0.01, 0, 0));
  Eigen::MatrixXd S2 = m.getMatrixXd();
  for (int i = 0; i < gd1.nAOs(); i++)
    for (int j = 0; j < gd2.nAOs(); j++)
      ASSERT_THAT(((S2 - S1)(i, j)), DoubleNear(0.01 * m.v1(i, j).derivatives().x(), 1e-5));

  m.get<Utils::derivOrder::one>().block(0, 0, 5, 5) =
      overlapBlock.getMatrixBlock(gd1, gd2, arbitraryVector + Eigen::Vector3d(0, 0.01, 0));
  S2 = m.getMatrixXd();
  for (int i = 0; i < gd1.nAOs(); i++)
    for (int j = 0; j < gd2.nAOs(); j++)
      ASSERT_THAT(((S2 - S1)(i, j)), DoubleNear(0.01 * m.v1(i, j).derivatives().y(), 1e-5));

  m.get<Utils::derivOrder::one>().block(0, 0, 5, 5) =
      overlapBlock.getMatrixBlock(gd1, gd2, arbitraryVector + Eigen::Vector3d(0, 0, 0.01));
  S2 = m.getMatrixXd();
  for (int i = 0; i < gd1.nAOs(); i++)
    for (int j = 0; j < gd2.nAOs(); j++)
      ASSERT_THAT(((S2 - S1)(i, j)), DoubleNear(0.01 * m.v1(i, j).derivatives().z(), 1e-5));
}
} // namespace Sparrow
} // namespace Scine
