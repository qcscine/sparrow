/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/AtomPairOverlap.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/GtoExpansion.h>
#include <Utils/DataStructures/SlaterToGaussian.h>
#include <gmock/gmock.h>
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class AAtomPairOverlapMatrix : public Test {
 public:
  AtomPairOverlap<Utils::derivOrder::zero> pairOverlap;
  Utils::GtoExpansion gs1, gs2, gp1, gp2, gd1, gd2;
  Utils::AtomicGtos aS1, aS2, aP1, aP2, aD1, aD2;
  Eigen::Vector3d arbitraryVector;

  void SetUp() override {
    arbitraryVector = Eigen::Vector3d(1.1, 2.0, 0.23);
    gs1 = Utils::SlaterToGaussian::getGTOExpansion(6, 1, 0);
    gs2 = Utils::SlaterToGaussian::getGTOExpansion(6, 4, 0, 1.42);
    gp1 = Utils::SlaterToGaussian::getGTOExpansion(6, 2, 1, 2.2);
    gp2 = Utils::SlaterToGaussian::getGTOExpansion(6, 3, 1);
    gd1 = Utils::SlaterToGaussian::getGTOExpansion(6, 5, 2, 0.92);
    gd2 = Utils::SlaterToGaussian::getGTOExpansion(6, 4, 2);

    aS1.setS(gs1);
    aS2.setS(gs2);
    aP1 = aS1;
    aP1.setP(gp1);
    aP2 = aS2;
    aP2.setP(gp2);
    aD1 = aP1;
    aD1.setD(gd1);
    aD2 = aP2;
    aD2.setD(gd2);
  }
};

TEST_F(AAtomPairOverlapMatrix, ReturnsAtomPairBlockOfCorrectSize) {
  ASSERT_THAT(pairOverlap.getMatrixBlock(aS1, aS2, arbitraryVector).size(), Eq(1));
  ASSERT_THAT(pairOverlap.getMatrixBlock(aD1, aD2, arbitraryVector).size(), Eq(81));
  ASSERT_THAT(pairOverlap.getMatrixBlock(aD1, aP2, arbitraryVector).rows(), Eq(9));
  ASSERT_THAT(pairOverlap.getMatrixBlock(aD1, aP2, arbitraryVector).cols(), Eq(4));
}

TEST_F(AAtomPairOverlapMatrix, ReturnsIdentityForOverlapOfAtomWithItself) {
  Eigen::MatrixXd blockS = pairOverlap.getMatrixBlock(aS1, aS1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockP = pairOverlap.getMatrixBlock(aP1, aP1, Eigen::Vector3d(0, 0, 0));
  Eigen::MatrixXd blockD = pairOverlap.getMatrixBlock(aD1, aD1, Eigen::Vector3d(0, 0, 0));

  ASSERT_TRUE(blockS.isApprox(Eigen::Matrix<double, 1, 1>::Identity(), 1e-7));
  ASSERT_TRUE(blockP.isApprox(Eigen::Matrix<double, 4, 4>::Identity(), 1e-7));
  ASSERT_TRUE(blockD.isApprox(Eigen::Matrix<double, 9, 9>::Identity(), 1e-7));
}
} // namespace Sparrow
} // namespace Scine
