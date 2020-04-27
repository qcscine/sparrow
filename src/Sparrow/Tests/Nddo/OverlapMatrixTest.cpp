/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/AtomPairOverlap.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/OverlapMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/GtoExpansion.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/DataStructures/SlaterToGaussian.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>
#include <Eigen/Core>
#include <memory>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace Utils;
using namespace nddo;

class AOverlapMatrix : public Test {
 public:
  AOverlapMatrix() : S(elementTypes_, positions_, aoIndexes_, elementParameters) {
  }

  Utils::AtomsOrbitalsIndexes aoIndexes_;
  ElementParameters elementParameters;
  OverlapMatrix S;
  std::unique_ptr<AtomicParameters> apH, apC, apV;
  AtomPairOverlap<Utils::derivOrder::zero> pairOverlap;
  Utils::GtoExpansion gs1, gs2, gp1, gp2, gd1, gd2;
  Utils::AtomicGtos aS1, aS2, aP1, aP2, aD1, aD2;
  Eigen::Vector3d arbitraryVector;
  Eigen::RowVector3d pos1, pos2, pos3, pos4;
  Utils::ElementTypeCollection elementTypes_;
  Utils::PositionCollection positions_;

  void SetUp() override {
    arbitraryVector = Eigen::Vector3d(1.1, 2.0, 0.23);
    pos1 = Utils::Position(0.0, 0.0, 0.0);
    pos2 = Utils::Position(0.1, -0.6, 0.5);
    pos3 = Utils::Position(-1.0, 0.2, -0.3);
    pos4 = Utils::Position(1e6, 1e6, 1e6);
    apH = std::make_unique<AtomicParameters>(Utils::ElementType::H);
    apC = std::make_unique<AtomicParameters>(Utils::ElementType::C);
    apV = std::make_unique<AtomicParameters>(Utils::ElementType::V);

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

    elementParameters.set(Utils::ElementType::H, std::move(apH));
    elementParameters.set(Utils::ElementType::C, std::move(apC));
    elementParameters.set(Utils::ElementType::V, std::move(apV));
  }
};

TEST_F(AOverlapMatrix, ReturnsIdentityForSingleAtom) {
  elementTypes_.push_back(Utils::ElementType::H);
  positions_.resize(1, 3);
  positions_.row(0) = pos1;
  aoIndexes_.setSize(1);
  aoIndexes_.addAtom(1);
  // S.newStructure... -> set matrix size, compute diagonal, etc.

  S.reset();
  S.calculateOverlap(derivOrder::zero);
  const auto& m = S.getOverlap();
  ASSERT_TRUE(m.get0().isApprox(Eigen::Matrix<double, 1, 1>::Identity(), 1e-7));
  elementTypes_.clear();
  elementTypes_.push_back(Utils::ElementType::C);
  positions_.row(0) = pos1;
  aoIndexes_.clear();
  aoIndexes_.setSize(1);
  aoIndexes_.addAtom(4);
  S.reset();
  S.calculateOverlap(derivOrder::zero);
  const auto& m2 = S.getOverlap();
  ASSERT_TRUE(m2.get0().isApprox(Eigen::Matrix<double, 4, 4>::Identity(), 1e-7));
}

TEST_F(AOverlapMatrix, IsNearToIdentityIfAtomsAreFarAway) {
  positions_.resize(2, 3);
  positions_.row(0) = pos1;
  positions_.row(1) = pos4;
  elementTypes_.push_back(Utils::ElementType::H);
  elementTypes_.push_back(Utils::ElementType::C);
  aoIndexes_.setSize(2);
  aoIndexes_.addAtom(1);
  aoIndexes_.addAtom(4);

  Utils::MatrixWithDerivatives m;
  S.reset();
  S.calculateOverlap(derivOrder::zero);
  m = S.getOverlap();
  ASSERT_TRUE(m.get0().isApprox(Eigen::Matrix<double, 5, 5>::Identity(), 1e-7));
}
} // namespace Sparrow
} // namespace Scine
