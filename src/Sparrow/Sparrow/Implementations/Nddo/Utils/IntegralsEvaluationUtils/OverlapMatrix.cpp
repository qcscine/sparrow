/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OverlapMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

OverlapMatrix::OverlapMatrix(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                             const Utils::AtomsOrbitalsIndexes& aoIndexes, const ElementParameters& elementParameters)
  : elementTypes_(elements), positions_(positions), aoIndexes_(aoIndexes), elementParameters_(elementParameters) {
}

void OverlapMatrix::reset() {
  nAOs_ = 0;
  nAtoms_ = static_cast<int>(elementTypes_.size());
  for (auto e : elementTypes_)
    nAOs_ += elementParameters_.get(e).nAOs();

  S_.setBaseMatrix(Eigen::MatrixXd::Identity(nAOs_, nAOs_));
}

void OverlapMatrix::calculateOverlap(Utils::derivOrder highestRequiredOrder) {
  auto order = Utils::AutomaticDifferentiation::getDerivativeOrder(highestRequiredOrder);
  S_.setOrder(highestRequiredOrder);
  if (nAOs_ == 0)
    return;
#pragma omp parallel for schedule(dynamic)
  for (int i = 1; i < nAtoms_; ++i) {
    auto rowIndex = aoIndexes_.getFirstOrbitalIndex(i);
    const auto& pA = elementParameters_.get(elementTypes_[i]);
    for (int j = 0; j < i; ++j) {
      auto colIndex = aoIndexes_.getFirstOrbitalIndex(j);
      const auto& pB = elementParameters_.get(elementTypes_[j]);
      const auto& GTOsA = pA.GTOs();
      const auto& GTOsB = pB.GTOs();

      const auto& Ri = positions_.row(i);
      const auto& Rj = positions_.row(j);
      Eigen::RowVector3d Rij = Rj - Ri;
      if (order == 0) {
        auto resultBlock = pairOverlapZeroOrder_.getMatrixBlock(GTOsA, GTOsB, Rij);
        S_.get<Utils::derivOrder::zero>().block(rowIndex, colIndex, resultBlock.rows(), resultBlock.cols()) = resultBlock;
      }
      else if (order == 1) {
        auto resultBlock = pairOverlapFirstOrder_.getMatrixBlock(GTOsA, GTOsB, Rij);
        S_.get<Utils::derivOrder::one>().block(rowIndex, colIndex, resultBlock.rows(), resultBlock.cols()) = resultBlock;
      }
      else if (order == 2) {
        auto resultBlock = pairOverlapSecondOrder_.getMatrixBlock(GTOsA, GTOsB, Rij);
        S_.get<Utils::derivOrder::two>().block(rowIndex, colIndex, resultBlock.rows(), resultBlock.cols()) = resultBlock;
      }
    }
  }
}

const Utils::MatrixWithDerivatives& OverlapMatrix::getOverlap() const {
  return S_;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
