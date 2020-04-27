/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomPairOverlap.h"
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

template<Utils::derivOrder O>
Eigen::Matrix<typename AtomPairOverlap<O>::Value3D, Eigen::Dynamic, Eigen::Dynamic>
AtomPairOverlap<O>::getMatrixBlock(const Utils::AtomicGtos& pA, const Utils::AtomicGtos& pB, Eigen::Vector3d Rab) const {
  auto localBlock = getInitialBlock(pA, pB);
  GTOOverlapMatrixBlock<O> block;

  if (pA.hasS()) {
    if (pB.hasS()) {
      localBlock.block(0, 0, 1, 1) = block.getMatrixBlock(pA.s(), pB.s(), Rab);
    }
    if (pB.hasP()) {
      localBlock.block(0, 1, 1, 3) = block.getMatrixBlock(pA.s(), pB.p(), Rab);
    }
    if (pB.hasD()) {
      localBlock.block(0, 4, 1, 5) = block.getMatrixBlock(pA.s(), pB.d(), Rab);
    }
  }
  if (pA.hasP()) {
    if (pB.hasS()) {
      localBlock.block(1, 0, 3, 1) = block.getMatrixBlock(pA.p(), pB.s(), Rab);
    }
    if (pB.hasP()) {
      localBlock.block(1, 1, 3, 3) = block.getMatrixBlock(pA.p(), pB.p(), Rab);
    }
    if (pB.hasD()) {
      localBlock.block(1, 4, 3, 5) = block.getMatrixBlock(pA.p(), pB.d(), Rab);
    }
  }
  if (pA.hasD()) {
    if (pB.hasS()) {
      localBlock.block(4, 0, 5, 1) = block.getMatrixBlock(pA.d(), pB.s(), Rab);
    }
    if (pB.hasP()) {
      localBlock.block(4, 1, 5, 3) = block.getMatrixBlock(pA.d(), pB.p(), Rab);
    }
    if (pB.hasD()) {
      localBlock.block(4, 4, 5, 5) = block.getMatrixBlock(pA.d(), pB.d(), Rab);
    }
  }

  return localBlock;
}

template<Utils::derivOrder O>
Eigen::Matrix<typename AtomPairOverlap<O>::Value3D, Eigen::Dynamic, Eigen::Dynamic>
AtomPairOverlap<O>::getInitialBlock(const Utils::AtomicGtos& pA, const Utils::AtomicGtos& pB) const {
  int dimA = (pA.hasS() ? 1 : 0) + (pA.hasP() ? 3 : 0) + (pA.hasD() ? 5 : 0);
  int dimB = (pB.hasS() ? 1 : 0) + (pB.hasP() ? 3 : 0) + (pB.hasD() ? 5 : 0);
  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic> initialBlock(dimA, dimB);
  return initialBlock;
}
template class AtomPairOverlap<Utils::derivOrder::one>;
template class AtomPairOverlap<Utils::derivOrder::two>;
template class AtomPairOverlap<Utils::derivOrder::zero>;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
