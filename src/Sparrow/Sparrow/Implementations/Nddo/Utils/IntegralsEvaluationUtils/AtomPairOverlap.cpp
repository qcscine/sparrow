/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AtomPairOverlap.h"
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

template<Utils::DerivativeOrder O>
Eigen::Matrix<typename AtomPairOverlap<O>::Value3D, Eigen::Dynamic, Eigen::Dynamic>
AtomPairOverlap<O>::getMatrixBlock(const Utils::AtomicGtos& pA, const Utils::AtomicGtos& pB, Eigen::Vector3d Rab) const {
  auto localBlock = getInitialBlock(pA, pB);
  GTOOverlapMatrixBlock<O> block;

  if (pA.s) {
    if (pB.s) {
      localBlock.block(0, 0, 1, 1) = block.getMatrixBlock(pA.s.value(), pB.s.value(), Rab);
    }
    if (pB.p) {
      localBlock.block(0, 1, 1, 3) = block.getMatrixBlock(pA.s.value(), pB.p.value(), Rab);
    }
    if (pB.d) {
      localBlock.block(0, 4, 1, 5) = block.getMatrixBlock(pA.s.value(), pB.d.value(), Rab);
    }
  }
  if (pA.p) {
    if (pB.s) {
      localBlock.block(1, 0, 3, 1) = block.getMatrixBlock(pA.p.value(), pB.s.value(), Rab);
    }
    if (pB.p) {
      localBlock.block(1, 1, 3, 3) = block.getMatrixBlock(pA.p.value(), pB.p.value(), Rab);
    }
    if (pB.d) {
      localBlock.block(1, 4, 3, 5) = block.getMatrixBlock(pA.p.value(), pB.d.value(), Rab);
    }
  }
  if (pA.d) {
    if (pB.s) {
      localBlock.block(4, 0, 5, 1) = block.getMatrixBlock(pA.d.value(), pB.s.value(), Rab);
    }
    if (pB.p) {
      localBlock.block(4, 1, 5, 3) = block.getMatrixBlock(pA.d.value(), pB.p.value(), Rab);
    }
    if (pB.d) {
      localBlock.block(4, 4, 5, 5) = block.getMatrixBlock(pA.d.value(), pB.d.value(), Rab);
    }
  }

  return localBlock;
}

template<Utils::DerivativeOrder O>
Eigen::Matrix<typename AtomPairOverlap<O>::Value3D, Eigen::Dynamic, Eigen::Dynamic>
AtomPairOverlap<O>::getInitialBlock(const Utils::AtomicGtos& pA, const Utils::AtomicGtos& pB) const {
  int dimA = (pA.s ? 1 : 0) + (pA.p ? 3 : 0) + (pA.d ? 5 : 0);
  int dimB = (pB.s ? 1 : 0) + (pB.p ? 3 : 0) + (pB.d ? 5 : 0);
  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic> initialBlock(dimA, dimB);
  return initialBlock;
}
template class AtomPairOverlap<Utils::DerivativeOrder::One>;
template class AtomPairOverlap<Utils::DerivativeOrder::Two>;
template class AtomPairOverlap<Utils::DerivativeOrder::Zero>;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
