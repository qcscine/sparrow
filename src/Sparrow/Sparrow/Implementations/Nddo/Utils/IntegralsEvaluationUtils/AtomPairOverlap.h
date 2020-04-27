/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ATOMPAIROVERLAP_H
#define SPARROW_ATOMPAIROVERLAP_H

#include "GTOOverlapMatrixBlock.h"
#include <Utils/Math/DerivOrderEnum.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
class AtomicGtos;
}

namespace Sparrow {

namespace nddo {

/**
 * @brief This class computes a block of the overlap matrix for two atoms.
 * The actual calculation is done by the GTOOverlapMatrixBlock class, here the blocks that need calculation are
 * identified and scheduled for calculation.
 */

template<Utils::derivOrder O>
class AtomPairOverlap {
 public:
  using Value3D = Utils::AutomaticDifferentiation::Value3DType<O>;

  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic>
  getMatrixBlock(const Utils::AtomicGtos& pA, const Utils::AtomicGtos& pB, Eigen::Vector3d Rab) const;

 private:
  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic> getInitialBlock(const Utils::AtomicGtos& pA,
                                                                         const Utils::AtomicGtos& pB) const;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ATOMPAIROVERLAP_H
