/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MNDOPairwiseRepulsion.h"

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {

MNDOPairwiseRepulsion::MNDOPairwiseRepulsion(const AtomicParameters& A, const AtomicParameters& B) : pA_(A), pB_(B) {
}
// No need of atomic writing as this class is local to a single atom pair
// -> no race condition
void MNDOPairwiseRepulsion::calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::derivOrder order) {
  double Rabs = R.norm();
  if (order == Utils::derivOrder::zero) {
    repulsionEnergy_ = calculateRepulsion<Utils::derivOrder::zero>(Rabs);
  }
  else if (order == Utils::derivOrder::one) {
    First1D rep = calculateRepulsion<Utils::derivOrder::one>(Rabs);
    repulsionEnergy_ = rep.value();
    repulsionGradient_ = Utils::Gradient(get3Dfrom1D<Utils::derivOrder::one>(rep, R).derivatives());
  }
  else if (order == Utils::derivOrder::two) {
    Second1D rep = calculateRepulsion<Utils::derivOrder::two>(Rabs);
    repulsionEnergy_ = rep.value();
    repulsionHessian_ = get3Dfrom1D<Utils::derivOrder::two>(rep, R);
  }
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
