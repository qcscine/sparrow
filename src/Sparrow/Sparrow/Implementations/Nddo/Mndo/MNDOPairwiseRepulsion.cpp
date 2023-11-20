/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
void MNDOPairwiseRepulsion::calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::DerivativeOrder order) {
  double Rabs = R.norm();
  if (order == Utils::DerivativeOrder::Zero) {
    repulsionEnergy_ = calculateRepulsion<Utils::DerivativeOrder::Zero>(Rabs);
  }
  else if (order == Utils::DerivativeOrder::One) {
    First1D rep = calculateRepulsion<Utils::DerivativeOrder::One>(Rabs);
    repulsionEnergy_ = rep.value();
    repulsionGradient_ = Utils::Gradient(get3Dfrom1D<Utils::DerivativeOrder::One>(rep, R).derivatives());
  }
  else if (order == Utils::DerivativeOrder::Two) {
    Second1D rep = calculateRepulsion<Utils::DerivativeOrder::Two>(Rabs);
    repulsionEnergy_ = rep.value();
    repulsionHessian_ = get3Dfrom1D<Utils::DerivativeOrder::Two>(rep, R);
  }
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
