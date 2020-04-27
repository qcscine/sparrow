/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PairwiseRepulsion.h"
#include "RepulsionParameters.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>

namespace Scine {
namespace Sparrow {

using std::exp;
using namespace Utils::AutomaticDifferentiation;

namespace dftb {

PairwiseRepulsion::PairwiseRepulsion(const RepulsionParameters& repulsionPars) : repulsionPars_(repulsionPars) {
}

void PairwiseRepulsion::calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::derivOrder order) {
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

template<Utils::derivOrder O>
Value1DType<O> PairwiseRepulsion::calculateRepulsion(double r) const {
  auto R = variableWithUnitDerivative<O>(r);
  if (r > repulsionPars_.cutoff)
    return constant1D<O>(0);
  if (r < repulsionPars_.splineStart[0])
    return exp(-repulsionPars_.a1 * R + repulsionPars_.a2) + repulsionPars_.a3;

  auto i = static_cast<int>((r - repulsionPars_.splineStart[0]) /
                            (repulsionPars_.cutoff - repulsionPars_.splineStart[0]) * repulsionPars_.nSplineInts);

  // If not the right bin, change it:
  if (repulsionPars_.splineStart[i] > r)
    while (repulsionPars_.splineStart[--i] > r) {
    }
  else if (repulsionPars_.splineEnd[i] < r)
    while (repulsionPars_.splineEnd[++i] < r) {
    }

  auto dr = R - repulsionPars_.splineStart[i];
  auto repulsion = repulsionPars_.c0[i] +
                   dr * (repulsionPars_.c1[i] +
                         dr * (repulsionPars_.c2[i] +
                               dr * (repulsionPars_.c3[i] + (i == repulsionPars_.nSplineInts - 1
                                                                 ? dr * (repulsionPars_.c4 + dr * repulsionPars_.c5)
                                                                 : constant1D<O>(0.0)))));

  return repulsion;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
