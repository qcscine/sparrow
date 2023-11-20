/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

void PairwiseRepulsion::calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::DerivativeOrder order) {
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

template<Utils::DerivativeOrder O>
Value1DType<O> PairwiseRepulsion::calculateRepulsion(double r) const {
  // TODO duplicated in SKPair

  auto R = variableWithUnitDerivative<O>(r);
  if (r > repulsionPars_.cutoff)
    return constant1D<O>(0);
  if (r < repulsionPars_.splines[0].start)
    return exp(-repulsionPars_.a1 * R + repulsionPars_.a2) + repulsionPars_.a3;

  auto i = static_cast<int>((r - repulsionPars_.splines[0].start) /
                            (repulsionPars_.cutoff - repulsionPars_.splines[0].start) * repulsionPars_.nSplineInts);

  // If not the right bin, change it:
  if (repulsionPars_.splines[i].start > r)
    while (repulsionPars_.splines[--i].start > r) {
    }
  else if (repulsionPars_.splines[i].end < r)
    while (repulsionPars_.splines[++i].end < r) {
    }

  const auto& spline = repulsionPars_.splines[i];

  auto dr = R - spline.start;
  auto repulsion =
      spline.c0 +
      dr * (spline.c1 + dr * (spline.c2 + dr * (spline.c3 + (i == repulsionPars_.nSplineInts - 1
                                                                 ? dr * (repulsionPars_.c4 + dr * repulsionPars_.c5)
                                                                 : constant1D<O>(0.0)))));

  return repulsion;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
