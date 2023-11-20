/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_PAIRWISEREPULSION_H
#define SPARROW_DFTB_PAIRWISEREPULSION_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>

namespace Scine {

namespace {
enum class DerivativeOrder;
}

namespace Sparrow {

namespace dftb {
struct RepulsionParameters;

class PairwiseRepulsion {
 public:
  explicit PairwiseRepulsion(const RepulsionParameters& repulsionPars);

  void calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::DerivativeOrder order);

  double getRepulsionEnergy() const;
  template<Utils::Derivative O>
  Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative() const;

 private:
  template<Utils::DerivativeOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> calculateRepulsion(double r) const;

  const RepulsionParameters& repulsionPars_;
  double repulsionEnergy_ = 0;
  Eigen::RowVector3d repulsionGradient_;
  Utils::AutomaticDifferentiation::Second3D repulsionHessian_;
};

inline double PairwiseRepulsion::getRepulsionEnergy() const {
  return repulsionEnergy_;
}

template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::Derivative::First>
PairwiseRepulsion::getDerivative<Utils::Derivative::First>() const {
  return repulsionGradient_;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::Derivative::SecondAtomic>
PairwiseRepulsion::getDerivative<Utils::Derivative::SecondAtomic>() const {
  return repulsionHessian_;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::Derivative::SecondFull>
PairwiseRepulsion::getDerivative<Utils::Derivative::SecondFull>() const {
  return repulsionHessian_;
}

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_PAIRWISEREPULSION_H
