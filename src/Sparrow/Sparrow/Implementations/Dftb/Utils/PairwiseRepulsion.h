/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_PAIRWISEREPULSION_H
#define SPARROW_DFTB_PAIRWISEREPULSION_H

#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/AutomaticDifferentiation/Second3D.h>

namespace Scine {

namespace {
enum class derivOrder;
}

namespace Sparrow {

namespace dftb {
struct RepulsionParameters;

class PairwiseRepulsion {
 public:
  explicit PairwiseRepulsion(const RepulsionParameters& repulsionPars);

  void calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::derivOrder order);

  double getRepulsionEnergy() const;
  template<Utils::derivativeType O>
  Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative() const;

 private:
  template<Utils::derivOrder O>
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
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::first>
PairwiseRepulsion::getDerivative<Utils::derivativeType::first>() const {
  return repulsionGradient_;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_atomic>
PairwiseRepulsion::getDerivative<Utils::derivativeType::second_atomic>() const {
  return repulsionHessian_;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_full>
PairwiseRepulsion::getDerivative<Utils::derivativeType::second_full>() const {
  return repulsionHessian_;
}

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_PAIRWISEREPULSION_H
