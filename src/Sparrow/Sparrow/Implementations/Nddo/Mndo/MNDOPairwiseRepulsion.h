/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MNDOPAIRWISEREPULSION_H
#define SPARROW_MNDOPAIRWISEREPULSION_H

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {
enum class DerivativeOrder;
}
namespace Sparrow {

namespace nddo {

// clang-format off
/**
 * @brief This class calculates the core-core repulsion between two atoms.
 */
class MNDOPairwiseRepulsion{
 public:
  MNDOPairwiseRepulsion(const AtomicParameters& A, const AtomicParameters&B);

  void calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::DerivativeOrder order);

  double getRepulsionEnergy() const { return repulsionEnergy_; }

  Eigen::RowVector3d getRepulsionGradient() const { return repulsionGradient_; }

  Utils::AutomaticDifferentiation::Second3D getRepulsionHessian() const { return repulsionHessian_; }
  template <Utils::Derivative O> Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative() const;

  template <Utils::DerivativeOrder O> Utils::AutomaticDifferentiation::Value1DType<O> calculateRepulsion(double R) const;
  template <Utils::DerivativeOrder O> Utils::AutomaticDifferentiation::Value1DType<O> baseTerm(double R) const;
  template <Utils::DerivativeOrder O> Utils::AutomaticDifferentiation::Value1DType<O> parenthesisValue(double R) const;
  template <Utils::DerivativeOrder O> Utils::AutomaticDifferentiation::Value1DType<O> standardParenthesis(double R) const;

 private:
  template <Utils::DerivativeOrder O> Utils::AutomaticDifferentiation::Value1DType<O> radius(double R) const;
  template <Utils::DerivativeOrder O> Utils::AutomaticDifferentiation::Value1DType<O> integral(double R) const;

  const AtomicParameters& pA_;
  const AtomicParameters& pB_;

  double repulsionEnergy_;
  Eigen::RowVector3d repulsionGradient_;
  Utils::AutomaticDifferentiation::Second3D repulsionHessian_;
};
// clang-format on

template<Utils::DerivativeOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> MNDOPairwiseRepulsion::calculateRepulsion(double R) const {
  return baseTerm<O>(R);
}

template<Utils::DerivativeOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> MNDOPairwiseRepulsion::baseTerm(double R) const {
  auto ssssIntegral = integral<O>(R);
  return pA_.coreCharge() * pB_.coreCharge() * ssssIntegral * parenthesisValue<O>(R);
}

template<Utils::DerivativeOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> MNDOPairwiseRepulsion::parenthesisValue(double R) const {
  return standardParenthesis<O>(R);
}

template<Utils::DerivativeOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> MNDOPairwiseRepulsion::standardParenthesis(double R) const {
  auto radiusDeriv = radius<O>(R);
  auto res = 1.0;
  if (pA_.element() == Utils::ElementType::H &&
      (pB_.element() == Utils::ElementType::N || pB_.element() == Utils::ElementType::O)) {
    return res + exp(-pA_.alpha() * radiusDeriv) +
           radiusDeriv * Utils::Constants::angstrom_per_bohr * exp(-pB_.alpha() * radiusDeriv);
  }
  else if (pB_.element() == Utils::ElementType::H &&
           (pA_.element() == Utils::ElementType::N || pA_.element() == Utils::ElementType::O)) {
    return res + radiusDeriv * Utils::Constants::angstrom_per_bohr * exp(-pA_.alpha() * radiusDeriv) +
           exp(-pB_.alpha() * radiusDeriv);
  }
  else {
    return res + exp(-pA_.alpha() * radiusDeriv) + exp(-pB_.alpha() * radiusDeriv);
  }
}

template<Utils::DerivativeOrder O>
inline Utils::AutomaticDifferentiation::Value1DType<O> MNDOPairwiseRepulsion::radius(double R) const {
  return Utils::AutomaticDifferentiation::variableWithUnitDerivative<O>(R);
}

template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::DerivativeOrder::Zero>
MNDOPairwiseRepulsion::integral<Utils::DerivativeOrder::Zero>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  return igr;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::DerivativeOrder::One>
MNDOPairwiseRepulsion::integral<Utils::DerivativeOrder::One>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  Utils::AutomaticDifferentiation::First1D integralD(igr, -igr * igr * igr * R);
  return integralD;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::DerivativeOrder::Two>
MNDOPairwiseRepulsion::integral<Utils::DerivativeOrder::Two>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  double igr3 = igr * igr * igr;
  Utils::AutomaticDifferentiation::Second1D integralD(igr, -igr3 * R, igr3 * igr * igr * (2 * R * R - pSum * pSum));
  return integralD;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::Derivative::First>
MNDOPairwiseRepulsion::getDerivative<Utils::Derivative::First>() const {
  return getRepulsionGradient();
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::Derivative::SecondAtomic>
MNDOPairwiseRepulsion::getDerivative<Utils::Derivative::SecondAtomic>() const {
  return getRepulsionHessian();
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::Derivative::SecondFull>
MNDOPairwiseRepulsion::getDerivative<Utils::Derivative::SecondFull>() const {
  return getRepulsionHessian();
}

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_PAIRWISEREPULSION_H
