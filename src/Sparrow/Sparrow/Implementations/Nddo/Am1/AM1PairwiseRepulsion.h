/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_AM1PAIRWISEREPULSION_H
#define SPARROW_AM1PAIRWISEREPULSION_H

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Eigen/Core>

namespace Scine {
namespace Utils {
enum class derivOrder;
}
namespace Sparrow {

namespace nddo {

// clang-format off
/*!
 * This class calculates the core-core repulsion between two atoms.
 */
class AM1PairwiseRepulsion{
 public:
  AM1PairwiseRepulsion(const AtomicParameters& A, const AtomicParameters&B);

  void calculate(const Eigen::Ref<Eigen::Vector3d>& R, Utils::derivOrder order);

  double getRepulsionEnergy() const { return repulsionEnergy_; }

  Eigen::RowVector3d getRepulsionGradient() const { return repulsionGradient_; }

  Utils::AutomaticDifferentiation::Second3D getRepulsionHessian() const { return repulsionHessian_; }
  template <Utils::derivativeType O> Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative() const;

  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> calculateRepulsion(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> baseTerm(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> parenthesisValue(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> standardParenthesis(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> gaussianRepulsionTerm(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> gaussianRepulsion(const AtomicParameters& P, double R) const;


 private:
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> radius(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> integral(double R) const;

  const AtomicParameters& pA_;
  const AtomicParameters& pB_;

  double repulsionEnergy_;
  Eigen::RowVector3d repulsionGradient_;
  Utils::AutomaticDifferentiation::Second3D repulsionHessian_;
};
// clang-format on

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::calculateRepulsion(double R) const {
  return baseTerm<O>(R) + gaussianRepulsionTerm<O>(R);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::baseTerm(double R) const {
  auto ssssIntegral = integral<O>(R);
  return pA_.coreCharge() * pB_.coreCharge() * ssssIntegral * parenthesisValue<O>(R);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::parenthesisValue(double R) const {
  return standardParenthesis<O>(R);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::standardParenthesis(double R) const {
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

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::gaussianRepulsionTerm(double R) const {
  auto RD = radius<O>(R);
  return (gaussianRepulsion<O>(pA_, R) + gaussianRepulsion<O>(pB_, R)) / RD *
         (pA_.coreCharge() * pB_.coreCharge() / Utils::Constants::ev_per_hartree);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::gaussianRepulsion(const AtomicParameters& P,
                                                                                        double R) const {
  if (P.hasGaussianRepulsionParameters()) {
    const auto& gaussianPar = P.getGaussianRepulsionParameters();
    auto gaussianRep = Utils::AutomaticDifferentiation::constant1D<O>(0);
    for (unsigned int i = 0; i < gaussianPar.size(); i += 1) {
      auto DistC = std::get<2>(gaussianPar[i]);
      auto RmC = radius<O>(R) - DistC;
      gaussianRep += std::get<0>(gaussianPar[i]) * exp(-std::get<1>(gaussianPar[i]) * RmC * RmC);
    };
    return gaussianRep;
  };

  return Utils::AutomaticDifferentiation::constant1D<O>(0);
}

template<Utils::derivOrder O>
inline Utils::AutomaticDifferentiation::Value1DType<O> AM1PairwiseRepulsion::radius(double R) const {
  return Utils::AutomaticDifferentiation::variableWithUnitDerivative<O>(R);
}

template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::zero>
AM1PairwiseRepulsion::integral<Utils::derivOrder::zero>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  return igr;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::one>
AM1PairwiseRepulsion::integral<Utils::derivOrder::one>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  Utils::AutomaticDifferentiation::First1D integralD(igr, -igr * igr * igr * R);
  return integralD;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::two>
AM1PairwiseRepulsion::integral<Utils::derivOrder::two>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  double igr3 = igr * igr * igr;
  Utils::AutomaticDifferentiation::Second1D integralD(igr, -igr3 * R, igr3 * igr * igr * (2 * R * R - pSum * pSum));
  return integralD;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::first>
AM1PairwiseRepulsion::getDerivative<Utils::derivativeType::first>() const {
  return getRepulsionGradient();
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_atomic>
AM1PairwiseRepulsion::getDerivative<Utils::derivativeType::second_atomic>() const {
  return getRepulsionHessian();
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_full>
AM1PairwiseRepulsion::getDerivative<Utils::derivativeType::second_full>() const {
  return getRepulsionHessian();
}

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_PAIRWISEREPULSION_H
