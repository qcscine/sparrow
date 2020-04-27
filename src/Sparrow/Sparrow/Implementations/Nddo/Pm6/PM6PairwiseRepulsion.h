/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_PAIRWISEREPULSION_H
#define SPARROW_PAIRWISEREPULSION_H

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Eigen/Core>

namespace Scine {

namespace Sparrow {

namespace nddo {

// clang-format off
/*!
 * This class calculates the core-core repulsion between two atoms.
 */
class PM6PairwiseRepulsion{
 public:
  PM6PairwiseRepulsion(const AtomicParameters& A, const AtomicParameters&B, const PM6DiatomicParameters& AB);

  void calculate(const Eigen::Vector3d& R, Utils::derivOrder order);

  double getRepulsionEnergy() const { return repulsionEnergy_; }

  Eigen::RowVector3d getRepulsionGradient() const { return repulsionGradient_; }

  Utils::AutomaticDifferentiation::Second3D getRepulsionHessian() const { return repulsionHessian_; }
  template <Utils::derivativeType O> Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative() const;

  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> calculateRepulsion(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> baseTerm(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> additionalTerm(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> gaussianRepulsionTerm(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> gaussianRepulsion(const AtomicParameters& P, double R) const;

  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> parenthesisValue(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> standardParenthesis(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> NHOHParenthesis(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> CCParenthesis(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> SiOParenthesis(double R) const;


 private:
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> radius(double R) const;
  template <Utils::derivOrder O> Utils::AutomaticDifferentiation::Value1DType<O> integral(double R) const;

  const double cOfAdditiveTerm; // NB: The value of c = 1e-8 described in the paper is transformed from Angstrom^12 to bohr^12.
  const double exponentCoefficient; // NB: The value of 0.0003 described in the paper is transformed from Angstrom^-5 to bohr^-5.
  const double furtherExponentCC; // NB: The value of 5.98 described in the paper is transformed from Angstrom^-1 to bohr^-1.
  const double distanceSiO; // NB: The value of 2.9 described in the paper is transformed from Angstrom to bohr.
  const double factorCC;
  const double factorSiO;

  const AtomicParameters& pA_;
  const AtomicParameters& pB_;
  const PM6DiatomicParameters& pAB_;

  double repulsionEnergy_;
  Eigen::RowVector3d repulsionGradient_;
  Utils::AutomaticDifferentiation::Second3D repulsionHessian_;
};
// clang-format on

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::calculateRepulsion(double R) const {
  return baseTerm<O>(R) + additionalTerm<O>(R) + gaussianRepulsionTerm<O>(R);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::baseTerm(double R) const {
  auto ssssIntegral = integral<O>(R);
  return pA_.coreCharge() * pB_.coreCharge() * ssssIntegral * parenthesisValue<O>(R);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::parenthesisValue(double R) const {
  if ((pA_.element() == Utils::ElementType::H && (pB_.element() == Utils::ElementType::O || pB_.element() == Utils::ElementType::N ||
                                                  pB_.element() == Utils::ElementType::C)) ||
      (pB_.element() == Utils::ElementType::H && (pA_.element() == Utils::ElementType::O || pA_.element() == Utils::ElementType::N ||
                                                  pA_.element() == Utils::ElementType::C)))
    return NHOHParenthesis<O>(R);
  if (pA_.element() == Utils::ElementType::C && pB_.element() == Utils::ElementType::C)
    return CCParenthesis<O>(R);
  if ((pA_.element() == Utils::ElementType::Si && pB_.element() == Utils::ElementType::O) ||
      (pA_.element() == Utils::ElementType::O && pB_.element() == Utils::ElementType::Si))
    return SiOParenthesis<O>(R);

  return standardParenthesis<O>(R);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::additionalTerm(double R) const {
  auto RD = radius<O>(R);
  double za13 = std::pow(Utils::ElementInfo::Z(pA_.element()), 1.0 / 3.0);
  double zb13 = std::pow(Utils::ElementInfo::Z(pB_.element()), 1.0 / 3.0);
  auto v = Utils::AutomaticDifferentiation::constant1D<O>(za13 + zb13);
  v /= RD;
  auto v2 = v * v;
  auto v6 = v2 * v2 * v2;
  return cOfAdditiveTerm * v6 * v6 / Utils::Constants::ev_per_hartree;
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::gaussianRepulsionTerm(double R) const {
  auto RD = radius<O>(R);
  return (gaussianRepulsion<O>(pA_, R) + gaussianRepulsion<O>(pB_, R)) / RD *
         (pA_.coreCharge() * pB_.coreCharge() / Utils::Constants::ev_per_hartree);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::gaussianRepulsion(const AtomicParameters& P,
                                                                                        double R) const {
  if (P.hasGaussianRepulsionParameters()) {
    const auto& gaussianPar = P.getGaussianRepulsionParameters();
    for (unsigned int i = 0; i < gaussianPar.size(); i += 1) {
      auto DistC = std::get<2>(gaussianPar[i]);
      auto RmC = radius<O>(R) - DistC;
      return std::get<0>(gaussianPar[i]) * exp(-std::get<1>(gaussianPar[i]) * RmC * RmC);
    };
  };

  return Utils::AutomaticDifferentiation::constant1D<O>(0);
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::standardParenthesis(double R) const {
  auto RD = radius<O>(R);
  auto R2 = RD * RD;
  auto R6 = R2 * R2 * R2;
  auto res = 1.0;
  return res + (2 * pAB_.x()) * exp(-pAB_.alpha() * (RD + exponentCoefficient * R6));
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::NHOHParenthesis(double R) const {
  auto RD = radius<O>(R);
  auto res = 1.0;
  return res + 2 * pAB_.x() * exp(-pAB_.alpha() * (RD * RD));
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::CCParenthesis(double R) const {
  return standardParenthesis<O>(R) + factorCC * exp(-furtherExponentCC * radius<O>(R));
}

template<Utils::derivOrder O>
Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::SiOParenthesis(double R) const {
  auto r1 = radius<O>(R) - distanceSiO;
  return standardParenthesis<O>(R) + factorSiO * exp(-r1 * r1);
}

template<Utils::derivOrder O>
inline Utils::AutomaticDifferentiation::Value1DType<O> PM6PairwiseRepulsion::radius(double R) const {
  return Utils::AutomaticDifferentiation::variableWithUnitDerivative<O>(R);
}

template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::zero>
PM6PairwiseRepulsion::integral<Utils::derivOrder::zero>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  return igr;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::one>
PM6PairwiseRepulsion::integral<Utils::derivOrder::one>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  Utils::AutomaticDifferentiation::First1D integralD(igr, -igr * igr * igr * R);
  return integralD;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::two>
PM6PairwiseRepulsion::integral<Utils::derivOrder::two>(double R) const {
  double pSum = pA_.pCore() + pB_.pCore();
  double igr = 1.0 / std::sqrt(R * R + pSum * pSum);
  double igr3 = igr * igr * igr;
  Utils::AutomaticDifferentiation::Second1D integralD(igr, -igr3 * R, igr3 * igr * igr * (2 * R * R - pSum * pSum));
  return integralD;
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::first>
PM6PairwiseRepulsion::getDerivative<Utils::derivativeType::first>() const {
  return getRepulsionGradient();
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_atomic>
PM6PairwiseRepulsion::getDerivative<Utils::derivativeType::second_atomic>() const {
  return getRepulsionHessian();
}
template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_full>
PM6PairwiseRepulsion::getDerivative<Utils::derivativeType::second_full>() const {
  return getRepulsionHessian();
}

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // PAIRWISEREPULSION_H
