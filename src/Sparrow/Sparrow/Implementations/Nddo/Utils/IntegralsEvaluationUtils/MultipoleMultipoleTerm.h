/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef NDDO_MULTIPOLE_MULTIPOLEMULTIPOLETERM_H
#define NDDO_MULTIPOLE_MULTIPOLEMULTIPOLETERM_H

#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/**
 * @brief This header-only class defines an object for the calculation of an interaction
 *        between two charges in a multipole.
 *
 * The total interaction between two electrons is approximated by a classical multipole expansion. This class
 * defines an object that calculates the interaction between two charges of said multipoles.
 * The charge interaction is calculated within the Klopman approximation to retrieve the correct 1-center, 2-electrons
 * interaction in the limit of vanishing distance between the charges.
 * The template functions allow for the analytical calculation of all the values up to the second derivative of the
 * repulsion energy.
 */

class MultipoleMultipoleTerm {
 public:
  /**
   * @brief constructor, sets the coordinates of the two point charges and the value of the two charges multiplied
   * @param f the result of the multiplication of the two charges
   * @param x1 x configuration relative to the first multipole center, a scaling factor for the charge separations
   * @param y1 y configuration relative to the first multipole center, a scaling factor for the charge separations
   * @param z1 z configuration relative to the first multipole center, a scaling factor for the charge separations
   * @param x2 x configuration relative to the second multipole center, a scaling factor for the charge separations
   * @param y2 y configuration relative to the second multipole center, a scaling factor for the charge separations
   * @param z2 z configuration relative to the second multipole center, a scaling factor for the charge separations
   */
  MultipoleMultipoleTerm(double f, double x1, double y1, double z1, double x2, double y2, double z2)
    : f_(f), x1_(x1), x2_(x2), y1_(y1), y2_(y2), z1_(z1), z2_(z2) {
  }
  /**
   * @brief calculates the interaction energy between two charges up to a certain derivative
   * @tparam O the required derivative order for the interaction energy
   * @param R the distance between the two multipole centers
   * @param D1 the charge separation parameter for the first multipole
   * @param d2 the charge separation parameter for the second multipole
   * @param squaredRhos the squared Klopman--Ohno parameters
   * @return an Utils::AutomaticDifferentiation::Value1DType<O> object, storing the value of the interaction up to the
   *         O-th derivative
   */
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> calculate(double R, double D1, double d2, double squaredRhos) const {
    double dx = d2 * x2_ - D1 * x1_;
    double dy = d2 * y2_ - D1 * y1_;
    double dz = R + d2 * z2_ - D1 * z1_;
    double invsqrt = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz + squaredRhos);
    return expr<O>(f_, dz, invsqrt);
  }
  /**
   * @brief returns the values up to the O-th derivative order with respect to the inter-multipole distance
   * @tparam O the required derivative order
   * @param f the result of the multiplication of the two charges
   * @param dz the distance between the two charges along the inter-multipole axis, for the 0-th derivative
   *        no dz argument is needed.
   * @param invsqrt the denominator of the interaction energy
   * @return an Utils::AutomaticDifferentiation::Value1DType<O> object containing the derivatives of the interaction
   * energy up to the O-th order
   */
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> expr(double f, double dz, double invsqrt) const;

 private:
  double f_;
  double x1_, x2_, y1_, y2_, z1_, z2_;
};

template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::zero>
MultipoleMultipoleTerm::expr<Utils::derivOrder::zero>(double f, double /*dz*/, double invsqrt) const {
  return f * invsqrt;
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::one>
MultipoleMultipoleTerm::expr<Utils::derivOrder::one>(double f, double dz, double invsqrt) const {
  return {f * invsqrt, -dz * f * invsqrt * invsqrt * invsqrt};
}
template<>
inline Utils::AutomaticDifferentiation::Value1DType<Utils::derivOrder::two>
MultipoleMultipoleTerm::expr<Utils::derivOrder::two>(double f, double dz, double invsqrt) const {
  double inv3 = invsqrt * invsqrt * invsqrt;
  return {f * invsqrt, -dz * f * inv3, f * inv3 * (3 * dz * dz * invsqrt * invsqrt - 1)};
}

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // NDDO_MULTIPOLE_MULTIPOLEMULTIPOLETERM_H
