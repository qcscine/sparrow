/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_VUVB_H
#define SPARROW_VUVB_H

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/Global2c2eMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/*!
 * This class returns the \f$V_{\mu\nu,B}\f$ terms needed in
 * semi-empirical methods.
 */

class VuvB {
 public:
  using orb_index_t = Global2c2eMatrix::orb_index_t;
  using orbPair_index_t = Global2c2eMatrix::orbPair_index_t;

  explicit VuvB(int l1) : g_(l1, 0, dist1_, dist2_, rho1_, rho2_) {
  }
  template<Utils::derivOrder O>
  void calculate(const Eigen::Vector3d& Rab, const ChargeSeparationParameter& D1, const KlopmanParameter& rho1,
                 double pcore, double ZB);
  double get(orb_index_t o1, orb_index_t o2) const {
    return -z_ * g_.get(o1, o2, static_cast<int>(GeneralTypes::orb_t::s), static_cast<int>(GeneralTypes::orb_t::s));
  }
  double get(orbPair_index_t op1) const {
    return -z_ * g_.get(op1, static_cast<int>(GeneralTypes::twoElIntegral_t::s_s));
  }
  template<Utils::derivativeType O>
  Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative(orb_index_t o1, orb_index_t o2) const {
    return -z_ * g_.getDerivative<O>(o1, o2, static_cast<int>(GeneralTypes::orb_t::s),
                                     static_cast<int>(GeneralTypes::orb_t::s));
  }

 private:
  ChargeSeparationParameter dist1_, dist2_;
  KlopmanParameter rho1_, rho2_;
  Global2c2eMatrix g_;
  double z_{0.0};
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_VUVB_H
