/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MULTIPOLECHARGEPAIR_H
#define SPARROW_MULTIPOLECHARGEPAIR_H

#include "MultipoleCharge.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

// d=dist, p=plus, m=minus
enum class ChargeDistanceSeparation {
  d00,
  d01,
  d10,
  d02,
  d20,
  d0s2,
  ds20,
  p11,
  m11,
  p12,
  p21,
  m12,
  m21,
  p1s2,
  ps21,
  m1s2,
  ms21,
  p22,
  m22,
  p2s2,
  ps22,
  m2s2,
  ms22,
  ps2s2,
  ms2s2
};

/*!
 * This class stores the information about the distance between two point charges
 * and about the product of their charges
 */

class MultipoleChargePair {
 public:
  MultipoleChargePair(const MultipoleCharge& c1, const MultipoleCharge& c2);

  ChargeDistanceSeparation getXDistance() const {
    return dx_;
  }
  ChargeDistanceSeparation getYDistance() const {
    return dy_;
  }
  // ChargeDistanceSeparation getZDistance() const { return dz_; }
  void setChargeProduct(double q) {
    q_ = q;
  }
  double getChargeProduct() const {
    return q_;
  }
  const MultipoleCharge& firstCharge() const {
    return c1_;
  }
  const MultipoleCharge& secondCharge() const {
    return c2_;
  }

 private:
  void calculateCharges();
  void calculateDistances();
  ChargeDistanceSeparation calculateXYDistance(ChargeDistance d1, ChargeDistance d2);

  MultipoleCharge c1_, c2_;
  ChargeDistanceSeparation dx_, dy_; //,dz_;
  double q_;
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MULTIPOLECHARGEPAIR_H
