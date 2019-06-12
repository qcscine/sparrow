/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MultipoleChargePair.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

MultipoleChargePair::MultipoleChargePair(const MultipoleCharge& c1, const MultipoleCharge& c2) : c1_(c1), c2_(c2) {
  calculateCharges();
  calculateDistances();
}

void MultipoleChargePair::calculateCharges() {
  q_ = c1_.q * c2_.q;
}

void MultipoleChargePair::calculateDistances() {
  dx_ = calculateXYDistance(c1_.x, c2_.x);
  dy_ = calculateXYDistance(c1_.y, c2_.y);
  //  dz_ = calculateDistance(c1_.z, c2_.z);
}

ChargeDistanceSeparation MultipoleChargePair::calculateXYDistance(ChargeDistance d1, ChargeDistance d2) {
  if (d1 == ChargeDistance::d0 && d2 == ChargeDistance::d0)
    return ChargeDistanceSeparation::d00;

  if (d1 == ChargeDistance::d0 && (d2 == ChargeDistance::dp1 || d2 == ChargeDistance::dm1))
    return ChargeDistanceSeparation::d01;
  if (d1 == ChargeDistance::d0 && (d2 == ChargeDistance::dp2 || d2 == ChargeDistance::dm2))
    return ChargeDistanceSeparation::d02;
  if (d1 == ChargeDistance::d0 && (d2 == ChargeDistance::dps2 || d2 == ChargeDistance::dms2))
    return ChargeDistanceSeparation::d0s2;
  //  if(d2==ChargeDistance::d0 && (d1==ChargeDistance::dp1 || d1==ChargeDistance::dm1)) return
  //  ChargeDistanceSeparation::d10; if(d2==ChargeDistance::d0 && (d1==ChargeDistance::dp2 || d1==ChargeDistance::dm2))
  //  return ChargeDistanceSeparation::d20; if(d2==ChargeDistance::d0 && (d1==ChargeDistance::dps2 ||
  //  d1==ChargeDistance::dms2)) return ChargeDistanceSeparation::ds20;
  if (d2 == ChargeDistance::d0 && (d1 == ChargeDistance::dp1 || d1 == ChargeDistance::dm1))
    return ChargeDistanceSeparation::d10;
  if (d2 == ChargeDistance::d0 && (d1 == ChargeDistance::dp2 || d1 == ChargeDistance::dm2))
    return ChargeDistanceSeparation::d20;
  if (d2 == ChargeDistance::d0 && (d1 == ChargeDistance::dps2 || d1 == ChargeDistance::dms2))
    return ChargeDistanceSeparation::ds20;

  if ((d1 == ChargeDistance::dp1 && d2 == ChargeDistance::dp1) || (d1 == ChargeDistance::dm1 && d2 == ChargeDistance::dm1))
    return ChargeDistanceSeparation::m11;
  if ((d1 == ChargeDistance::dp2 && d2 == ChargeDistance::dp2) || (d1 == ChargeDistance::dm2 && d2 == ChargeDistance::dm2))
    return ChargeDistanceSeparation::m22;
  if ((d1 == ChargeDistance::dps2 && d2 == ChargeDistance::dps2) || (d1 == ChargeDistance::dms2 && d2 == ChargeDistance::dms2))
    return ChargeDistanceSeparation::ms2s2;
  if ((d1 == ChargeDistance::dp1 && d2 == ChargeDistance::dm1) || (d1 == ChargeDistance::dm1 && d2 == ChargeDistance::dp1))
    return ChargeDistanceSeparation::p11;
  if ((d1 == ChargeDistance::dp2 && d2 == ChargeDistance::dm2) || (d1 == ChargeDistance::dm2 && d2 == ChargeDistance::dp2))
    return ChargeDistanceSeparation::p22;
  if ((d1 == ChargeDistance::dps2 && d2 == ChargeDistance::dms2) || (d1 == ChargeDistance::dms2 && d2 == ChargeDistance::dps2))
    return ChargeDistanceSeparation::ps2s2;

  if ((d1 == ChargeDistance::dp1 && d2 == ChargeDistance::dp2) || (d1 == ChargeDistance::dm1 && d2 == ChargeDistance::dm2))
    return ChargeDistanceSeparation::m12;
  if ((d1 == ChargeDistance::dp1 && d2 == ChargeDistance::dps2) || (d1 == ChargeDistance::dm1 && d2 == ChargeDistance::dms2))
    return ChargeDistanceSeparation::m1s2;
  if ((d1 == ChargeDistance::dp2 && d2 == ChargeDistance::dp1) || (d1 == ChargeDistance::dm2 && d2 == ChargeDistance::dm1))
    return ChargeDistanceSeparation::m21;
  if ((d1 == ChargeDistance::dps2 && d2 == ChargeDistance::dp1) || (d1 == ChargeDistance::dms2 && d2 == ChargeDistance::dm1))
    return ChargeDistanceSeparation::ms21;

  if ((d1 == ChargeDistance::dp1 && d2 == ChargeDistance::dm2) || (d1 == ChargeDistance::dm1 && d2 == ChargeDistance::dp2))
    return ChargeDistanceSeparation::p12;
  if ((d1 == ChargeDistance::dp1 && d2 == ChargeDistance::dms2) || (d1 == ChargeDistance::dm1 && d2 == ChargeDistance::dps2))
    return ChargeDistanceSeparation::p1s2;
  if ((d1 == ChargeDistance::dp2 && d2 == ChargeDistance::dm1) || (d1 == ChargeDistance::dm2 && d2 == ChargeDistance::dp1))
    return ChargeDistanceSeparation::p21;
  if ((d1 == ChargeDistance::dps2 && d2 == ChargeDistance::dm1) || (d1 == ChargeDistance::dms2 && d2 == ChargeDistance::dp1))
    return ChargeDistanceSeparation::ps21;

  if ((d1 == ChargeDistance::dp2 && d2 == ChargeDistance::dps2) || (d1 == ChargeDistance::dm2 && d2 == ChargeDistance::dms2))
    return ChargeDistanceSeparation::m2s2;
  if ((d1 == ChargeDistance::dps2 && d2 == ChargeDistance::dp2) || (d1 == ChargeDistance::dms2 && d2 == ChargeDistance::dm2))
    return ChargeDistanceSeparation::ms22;
  if ((d1 == ChargeDistance::dp2 && d2 == ChargeDistance::dms2) || (d1 == ChargeDistance::dm2 && d2 == ChargeDistance::dps2))
    return ChargeDistanceSeparation::p2s2;
  if ((d1 == ChargeDistance::dps2 && d2 == ChargeDistance::dm2) || (d1 == ChargeDistance::dms2 && d2 == ChargeDistance::dp2))
    return ChargeDistanceSeparation::ps22;

  return ChargeDistanceSeparation::d00; // Shouldn't come here
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
