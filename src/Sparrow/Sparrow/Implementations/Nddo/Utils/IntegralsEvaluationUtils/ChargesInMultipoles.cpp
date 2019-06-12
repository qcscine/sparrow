/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ChargesInMultipoles.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

const std::vector<MultipoleCharge>& ChargesInMultipoles::getChargeConfiguration(multipole_t t) {
  static MultipoleChargesArray chargeConfigurations = createChargeConfigurations();
  return chargeConfigurations[static_cast<int>(t)];
}

ChargesInMultipoles::MultipoleChargesArray ChargesInMultipoles::createChargeConfigurations() {
  MultipoleChargesArray chargeConfigurations;

  chargeConfigurations[static_cast<int>(M00)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::d0, 1);

  chargeConfigurations[static_cast<int>(Qxx)].emplace_back(ChargeDistance::dm2, ChargeDistance::d0, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(Qxx)].emplace_back(ChargeDistance::dp2, ChargeDistance::d0, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(Qxx)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::d0, -0.50);

  chargeConfigurations[static_cast<int>(Qyy)].emplace_back(ChargeDistance::d0, ChargeDistance::dm2, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(Qyy)].emplace_back(ChargeDistance::d0, ChargeDistance::dp2, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(Qyy)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::d0, -0.50);

  chargeConfigurations[static_cast<int>(Qzz)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dm2, 0.25);
  chargeConfigurations[static_cast<int>(Qzz)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dp2, 0.25);
  chargeConfigurations[static_cast<int>(Qzz)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::d0, -0.50);

  chargeConfigurations[static_cast<int>(M2m2)].emplace_back(ChargeDistance::dp1, ChargeDistance::dp1, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(M2m2)].emplace_back(ChargeDistance::dm1, ChargeDistance::dm1, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(M2m2)].emplace_back(ChargeDistance::dp1, ChargeDistance::dm1, ChargeDistance::d0, -0.25);
  chargeConfigurations[static_cast<int>(M2m2)].emplace_back(ChargeDistance::dm1, ChargeDistance::dp1, ChargeDistance::d0, -0.25);

  chargeConfigurations[static_cast<int>(M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dp1, ChargeDistance::dp1, 0.25);
  chargeConfigurations[static_cast<int>(M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dm1, ChargeDistance::dm1, 0.25);
  chargeConfigurations[static_cast<int>(M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dp1, ChargeDistance::dm1, -0.25);
  chargeConfigurations[static_cast<int>(M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dm1, ChargeDistance::dp1, -0.25);

  chargeConfigurations[static_cast<int>(M21)].emplace_back(ChargeDistance::dp1, ChargeDistance::d0, ChargeDistance::dp1, 0.25);
  chargeConfigurations[static_cast<int>(M21)].emplace_back(ChargeDistance::dm1, ChargeDistance::d0, ChargeDistance::dm1, 0.25);
  chargeConfigurations[static_cast<int>(M21)].emplace_back(ChargeDistance::dp1, ChargeDistance::d0, ChargeDistance::dm1, -0.25);
  chargeConfigurations[static_cast<int>(M21)].emplace_back(ChargeDistance::dm1, ChargeDistance::d0, ChargeDistance::dp1, -0.25);

  chargeConfigurations[static_cast<int>(M20)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dps2, 0.250);
  chargeConfigurations[static_cast<int>(M20)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dms2, 0.250);
  chargeConfigurations[static_cast<int>(M20)].emplace_back(ChargeDistance::dps2, ChargeDistance::d0, ChargeDistance::d0, -0.375);
  chargeConfigurations[static_cast<int>(M20)].emplace_back(ChargeDistance::dms2, ChargeDistance::d0, ChargeDistance::d0, -0.375);
  chargeConfigurations[static_cast<int>(M20)].emplace_back(ChargeDistance::d0, ChargeDistance::dps2, ChargeDistance::d0, 0.125);
  chargeConfigurations[static_cast<int>(M20)].emplace_back(ChargeDistance::d0, ChargeDistance::dms2, ChargeDistance::d0, 0.125);

  chargeConfigurations[static_cast<int>(M22)].emplace_back(ChargeDistance::dps2, ChargeDistance::d0, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(M22)].emplace_back(ChargeDistance::dms2, ChargeDistance::d0, ChargeDistance::d0, 0.25);
  chargeConfigurations[static_cast<int>(M22)].emplace_back(ChargeDistance::d0, ChargeDistance::dps2, ChargeDistance::d0, -0.25);
  chargeConfigurations[static_cast<int>(M22)].emplace_back(ChargeDistance::d0, ChargeDistance::dms2, ChargeDistance::d0, -0.25);

  chargeConfigurations[static_cast<int>(Qzx)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dps2, 0.25);
  chargeConfigurations[static_cast<int>(Qzx)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dms2, 0.25);
  chargeConfigurations[static_cast<int>(Qzx)].emplace_back(ChargeDistance::dps2, ChargeDistance::d0, ChargeDistance::d0, -0.25);
  chargeConfigurations[static_cast<int>(Qzx)].emplace_back(ChargeDistance::dms2, ChargeDistance::d0, ChargeDistance::d0, -0.25);

  chargeConfigurations[static_cast<int>(M10)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dm1, -0.5);
  chargeConfigurations[static_cast<int>(M10)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::dp1, 0.5);
  chargeConfigurations[static_cast<int>(M11)].emplace_back(ChargeDistance::dm1, ChargeDistance::d0, ChargeDistance::d0, -0.5);
  chargeConfigurations[static_cast<int>(M11)].emplace_back(ChargeDistance::dp1, ChargeDistance::d0, ChargeDistance::d0, 0.5);
  chargeConfigurations[static_cast<int>(M1m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dm1, ChargeDistance::d0, -0.5);
  chargeConfigurations[static_cast<int>(M1m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dp1, ChargeDistance::d0, 0.5);

  return chargeConfigurations;
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
