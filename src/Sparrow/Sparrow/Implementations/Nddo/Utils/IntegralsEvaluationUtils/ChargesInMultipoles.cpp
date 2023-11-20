/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ChargesInMultipoles.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

template<typename T>
auto underlying(T t) {
  return static_cast<typename std::underlying_type<T>::type>(t);
}

const std::vector<MultipoleCharge>& ChargesInMultipoles::getChargeConfiguration(Multipole t) {
  static MultipoleChargesArray chargeConfigurations = createChargeConfigurations();
  return chargeConfigurations[underlying(t)];
}

ChargesInMultipoles::MultipoleChargesArray ChargesInMultipoles::createChargeConfigurations() {
  MultipoleChargesArray chargeConfigurations;

  chargeConfigurations[underlying(Multipole::M00)].emplace_back(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::d0, 1);

  chargeConfigurations[underlying(Multipole::Qxx)].emplace_back(ChargeDistance::dm2, ChargeDistance::d0,
                                                                ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::Qxx)].emplace_back(ChargeDistance::dp2, ChargeDistance::d0,
                                                                ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::Qxx)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.50);

  chargeConfigurations[underlying(Multipole::Qyy)].emplace_back(ChargeDistance::d0, ChargeDistance::dm2,
                                                                ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::Qyy)].emplace_back(ChargeDistance::d0, ChargeDistance::dp2,
                                                                ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::Qyy)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.50);

  chargeConfigurations[underlying(Multipole::Qzz)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dm2, 0.25);
  chargeConfigurations[underlying(Multipole::Qzz)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dp2, 0.25);
  chargeConfigurations[underlying(Multipole::Qzz)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.50);

  chargeConfigurations[underlying(Multipole::M2m2)].emplace_back(ChargeDistance::dp1, ChargeDistance::dp1,
                                                                 ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::M2m2)].emplace_back(ChargeDistance::dm1, ChargeDistance::dm1,
                                                                 ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::M2m2)].emplace_back(ChargeDistance::dp1, ChargeDistance::dm1,
                                                                 ChargeDistance::d0, -0.25);
  chargeConfigurations[underlying(Multipole::M2m2)].emplace_back(ChargeDistance::dm1, ChargeDistance::dp1,
                                                                 ChargeDistance::d0, -0.25);

  chargeConfigurations[underlying(Multipole::M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dp1,
                                                                 ChargeDistance::dp1, 0.25);
  chargeConfigurations[underlying(Multipole::M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dm1,
                                                                 ChargeDistance::dm1, 0.25);
  chargeConfigurations[underlying(Multipole::M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dp1,
                                                                 ChargeDistance::dm1, -0.25);
  chargeConfigurations[underlying(Multipole::M2m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dm1,
                                                                 ChargeDistance::dp1, -0.25);

  chargeConfigurations[underlying(Multipole::M21)].emplace_back(ChargeDistance::dp1, ChargeDistance::d0,
                                                                ChargeDistance::dp1, 0.25);
  chargeConfigurations[underlying(Multipole::M21)].emplace_back(ChargeDistance::dm1, ChargeDistance::d0,
                                                                ChargeDistance::dm1, 0.25);
  chargeConfigurations[underlying(Multipole::M21)].emplace_back(ChargeDistance::dp1, ChargeDistance::d0,
                                                                ChargeDistance::dm1, -0.25);
  chargeConfigurations[underlying(Multipole::M21)].emplace_back(ChargeDistance::dm1, ChargeDistance::d0,
                                                                ChargeDistance::dp1, -0.25);

  chargeConfigurations[underlying(Multipole::M20)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dps2, 0.250);
  chargeConfigurations[underlying(Multipole::M20)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dms2, 0.250);
  chargeConfigurations[underlying(Multipole::M20)].emplace_back(ChargeDistance::dps2, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.375);
  chargeConfigurations[underlying(Multipole::M20)].emplace_back(ChargeDistance::dms2, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.375);
  chargeConfigurations[underlying(Multipole::M20)].emplace_back(ChargeDistance::d0, ChargeDistance::dps2,
                                                                ChargeDistance::d0, 0.125);
  chargeConfigurations[underlying(Multipole::M20)].emplace_back(ChargeDistance::d0, ChargeDistance::dms2,
                                                                ChargeDistance::d0, 0.125);

  chargeConfigurations[underlying(Multipole::M22)].emplace_back(ChargeDistance::dps2, ChargeDistance::d0,
                                                                ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::M22)].emplace_back(ChargeDistance::dms2, ChargeDistance::d0,
                                                                ChargeDistance::d0, 0.25);
  chargeConfigurations[underlying(Multipole::M22)].emplace_back(ChargeDistance::d0, ChargeDistance::dps2,
                                                                ChargeDistance::d0, -0.25);
  chargeConfigurations[underlying(Multipole::M22)].emplace_back(ChargeDistance::d0, ChargeDistance::dms2,
                                                                ChargeDistance::d0, -0.25);

  chargeConfigurations[underlying(Multipole::Qzx)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dps2, 0.25);
  chargeConfigurations[underlying(Multipole::Qzx)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dms2, 0.25);
  chargeConfigurations[underlying(Multipole::Qzx)].emplace_back(ChargeDistance::dps2, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.25);
  chargeConfigurations[underlying(Multipole::Qzx)].emplace_back(ChargeDistance::dms2, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.25);

  chargeConfigurations[underlying(Multipole::M10)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dm1, -0.5);
  chargeConfigurations[underlying(Multipole::M10)].emplace_back(ChargeDistance::d0, ChargeDistance::d0,
                                                                ChargeDistance::dp1, 0.5);
  chargeConfigurations[underlying(Multipole::M11)].emplace_back(ChargeDistance::dm1, ChargeDistance::d0,
                                                                ChargeDistance::d0, -0.5);
  chargeConfigurations[underlying(Multipole::M11)].emplace_back(ChargeDistance::dp1, ChargeDistance::d0,
                                                                ChargeDistance::d0, 0.5);
  chargeConfigurations[underlying(Multipole::M1m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dm1,
                                                                 ChargeDistance::d0, -0.5);
  chargeConfigurations[underlying(Multipole::M1m1)].emplace_back(ChargeDistance::d0, ChargeDistance::dp1,
                                                                 ChargeDistance::d0, 0.5);

  return chargeConfigurations;
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
