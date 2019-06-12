/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MMTermCreator.h"
#include "MultipoleCharge.h"
#include <cmath>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

void MMTermCreator::add(const MultipoleChargePair& cp) {
  for (auto& t : terms_) {
    if (equivalent(t, cp)) {
      t.setChargeProduct(t.getChargeProduct() + cp.getChargeProduct());
      return;
    }
  }

  // was not already in list
  terms_.push_back(cp);
}

MultipoleMultipoleInteraction MMTermCreator::computeList() const {
  MultipoleMultipoleInteraction m;

  for (const auto& t : terms_) {
    if (std::abs(t.getChargeProduct()) < 1e-6)
      continue; // "==0" would be a bit risky

    MultipoleCharge c1 = t.firstCharge();
    MultipoleCharge c2 = t.secondCharge();
    MultipoleMultipoleTerm term(t.getChargeProduct(), fromEnum(c1.x), fromEnum(c1.y), fromEnum(c1.z), fromEnum(c2.x),
                                fromEnum(c2.y), fromEnum(c2.z));
    m.add(term);
  }

  return m;
}

bool MMTermCreator::equivalent(const MultipoleChargePair& p1, const MultipoleChargePair& p2) {
  if ((p1.firstCharge().z != p2.firstCharge().z) || (p1.secondCharge().z != p2.secondCharge().z))
    return false;

  if (p1.getXDistance() == p2.getXDistance() && p1.getYDistance() == p2.getYDistance())
    return true;
  if (p1.getXDistance() == p2.getYDistance() && p1.getYDistance() == p2.getXDistance())
    return true;

  return false;
}

double MMTermCreator::fromEnum(ChargeDistance d) const {
  switch (d) {
    case ChargeDistance::d0:
      return 0.0;
    case ChargeDistance::dp1:
      return 1.0;
    case ChargeDistance::dm1:
      return -1.0;
    case ChargeDistance::dp2:
      return 2.0;
    case ChargeDistance::dm2:
      return -2.0;
    case ChargeDistance::dps2:
      return sqrt(2.0);
    case ChargeDistance::dms2:
      return -sqrt(2.0);
    default:
      return 0.0;
  }
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
