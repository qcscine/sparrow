/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MULTIPOLECHARGE_H
#define SPARROW_MULTIPOLECHARGE_H

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

// All multipliers for the charge separation parameters possible: 0, +-1, +-2, +-sqrt(2)
enum class ChargeDistance { d0, dp1, dm1, dp2, dm2, dps2, dms2 };

/*!
 * This class defines an object containing the position and charge of a point charge.
 */
class MultipoleCharge {
 public:
  MultipoleCharge() = default;

  MultipoleCharge(ChargeDistance dx, ChargeDistance dy, ChargeDistance dz, double c) : x(dx), y(dy), z(dz), q(c) {
  }

  ChargeDistance x{ChargeDistance::d0}, y{ChargeDistance::d0}, z{ChargeDistance::d0};
  double q{0};
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MULTIPOLECHARGE_H
