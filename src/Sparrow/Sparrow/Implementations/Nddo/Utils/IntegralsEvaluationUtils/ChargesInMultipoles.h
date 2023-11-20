/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_CHARGESINMULTIPOLES_H
#define SPARROW_CHARGESINMULTIPOLES_H

#include "MultipoleCharge.h"
#include "multipoleTypes.h"
#include <array>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/*!
 * This class defines the point charges of the different multipole.
 */

class ChargesInMultipoles {
 public:
  using MultipoleCharges = std::vector<MultipoleCharge>;
  using MultipoleChargesArray = std::array<MultipoleCharges, 13>;
  static_assert(static_cast<std::underlying_type<Multipole>::type>(Multipole::Qzx) == 12,
                "multipole type enum layout has changed");

  static const std::vector<MultipoleCharge>& getChargeConfiguration(Multipole t);

 private:
  static MultipoleChargesArray createChargeConfigurations();
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_CHARGESINMULTIPOLES_H
