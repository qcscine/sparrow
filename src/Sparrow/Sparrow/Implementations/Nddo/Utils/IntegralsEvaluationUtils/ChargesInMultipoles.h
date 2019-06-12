/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

  static const std::vector<MultipoleCharge>& getChargeConfiguration(multipole_t t);

 private:
  static MultipoleChargesArray createChargeConfigurations();
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_CHARGESINMULTIPOLES_H
