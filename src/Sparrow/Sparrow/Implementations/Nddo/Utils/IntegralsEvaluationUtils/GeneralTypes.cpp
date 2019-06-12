/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GeneralTypes.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace GeneralTypes {

int orbitalQN(orb_t o) {
  switch (o) {
    case orb_t::s:
      return 0;
    case orb_t::x:
      return 1;
    case orb_t::y:
      return 1;
    case orb_t::z:
      return 1;
    case orb_t::z2:
      return 2;
    case orb_t::xz:
      return 2;
    case orb_t::yz:
      return 2;
    case orb_t::x2y2:
      return 2;
    case orb_t::xy:
      return 2;
    default:
      return 0;
  }
}

std::pair<orb_t, orb_t> separatePair(twoElIntegral_t t) {
  orb_t l1, l2;
  switch (t) {
    case twoElIntegral_t::s_s:
      l1 = orb_t::s;
      l2 = orb_t::s;
      break;
    case twoElIntegral_t::s_z:
      l1 = orb_t::s;
      l2 = orb_t::z;
      break;
    case twoElIntegral_t::s_x:
      l1 = orb_t::s;
      l2 = orb_t::x;
      break;
    case twoElIntegral_t::s_y:
      l1 = orb_t::s;
      l2 = orb_t::y;
      break;
    case twoElIntegral_t::z_z:
      l1 = orb_t::z;
      l2 = orb_t::z;
      break;
    case twoElIntegral_t::x_z:
      l1 = orb_t::x;
      l2 = orb_t::z;
      break;
    case twoElIntegral_t::y_z:
      l1 = orb_t::y;
      l2 = orb_t::z;
      break;
    case twoElIntegral_t::x_x:
      l1 = orb_t::x;
      l2 = orb_t::x;
      break;
    case twoElIntegral_t::x_y:
      l1 = orb_t::x;
      l2 = orb_t::y;
      break;
    case twoElIntegral_t::y_y:
      l1 = orb_t::y;
      l2 = orb_t::y;
      break;
    case twoElIntegral_t::s_z2:
      l1 = orb_t::s;
      l2 = orb_t::z2;
      break;
    case twoElIntegral_t::s_xz:
      l1 = orb_t::s;
      l2 = orb_t::xz;
      break;
    case twoElIntegral_t::s_yz:
      l1 = orb_t::s;
      l2 = orb_t::yz;
      break;
    case twoElIntegral_t::s_x2y2:
      l1 = orb_t::s;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::s_xy:
      l1 = orb_t::s;
      l2 = orb_t::xy;
      break;
    case twoElIntegral_t::z_z2:
      l1 = orb_t::z;
      l2 = orb_t::z2;
      break;
    case twoElIntegral_t::z_xz:
      l1 = orb_t::z;
      l2 = orb_t::xz;
      break;
    case twoElIntegral_t::z_yz:
      l1 = orb_t::z;
      l2 = orb_t::yz;
      break;
    case twoElIntegral_t::x_z2:
      l1 = orb_t::x;
      l2 = orb_t::z2;
      break;
    case twoElIntegral_t::x_xz:
      l1 = orb_t::x;
      l2 = orb_t::xz;
      break;
    case twoElIntegral_t::x_x2y2:
      l1 = orb_t::x;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::x_xy:
      l1 = orb_t::x;
      l2 = orb_t::xy;
      break;
    case twoElIntegral_t::y_z2:
      l1 = orb_t::y;
      l2 = orb_t::z2;
      break;
    case twoElIntegral_t::y_yz:
      l1 = orb_t::y;
      l2 = orb_t::yz;
      break;
    case twoElIntegral_t::y_x2y2:
      l1 = orb_t::y;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::y_xy:
      l1 = orb_t::y;
      l2 = orb_t::xy;
      break;
    case twoElIntegral_t::z2_z2:
      l1 = orb_t::z2;
      l2 = orb_t::z2;
      break;
    case twoElIntegral_t::z2_xz:
      l1 = orb_t::z2;
      l2 = orb_t::xz;
      break;
    case twoElIntegral_t::z2_yz:
      l1 = orb_t::z2;
      l2 = orb_t::yz;
      break;
    case twoElIntegral_t::z2_x2y2:
      l1 = orb_t::z2;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::z2_xy:
      l1 = orb_t::z2;
      l2 = orb_t::xy;
      break;
    case twoElIntegral_t::xz_xz:
      l1 = orb_t::xz;
      l2 = orb_t::xz;
      break;
    case twoElIntegral_t::xz_yz:
      l1 = orb_t::xz;
      l2 = orb_t::yz;
      break;
    case twoElIntegral_t::xz_x2y2:
      l1 = orb_t::xz;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::xz_xy:
      l1 = orb_t::xz;
      l2 = orb_t::xy;
      break;
    case twoElIntegral_t::yz_yz:
      l1 = orb_t::yz;
      l2 = orb_t::yz;
      break;
    case twoElIntegral_t::yz_x2y2:
      l1 = orb_t::yz;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::yz_xy:
      l1 = orb_t::yz;
      l2 = orb_t::xy;
      break;
    case twoElIntegral_t::x2y2_x2y2:
      l1 = orb_t::x2y2;
      l2 = orb_t::x2y2;
      break;
    case twoElIntegral_t::xy_xy:
      l1 = orb_t::xy;
      l2 = orb_t::xy;
      break;
    default:
      throw InvalidOrbitalPairException();
  }
  return std::make_pair(l1, l2);
}

} // namespace GeneralTypes

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
