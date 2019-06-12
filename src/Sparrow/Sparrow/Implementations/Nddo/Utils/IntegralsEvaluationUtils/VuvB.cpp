/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "VuvB.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

template<Utils::derivOrder O>
void VuvB::calculate(const Eigen::Vector3d& Rab, const ChargeSeparationParameter& D1, const KlopmanParameter& rho1,
                     double pcore, double ZB) {
  z_ = ZB;
  ChargeSeparationParameter dCore;
  KlopmanParameter rhoCore;
  rhoCore.set(multipolePair_t::ss0, pcore);

  dist1_ = D1;
  dist2_ = dCore;
  rho1_ = rho1;
  rho2_ = rhoCore;

  g_.calculate<O>(Rab);
}

template void VuvB::calculate<Utils::derivOrder::zero>(const Eigen::Vector3d&, const ChargeSeparationParameter&,
                                                       const KlopmanParameter&, double, double);
template void VuvB::calculate<Utils::derivOrder::one>(const Eigen::Vector3d&, const ChargeSeparationParameter&,
                                                      const KlopmanParameter&, double, double);
template void VuvB::calculate<Utils::derivOrder::two>(const Eigen::Vector3d&, const ChargeSeparationParameter&,
                                                      const KlopmanParameter&, double, double);

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
