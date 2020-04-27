/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "oneCenterIntegralContainer.h"
#include "oneCenterTwoElectronIntegrals.h"
#include <Utils/Geometry/ElementTypes.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

OneCenterIntegralContainer::OneCenterIntegralContainer() : matrices_(Container(110)) {
}

OneCenterIntegralContainer::~OneCenterIntegralContainer() = default;

void OneCenterIntegralContainer::clear() {
  std::fill(matrices_.begin(), matrices_.end(), nullptr);
}

void OneCenterIntegralContainer::set(Utils::ElementType e, integralMatrix_t mat) {
  matrices_[Utils::ElementInfo::Z(e)] = std::move(mat);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
