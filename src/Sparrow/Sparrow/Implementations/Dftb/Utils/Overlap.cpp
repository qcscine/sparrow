/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Overlap.h"
#include "ZeroOrderMatricesCalculator.h"

namespace Scine {
namespace Sparrow {

namespace dftb {

Overlap::Overlap(ZeroOrderMatricesCalculator& matricesCalculator) : matricesCalculator_(matricesCalculator) {
}

void Overlap::calculateOverlap(Utils::derivOrder highestRequiredOrder) {
  matricesCalculator_.calculateOverlap(highestRequiredOrder);
}

const Utils::MatrixWithDerivatives& Overlap::getOverlap() const {
  return matricesCalculator_.getOverlap();
}

void Overlap::reset() {
  matricesCalculator_.resetOverlap();
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
