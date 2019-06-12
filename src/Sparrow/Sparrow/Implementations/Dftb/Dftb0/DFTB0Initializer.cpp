/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTB0Initializer.h"
#include "DFTB0.h"
#include <Utils/IO/NativeFilenames.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

DFTB0Initializer::DFTB0Initializer(Utils::SinglePointMethod* method)
  : DFTBTypeInitializer("3ob-2-1"), dftb0_(*dynamic_cast<DFTB0*>(method)) {
}

void DFTB0Initializer::initialize(const Utils::ElementTypeCollection& elementTypes) {
  dftb0_.initializeStructure(elementTypes);
  dftb0_.initializeFromParameterPath(Utils::NativeFilenames::combinePathSegments(getParameterFolder(), parameterSet_));
}

void DFTB0Initializer::initialize(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes) {
  dftb0_.initializeStructure(elementTypes, positions);
  dftb0_.initializeFromParameterPath(Utils::NativeFilenames::combinePathSegments(getParameterFolder(), parameterSet_));
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
