/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTB3Initializer.h"
#include "DFTB3.h"
#include <Utils/IO/NativeFilenames.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

DFTB3Initializer::DFTB3Initializer(Utils::SinglePointMethod* method)
  : DFTBTypeInitializer("3ob-2-1"), dftb3_(*dynamic_cast<DFTB3*>(method)) {
}

void DFTB3Initializer::initialize(const Utils::ElementTypeCollection& elementTypes) {
  dftb3_.initializeStructure(elementTypes);
  auto path = Utils::NativeFilenames::combinePathSegments(getParameterFolder(), parameterSet_);
  dftb3_.initializeFromParameterPath(path);
}

void DFTB3Initializer::initialize(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes) {
  dftb3_.initializeStructure(elementTypes, positions);
  auto path = Utils::NativeFilenames::combinePathSegments(getParameterFolder(), parameterSet_);
  dftb3_.initializeFromParameterPath(path);
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
