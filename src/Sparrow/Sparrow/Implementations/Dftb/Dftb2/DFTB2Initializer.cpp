/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTB2Initializer.h"
#include "DFTB2.h"
#include <Utils/IO/NativeFilenames.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

DFTB2Initializer::DFTB2Initializer(Utils::SinglePointMethod* method)
  : DFTBTypeInitializer("mio-1-1"), dftb2_(*dynamic_cast<DFTB2*>(method)) {
}

void DFTB2Initializer::initialize(const Utils::ElementTypeCollection& elementTypes) {
  dftb2_.initializeStructure(elementTypes);
  auto path = Utils::NativeFilenames::combinePathSegments(getParameterFolder(), parameterSet_);
  dftb2_.initializeFromParameterPath(path);
}

void DFTB2Initializer::initialize(const Utils::PositionCollection& positions, const Utils::ElementTypeCollection& elementTypes) {
  dftb2_.initializeStructure(elementTypes, positions);
  auto path = Utils::NativeFilenames::combinePathSegments(getParameterFolder(), parameterSet_);
  dftb2_.initializeFromParameterPath(path);
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
