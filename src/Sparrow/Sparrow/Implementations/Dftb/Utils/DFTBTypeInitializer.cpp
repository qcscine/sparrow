/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTBTypeInitializer.h"
#include <Utils/IO/NativeFilenames.h>
#include <utility>

namespace Scine {
namespace Sparrow {

namespace dftb {

DFTBTypeInitializer::DFTBTypeInitializer(std::string set) : parameterSet_(std::move(set)) {
}

std::string DFTBTypeInitializer::getParameterFolder() const {
  return Utils::NativeFilenames::combinePathSegments(rootParameterFolder_, "dftb");
}

void DFTBTypeInitializer::setParameterSet(std::string set) {
  parameterSet_ = std::move(set);
}

std::string DFTBTypeInitializer::getParameterSet() const {
  return parameterSet_;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
