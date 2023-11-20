/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_RESOURCES_DFTB_PARAMETER_SETS_H
#define INCLUDE_SPARROW_RESOURCES_DFTB_PARAMETER_SETS_H

#include "Sparrow/Implementations/Dftb/ParameterSet.h"

namespace Scine {
namespace Sparrow {
namespace dftb {

boost::optional<ParameterSet> embeddedParameters(const std::string& name, const std::vector<int>& elements);

} // namespace dftb
} // namespace Sparrow
} // namespace Scine

#endif
