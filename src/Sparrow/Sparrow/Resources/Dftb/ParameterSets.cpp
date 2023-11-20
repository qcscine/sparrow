/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Sparrow/Resources/Dftb/ParameterSets.h"
#include "Sparrow/Resources/Dftb/3ob-2-1/Parameters.h"
#include "Sparrow/Resources/Dftb/3ob-3-1/Parameters.h"
#include "Sparrow/Resources/Dftb/borg-0-1/Parameters.h"
#include "Sparrow/Resources/Dftb/hyb-0-2/Parameters.h"
#include "Sparrow/Resources/Dftb/mio-0-1/Parameters.h"
#include "Sparrow/Resources/Dftb/mio-1-1/Parameters.h"
#include "Sparrow/Resources/Dftb/trans3d-0-1/Parameters.h"
#include "Sparrow/Resources/Dftb/znorg-0-1/Parameters.h"

namespace Scine {
namespace Sparrow {
namespace dftb {

boost::optional<ParameterSet> embeddedParameters(const std::string& name, const std::vector<int>& elements) {
  if (name == "3ob-3-1") {
    return params_3ob_3_1(elements);
  }

  if (name == "3ob-2-1") {
    return params_3ob_2_1(elements);
  }

  if (name == "mio-1-1") {
    return params_mio_1_1(elements);
  }

  if (name == "mio-0-1") {
    return params_mio_0_1(elements);
  }

  if (name == "znorg-0-1") {
    return params_mio_0_1(elements).patch(params_znorg_0_1(elements));
  }

  if (name == "hyb-0-2") {
    return params_mio_0_1(elements).patch(params_hyb_0_2(elements));
  }

  if (name == "trans3d-0-1") {
    return params_mio_1_1(elements).patch(params_trans3d_0_1(elements));
  }

  if (name == "borg-0-1") {
    return params_mio_1_1(elements).patch(params_borg_0_1(elements));
  }

  return boost::none;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
