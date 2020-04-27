/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ELEMENTPARAMETERS_H
#define SPARROW_ELEMENTPARAMETERS_H

#include "AtomicParameters.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/ElementTypes.h>
#include <memory>

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * This class holds the pointers to the element-specific
 * runtime parameters
 */

class ElementParameters {
 public:
  using par_t = std::unique_ptr<AtomicParameters>;
  using PPContainer = std::vector<par_t>;

  ElementParameters() : parameters_(110) {
  }
  void clear() {
    std::fill(parameters_.begin(), parameters_.end(), nullptr);
  }
  bool isSet(Utils::ElementType e) const {
    return parameters_[Utils::ElementInfo::Z(e)] != nullptr;
  }
  void set(Utils::ElementType e, par_t parameters) {
    parameters_[Utils::ElementInfo::Z(e)] = std::move(parameters);
  } // TODO: elementtype not needed, included in parameters
  const AtomicParameters& get(Utils::ElementType e) const {
    return *parameters_[Utils::ElementInfo::Z(e)];
  }

 private:
  PPContainer parameters_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ELEMENTPARAMETERS_H
