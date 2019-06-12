/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ELEMENTPAIRPARAMETERS_H
#define SPARROW_ELEMENTPAIRPARAMETERS_H

#include "PM6DiatomicParameters.h"
#include <Utils/Geometry/ElementTypes.h>
#include <Utils/MethodEssentials/util/MethodExceptions.h>
#include <array>
#include <cassert>
#include <memory>

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * This class holds the all the pointers to the parameters for
 * element pairs.
 */
class ElementPairParameters {
 public:
  using PM6DiatomicParametersPtr = std::unique_ptr<PM6DiatomicParameters>;
  using ParameterContainer = std::array<PM6DiatomicParametersPtr, 110>;
  using PPContainer = std::array<ParameterContainer, 110>;

  /*! Delete all parameters. */
  void clear() {
    parameters_ = PPContainer{};
  }

  bool isSet(Utils::ElementType e1, Utils::ElementType e2) const {
    if (e1 < e2)
      std::swap(e1, e2);
    return parameters_[static_cast<int>(e1)][static_cast<int>(e2)] != nullptr;
  }
  void set(Utils::ElementType e1, Utils::ElementType e2, std::unique_ptr<PM6DiatomicParameters>&& parameters) {
    if (e1 < e2)
      std::swap(e1, e2);
    parameters_[static_cast<int>(e1)][static_cast<int>(e2)] = std::move(parameters);
  }
  const PM6DiatomicParameters& get(Utils::ElementType e1, Utils::ElementType e2) const {
    if (isSet(e1, e2)) {
      if (e1 < e2)
        std::swap(e1, e2);
      return *parameters_[static_cast<int>(e1)][static_cast<int>(e2)];
    }
    throw Utils::Methods::ParametersDoNotExistForElementPairException(e1, e2);
  }

 private:
  PPContainer parameters_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ELEMENTPAIRPARAMETERS_H
