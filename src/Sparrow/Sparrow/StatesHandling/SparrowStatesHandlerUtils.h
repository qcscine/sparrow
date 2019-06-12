/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_SPARROWSTATESHANDLERUTILS_H
#define SPARROW_SPARROWSTATESHANDLERUTILS_H

#include <exception>

namespace Scine {
namespace Utils {
class DensityMatrix;
class LCAOMethod;
} // namespace Utils
namespace Sparrow {
class SparrowState;

class StateNotCompatibleException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "State is not compatible with SparrowState.";
  }
};

class SparrowStatesHandlerUtils {
 public:
  static void loadDensityMatrix(SparrowState& state, Utils::LCAOMethod& method);
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_SPARROWSTATESHANDLERUTILS_H
