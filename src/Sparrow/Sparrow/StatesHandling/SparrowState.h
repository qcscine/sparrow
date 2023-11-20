/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_SPARROWSTATE_H
#define SPARROW_SPARROWSTATE_H

#include <Core/BaseClasses/StateHandableObject.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <exception>

namespace Scine {
namespace Sparrow {

class GenericMethodWrapper;

class IncompatibleStateException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Sparrow state is not compatible with calculation (wrong number of electrons!).";
  }
};

/**
 * @brief Definition of a calculation state for methods implemented in Sparrow.
 * The calculation state is defined as the density matrix and the coefficient matrix at some point.
 * If the density matrix is empty, an exception is thrown at loading time.
 */
struct SparrowState final : public Core::State {
  explicit SparrowState(Utils::DensityMatrix densityMatrix) : densityMatrix_(std::move(densityMatrix)){};
  ~SparrowState() final = default;

  const Utils::DensityMatrix& getDensityMatrix() const {
    return densityMatrix_;
  }

 private:
  Utils::DensityMatrix densityMatrix_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_NDDOSTATE_H
