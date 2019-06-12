/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DIPOLEMOMENTCALCULATOR_H
#define SPARROW_DIPOLEMOMENTCALCULATOR_H

#include <Utils/Typenames.h>

namespace Scine {
namespace Sparrow {

/**
 * @brief Interface for the calculation of the electrical dipole moment in a semiempirical method.
 */

class DipoleMomentCalculator {
 public:
  virtual Eigen::RowVector3d calculate() const = 0;
  virtual ~DipoleMomentCalculator() = default;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DIPOLEMOMENTCALCULATOR_H
