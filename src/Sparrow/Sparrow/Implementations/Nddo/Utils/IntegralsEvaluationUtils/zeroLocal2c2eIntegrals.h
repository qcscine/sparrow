/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ZEROLOCAL2C2EINTEGRALS_H
#define SPARROW_ZEROLOCAL2C2EINTEGRALS_H

#include "multipoleTypes.h"
#include <array>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/**
 * @brief Class that specifies which local two-center two-electron integrals are equal to zero
 *        in the semi-empirical approximation.
 */
class ZeroLocal2c2eIntegrals {
 public:
  /**
   * @brief Returns whether the two two electron integral types have as result 0 because of the NDDO approximation.
   * @param t1 first two electron integral type.
   * @param t2 second two electron integral type.
   * @return true if the integral <t1|1/r|t2> is zero due to the NDDO approximation.
   */
  static bool isZero(GeneralTypes::twoElIntegral_t t1, GeneralTypes::twoElIntegral_t t2);

 private:
  constexpr static int numberTwoElectronIntegrals = 40;
  using Array = std::array<std::array<bool, numberTwoElectronIntegrals>, numberTwoElectronIntegrals>;

  static Array getArrayForZeroIntegrals();
  static void setElementToNonZero(Array& array, GeneralTypes::twoElIntegral_t t1, GeneralTypes::twoElIntegral_t t2);
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ZEROLOCAL2C2EINTEGRALS_H
