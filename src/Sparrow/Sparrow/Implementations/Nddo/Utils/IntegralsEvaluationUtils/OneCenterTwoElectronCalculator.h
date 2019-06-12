/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ONECENTERTWOELECTRONCALCULATOR_H
#define SPARROW_ONECENTERTWOELECTRONCALCULATOR_H

#include "OneCenterTwoElectronIntegralExpression.h"

namespace Scine {
namespace Sparrow {

namespace nddo {
class SlaterCondonParameters;

/*!
 * This class implements formulas for the calculation of the
 * one-center two-electron integrals from 14 Slater-Condon parameters
 * and three radial parameters. The main reference is
 * Pelikan, P; Turi Nagy L., Chemical Papers, 1974, 28, 594-598.
 * The indexes used here are the same as in the reference minus one.
 * In formula 17 of the reference, sqrt(12) is used instead of 12.
 * In formula 54, R1sppd would be correct instead of R1spdd
 * In formulas 51, 53, 56, 57, R2sppd is R2sdpp
 */
class OneCenterTwoElectronCalculator {
 public:
  //! @brief constructor
  OneCenterTwoElectronCalculator();

  /**
   * @brief Generates the set of used indexes for the one center two electrons statically
   * This functions defines the first time the indexes through the private method setUniqueIndexes().
   * The generated index array (containing all the possible orbital combination, so a 9x9x9x9 array) contains
   * for each element the index of the corresponding integral value.
   */
  static void setIndexes();
  static int getIndex(int i, int j, int k, int l);
  double getIntegral(int ind, const SlaterCondonParameters* p);

 private:
  static void setUniqueIndexes();
  static void setIndex(int i, int j, int k, int l, int* p);
  // Is the index array already set?
  static bool indexesSet_;
  // The index array containing for each orbital combination an index to the corresponding integral value
  static int index[9][9][9][9];
  // The expression of the integral values, accessed through the indexes in index[9][9][9][9]
  static OneCenterTwoElectronIntegralExpression expr[58];
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ONECENTERTWOELECTRONCALCULATOR_H
