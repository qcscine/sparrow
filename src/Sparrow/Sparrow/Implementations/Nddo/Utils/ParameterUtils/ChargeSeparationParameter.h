/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_CHARGESEPARATIONPARAMETER_H
#define SPARROW_CHARGESEPARATIONPARAMETER_H

#include "Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/multipoleTypes.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/*!
 * @brief   Charge separation D of semi-empirical models. It describes the separation
 *          between two charges of opposite sign in a multipole.
 *
 * The calculation of the charge separation is described in Thiel, Voityuk, Theor Chim Acta, 1992, 81, 391.
 * The charge separation are stored ad a static c-array of size 5. The charge separations are ordered as in the
 * nddo::multipole::multipolePair_t enum, i.e. sp1, pd1, pp2, sd2, dd2.
 * ss0, pp0, dd0 are not present as there is no charge separation for them (they are monopoles).
 */

class ChargeSeparationParameter {
 public:
  /** @brief constructor, resets the charge separation array to zero. */
  ChargeSeparationParameter();

  /**
   * @brief calculates the charge separation from the main quantum numbers and the orbital exponents of s and p orbitals
   * @param ns main quantum number of the s orbital
   * @param np main quantum number of the p orbital
   * @param zs orbital exponent of the s orbital
   * @param zp orbital exponent of the p orbital
   */
  void computeFromExponents(unsigned int ns, unsigned int np, double zs, double zp);
  /**
   * @brief calculate the charge separation from the main quantum numbers and the orbital exponents of s, p and d
   *        orbitals
   * @param ns main quantum number of the s orbital
   * @param np main quantum number of the p orbital
   * @param nd main quantum number of the d orbital
   * @param zs orbital exponent of the s orbital
   * @param zp orbital exponent of the p orbital
   * @param zd orbital exponent of the d orbital
   */
  void computeFromExponents(unsigned int ns, unsigned int np, unsigned int nd, double zs, double zp, double zd);
  /**
   * @brief sets the value of a charge separation for a given multipole type
   * @param type one of the first 5 elements of the nddo::multipole:multipolePair_t enum, corresponding to the
   *        dipole and the quadrupole charge configurations. The assignment is not checked.
   * @param value the value of the charge separation parameter
   */
  void set(multipolePair_t type, double value);
  /**
   * @brief gets the value corresponding to one of the multipolescharge configuration.
   *
   * No check is performed to ensure that no monopolar charge configuration is requested.
   */
  double get(multipolePair_t type) const;

 private:
  //! Calculates the factorial recursively
  long long factorial(long long n) const;
  //! Calculates the A factor. It is necessary to calculate the charge separations.
  double A(unsigned int n1, unsigned int n2, double z1, double z2, int L) const;
  double D_[5];
};

inline void ChargeSeparationParameter::set(nddo::multipole::multipolePair_t type, double value) {
  D_[type] = value;
}

inline double ChargeSeparationParameter::get(nddo::multipole::multipolePair_t type) const {
  return D_[type];
}

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_CHARGESEPARATIONPARAMETER_H
