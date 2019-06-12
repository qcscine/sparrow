/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_KLOPMANPARAMETER_H
#define SPARROW_KLOPMANPARAMETER_H

#include "Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/multipoleTypes.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/*!
 * @brief This class is the container for the Klopman-Ohno parameters
 *        used for the evaluation of the multipoles
 *
 * The Klopman-Ohno are used in the calculation of the two-center ERIs in the NDDO formalism.
 * \f$ U\left(\Theta^{\mu\nu}_t, \Theta^{\lambda\sigma}_s \right) = \sum^{C_t}_{c=1}\sum^{C_s}_{d=1}
 * \frac{q_c^Iq_d^J}{\sqrt{|r_c^I - r_d^J|^2 +
 * \left( \theta_c^I(\chi_\mu^I\chi_\nu^I) + \theta_d^J(\chi_\lambda^J\chi_\sigma^J) \right)}} \f$
 */

class KlopmanParameter {
 public:
  /**
   * @brief Constructor, calls the reset() function.
   */
  KlopmanParameter();

  /** @brief Sets the values to zero. */
  void reset();

  /** @brief stores the value of the Klopman--Ohno parameter for a multipole type
   *         at the index given by the cast multipolePair_t. */
  void set(multipolePair_t type, double value) {
    K_[type] = value;
  }

  /** @brief returns the value correspoding to the multipolePair_t typed multipole */
  double get(multipolePair_t type) const {
    return K_[type];
  }

  /**
   * @brief generates the Klopman--Ohno parameters up to the s orbitals
   * @param gss The one center repulsion parameter \f$\gamma_{ss}\f$ parametrizing the \f$ \langle ss|ss \rangle\f$
   *        integral
   * @return The Klopman--Ohno parameter for the s orbital
   */
  void generateUpToS(double gss);
  /**
   * @brief  generates the Klopman--Ohno parameters up to the p orbitals
   * @param gss the one center repulsion parameter \f$\gamma_{ss}\f$ parametrizing the \f$ \langle ss|ss \rangle\f$
   *        integral
   * @param hsp the one center repulsion parameter \f$\~{\gamma}_{sp}\f$ parametrizing the \f$ \langle sp|sp \rangle\f$
   *        integral
   * @param D1sp the charge separation parameter describing the charge distance on the dipoles formed by an s and a p
   *        orbitals
   * @param hpp the one center repulsion parameter \f$ \gamma_{pp} \f$ parametrizing the \f$ \langle pp|pp \rangle \f$
   *        integral
   * @param D2pp the charge separation parameter describing the charge distance on the quadrupoles formed by two p
   *        orbitals
   *
   * The parameter is calculated in an iterative fashion with the Newton's method. The equation of which the root
   * is searched are listed in
   * T. Husch, A. C. Vaucher, M. Reiher, Semiempirical Molecular Orbital Models based on the Neglect of Diatomic
   * Differential Overlap Approximation, doi: 1806.06147v2
   */
  void generateUpToP(double gss, double hsp, double D1sp, double hpp, double D2pp);
  /**
   * @brief  generates the Klopman--Ohno parameters up to the d orbitals
   * @param gss the one center repulsion parameter \f$\gamma_{ss}\f$ parametrizing the \f$ \langle ss|ss \rangle\f$
   *        integral
   * @param hsp the one center repulsion parameter \f$\~{\gamma}_{sp}\f$ parametrizing the \f$ \langle sp|sp \rangle\f$
   *        integral
   * @param D1sp the charge separation parameter describing the charge distance on the dipoles formed by an s and a p
   *        orbitals
   * @param hpp the one center repulsion parameter \f$ \gamma_{pp} \f$ parametrizing the \f$ \langle pp|pp \rangle \f$
   *        integral
   * @param D2pp the charge separation parameter describing the charge distance on the quadrupoles formed by two p
   *        orbitals
   * @param F0dd the parametrized one-center radial integral corresponding to the zero-th term of equation 116 in the
   *        reference listed below between two d orbitals in the order given by equation 118
   *        in the reference listed below.
   * @param G1pd the parametrized one-center radial integral corresponding to the first term of equation 116 in the
   *        reference listed below between an s and a p orbital in the order given by equation 119.
   * @param D1pd the charge separation parameter describing the charge distance on the quadrupoles formed by a p and
   *        a d orbital
   * @param G2sd the parametrized one-center radial integral corresponding to the second term of equation 116 in the
   *        reference listed below between an s and a d orbital in the order given by equation 119
   * @param D2sd the charge separation parameter describing the charge distance on the quadrupoles formed by an s and
   *        a d orbital
   * @param F2dd the parametrized one-center radial integral corresponding to the second term of quation 116 in the
   *        reference listed below between two d orbitals in the order given by equation 118
   * @param D2dd the charge separation parameter describing the charge distance on the quadrupoles formed by two d
   *        orbitals
   *
   * The parameter is calculated in an iterative fashion with the Newton's method. The equation of which the root
   * is searched are listed in
   * T. Husch, A. C. Vaucher, M. Reiher, Semiempirical Molecular Orbital Models based on the Neglect of Diatomic
   * Differential Overlap Approximation, doi: 1806.06147v2
   */
  void generateUpToD(double gss, double hsp, double D1sp, double hpp, double D2pp, double F0dd, double G1pd,
                     double D1pd, double G2sd, double D2sd, double F2dd, double D2dd);

 private:
  double findRootVar1(double D, double shift);
  double findRootVar2(double D, double shift);
  /*! This formula is to be used when the D's given out by mopac are used. Some of them are too large by a factor
   * sqrt(2). */
  double findRootVar2Wrong(double D, double shift);

  double K_[8];
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_KLOPMANPARAMETER_H
