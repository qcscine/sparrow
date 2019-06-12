/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MULTIPOLE_MULTIPOLETYPES_H
#define SPARROW_MULTIPOLE_MULTIPOLETYPES_H

#include "GeneralTypes.h"

namespace Scine {
namespace Sparrow {

/*!
 * \file generalTypes.h
 * This file defines multiple functions for the calculation of
 * semi-empirical two-center two-electron integrals
 */

namespace nddo {

namespace multipole {

/** @brief enum listing the possible charge configurations of a multipole.
 * It is possible i.e. to infer the charge separation D_{sp1} from it.
 * They are separated in monopole (l = 0), dipole (l = 1) and quadrupole (l = 2).
 */

enum multipolePair_t { sp1, pd1, pp2, sd2, dd2, ss0, pp0, dd0 };
/**
 * @brief Multipole types used in the calculation of the ERI with the multipole expansion approximation.
 *
 * Confusion might arise by the use of two merged formalisms: the one for just s and p orbitals
 * and the one for s, p and d orbitals.
 *
 * l: orbital quantum number
 * m: magnetic quantum number
 * \f$ M00 = M_{0,0} = q^I\f$ a monopole with l = 0, m = 0
 * \f$ M1m1 = M_{1,-1} = \mu_y \f$ a dipole in y direction with l = 1, m = -1
 * \f$ M10 = M_{1,0} = \mu_z \f$ a dipole in z direction with l = 1, m = 0
 * \f$ M11 = M_{1,1} = \mu_x \f$ a dipole in x direction with l = 1, m = 1
 * \f$ Qxx = Q_{x,x} \f$ a linear quadrupole in x direction with l = 2, m = 0
 * \f$ Qyy = Q_{y,y} \f$ a linear quadrupole in y direction with l = 2, m = 0
 * \f$ Qzz = Q_{z,z} \f$ a linear quadrupole in z direction with l = 2, m = 0
 * \f$ M2m2 = M_{2,-2} = Q_{x,y} \f$ a x,y square quadrupole with l = 2, m = -2
 * \f$ M2m1 = M_{2,-1} = Q_{y,z} \f$ a y,z square quadrupole with l = 2, m = -1
 * \f$ M20 = M_{2,0} = -\~{Q}_{x,z} - \frac{1}{2}\~{Q}_{x,y} \f$ a quadrupole with charges along each axis
 * at \f$ \sqrt{2} \f$ distance from the origin with l = 2, m = 0
 * \f$ M21 = M_{2,1} = Q_{x,z}\f$ a x,z square quadrupole with l = 2, m = 1
 * \f$ M22 = M_{2,2} = \~{Q}_{x,y} \f$ a square quadrupole with charges along the x,y axes at \f$ \sqrt{2} \f$
 * distance from the origin with l = 2, m = 2
 * \f$ Qzx = \~{Q}_{z,x} = -\~{Q}_{x,z} \f$ a square quadrupole with charges along the x,z axes at \f$ \sqrt{2} \f$
 * distance from the origin with l = 2, m = 2
 */
enum multipole_t { M00, Qxx, Qyy, Qzz, M1m1, M10, M11, M2m2, M2m1, M20, M21, M22, Qzx };
/**
 * @brief Given 2 orbitals, gives the corresponding orbital pair. Throws InvalidOrbitalPairException() if the orbital
 *        types given are invalid. Order matters.
 * @param o1 first orbital
 * @param o2 second orbital
 * @return a GeneralTypes::rotationOrbitalPair corresponding to the input orbitals
 */
GeneralTypes::rotationOrbitalPair getRotPairType(GeneralTypes::orb_t o1, GeneralTypes::orb_t o2);

/** @brief Returns the magnetic quantum number m of a multipole, i.e. -1 for a dipole in y direction,
 *         0 for a dipole in z direction and 1 for a dipole in x direction. Throws InvalidMultipoleException() if the
 *         multipole is not valid.
 */
int MQuantumNumber(multipole_t m);

/** @brief Returns the orbital quantum number l of an orbital, i.e. 0 for s, 1 for p and 2 for d type orbitals.
 *         Throws InvalidMultipoleException() if the multipole is not a valid one.
 */
int LQuantumNumber(multipole_t m);

/**
 * @brief Function to infer the charge configuration of a multipole
 * @param l1 the orbital quantum number of the first orbital
 * @param l2 the orbital quantum number of the second orbital
 * @param l  the multipole orbital quantum number
 * @return throws InvalidQuantumNumbersException() if the quantum number is invalid. Otherwise returns a multipolePair_t.
 */
multipolePair_t pairType(int l1, int l2, int l);

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MULTIPOLE_MULTIPOLETYPES_H
