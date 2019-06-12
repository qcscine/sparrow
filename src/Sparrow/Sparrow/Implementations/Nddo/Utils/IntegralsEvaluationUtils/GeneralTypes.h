/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_GENERALTYPES_H
#define SPARROW_GENERALTYPES_H

#include <exception>
#include <utility>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace GeneralTypes {

class InvalidMultipoleException : public std::exception {};
class InvalidOrbitalPairException : public std::exception {};
class InvalidQuantumNumbersException : public std::exception {};

/**
 *  @brief enum indicating the possible orbitals
 *
 *  orb_t::s A s orbital
 *  orb_t::x A p_x orbital
 *  orb_t::y A p_y orbital
 *  orb_t::z A p_z orbital
 *  orb_t::x2y2 A d_{x^2-y^2} orbital
 *  orb_t::xz A d_{xz} orbital
 *  orb_t::z2 A d_{z^2} orbital
 *  orb_t::yz A d_{yz} orbital
 *  orb_t::xy A d_{xy} orbital
 */
enum orb_t { s, x, y, z, x2y2, xz, z2, yz, xy };

// clang-format off
/**
 * @brief enum listing all of the orbital pairs giving rise to a valid charge distribution
 */
enum class twoElIntegral_t {//1 ss, 9 sp/pp, 30 sd/pd/dd
                      s_s,
                      s_x,x_x,s_y,x_y,y_y,s_z,x_z,y_z,z_z,
                      s_z2,s_xz,s_yz,s_x2y2,s_xy,
                      x_z2,x_xz,x_x2y2,x_xy, //x_yz would always be zero 
                      y_z2,y_yz,y_x2y2,y_xy, //y_xz would always be zero
                      z_z2,z_xz,z_yz,//z_x2y2 and z_xy would always be zero   
                      z2_z2,z2_xz,z2_yz,z2_x2y2,z2_xy,
                      xz_xz,xz_yz,xz_x2y2,xz_xy,
                      yz_yz,yz_x2y2,yz_xy,
                      x2y2_x2y2,xy_xy //x2y2_xy would always be zero
                    };
enum class rotationOrbitalPair {s_s, 
                          x_x, x_y, x_z, y_x, y_y, y_z, z_x, z_y, z_z,
                          x2y2_x2y2, x2y2_xz, x2y2_z2, x2y2_yz, x2y2_xy, 
                          xz_x2y2, xz_xz, xz_z2, xz_yz, xz_xy, 
                          z2_x2y2, z2_xz, z2_z2, z2_yz, z2_xy, 
                          yz_x2y2, yz_xz, yz_z2, yz_yz, yz_xy, 
                          xy_x2y2, xy_xz, xy_z2, xy_yz, xy_xy};
// clang-format on
/**
 * @brief gets the orbital quantum number ("l") of the argument, i.e. 0 for s, 1 for px,py,pz, 2 for all d orbitals.
 * @return Returns 0 if the orbital type is invalid.
 */
int orbitalQN(orb_t o);

/**
 * @brief separates an orbital pair into its orbital components, throws InvalidOrbitalPairException if pair not valid
 * @param an orbital pair, i.e. s_x, corresponding to the s and the p_x orbitals
 * @return a std::pair with orbital type elements, i.e. for s_x a pair with an orb_t::s and a orb_t::x
 */
std::pair<orb_t, orb_t> separatePair(twoElIntegral_t t);

} // namespace GeneralTypes

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_GENERALTYPES_H
