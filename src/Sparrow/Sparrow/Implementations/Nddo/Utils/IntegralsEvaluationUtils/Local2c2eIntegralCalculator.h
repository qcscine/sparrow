/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_LOCAL2C2EINTEGRALCALCULATOR_H
#define SPARROW_LOCAL2C2EINTEGRALCALCULATOR_H

#include "MultipoleMultipoleInteractionContainer.h"
#include "multipoleTypes.h"
#include <array>
#include <list>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {
class ChargeSeparationParameter;
class KlopmanParameter;

/**
 * @brief This class is responsible for the calculation of the 2-center-2-electron integrals in the local coordinate
 *        system.
 *
 */

class Local2c2eIntegralCalculator {
 public:
  /**
   * @brief Struct defining an interaction between two multipoles.
   * It consists of a prefactor, f, the charge distributions and the corresponding multipoles.
   * It thus uniquely defines an interaction.
   */
  struct LocalTerm {
    double f;
    multipolePair_t p1, p2;
    multipole_t m1, m2;
  };

  //! Each charge distribution might be expanded in more than just one multipole (so, more than one LocalTerm possible).
  using LocalTerms = std::list<LocalTerm>;
  //! There are 40 possible charge distributions.
  using LocalTermArray = std::array<std::array<LocalTerms, 40>, 40>;

  //! Wrapper around the next overload. It enables the reference of the twoElIntegral_t without the need to cast them.
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> static getIntegral(
      GeneralTypes::twoElIntegral_t t1, GeneralTypes::twoElIntegral_t t2, double R, const ChargeSeparationParameter& D1,
      const ChargeSeparationParameter& d2, const KlopmanParameter& rho1, const KlopmanParameter& rho2);

  /**
   * @brief Calculates the integral value (up to the derivative order O) of the repulsion integral.
   * @tparam O the derivative order, defines up to which derivative the value must be computed.
   * @param t1 the index referring to the first GeneralTypes::twoElIntegral_t orbital pair.
   * @param t2 the index referring to the second GeneralTypes::twoElIntegral_t orbital pair.
   * @param R the distance between the two centers at which the multipoles are located.
   * @param D1 the charge separation parameter belonging to the first multipole.
   * @param d2 the charge separation parameter belonging to the second multipole.
   * @param rho1 the Klopman parameter for the first center.
   * @param rho2 the Klopman parameter for the second center.
   * @return The 2-center-2-electron value up to the O-th derivative in the local coordinates system.
   */
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> static getIntegral(int t1, int t2, double R,
                                                                     const ChargeSeparationParameter& D1,
                                                                     const ChargeSeparationParameter& d2,
                                                                     const KlopmanParameter& rho1,
                                                                     const KlopmanParameter& rho2);

 private:
  /// Generates the possible terms that define the multipole/multipole interactions.
  static LocalTermArray setUpTerms();
  /// Generates a multipole corresponding to the two electron integral type t.
  /// hasD makes sure that charge configuration corresponding to different multipoles in presence or absence
  /// of d basis functions are correctly represented
  static std::list<std::pair<double, multipole_t>> getMultipoles(GeneralTypes::twoElIntegral_t t, bool hasD);
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_LOCAL2C2EINTEGRALCALCULATOR_H
