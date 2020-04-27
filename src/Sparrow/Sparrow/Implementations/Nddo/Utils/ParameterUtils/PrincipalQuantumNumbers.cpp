/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "PrincipalQuantumNumbers.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Geometry/ElementTypes.h>
#include <cassert>

/* NOTE
 * - All of the calls to ElementInfo::Z in the right hand of comparisons are
 *   optimized away (the function is constexpr-evaluable), no runtime cost.
 */

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace PM6Elements {

unsigned int getQuantumNumberForSOrbital(Utils::ElementType e) {
  const unsigned Z = Utils::ElementInfo::Z(e);

  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::He))
    return 1;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::F))
    return 2;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Cl))
    return 3;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Br))
    return 4;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::I))
    return 5;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Bi))
    return 6;

  return 7;
}

unsigned int getQuantumNumberForPOrbital(Utils::ElementType e) {
  const unsigned Z = Utils::ElementInfo::Z(e);

  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ne))
    return 2;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ar))
    return 3;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Kr))
    return 4;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Xe))
    return 5;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Bi))
    return 6;

  return 7;
}

unsigned int getQuantumNumberForDOrbital(Utils::ElementType e) {
  const unsigned Z = Utils::ElementInfo::Z(e);

  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ge))
    return 3;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Sn))
    return 4;
  if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Bi))
    return 5;

  return 6;
}

unsigned int getNumberOfAOs(Utils::ElementType e, BasisFunctions basisFunctions) {
  unsigned int n = 0;
  const unsigned Z = Utils::ElementInfo::Z(e);

  if (e == Utils::ElementType::none)
    n = 0;
  else if (Z == Utils::ElementInfo::Z(Utils::ElementType::H))
    n = 1;
  else if (basisFunctions == BasisFunctions::sp)
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Mg))
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Cl))
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ca))
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Cu))
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ge))
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::As)) // CHANGED LATER
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Se)) // CHANGED LATER
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Br))
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Sr))
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ag))
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Sn))
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Sb)) // CHANGED LATER
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Te)) // CHANGED LATER
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::I))
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Ba))
    n = 4;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Au))
    n = 9;
  else if (Z <= Utils::ElementInfo::Z(Utils::ElementType::Bi))
    n = 4;

  return n;
}

unsigned int getNumberOneCenterTwoElectronIntegrals(Utils::ElementType e, BasisFunctions basisFunctions) {
  unsigned int n = getNumberOfAOs(e, basisFunctions);

  if (n == 1)
    return 1;
  if (n == 4)
    return 6;
  if (n == 9)
    return 58;

  return 0;
}

double getCoreCharge(Utils::ElementType elementType) {
  const unsigned Z = Utils::ElementInfo::Z(elementType);
  unsigned coreElectrons = 0;

  // Rare gases
  if (elementType == Utils::ElementType::Ne || elementType == Utils::ElementType::Ar ||
      elementType == Utils::ElementType::Kr || elementType == Utils::ElementType::Xe)
    return 6;

  if (Z > 86)
    coreElectrons = 86;
  else if (Z > 79)
    coreElectrons = 78;
  else if (Z > 70)
    coreElectrons = 68;
  else if (Z > 54)
    coreElectrons = 54;
  else if (Z > 47)
    coreElectrons = 46;
  else if (Z > 36)
    coreElectrons = 36;
  else if (Z > 29)
    coreElectrons = 28;
  else if (Z > 18)
    coreElectrons = 18;
  else if (Z > 10)
    coreElectrons = 10;
  else if (Z > 2)
    coreElectrons = 2;

  assert(Z >= coreElectrons && "Unsigned underflow!");
  return Z - coreElectrons;
}

} // namespace PM6Elements

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
