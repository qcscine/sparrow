/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_PM6PRINCIPALQUANTUMNUMBERS_H
#define SPARROW_PM6PRINCIPALQUANTUMNUMBERS_H

#include <Utils/Geometry/ElementTypes.h>

namespace Scine {
namespace Sparrow {

/*! \file PrincipalQuantumNumbers.h
 * Functions for NDDO-specific values concerning the elements.
 */

namespace nddo {

enum class BasisFunctions { sp, spd };

namespace PM6Elements {

unsigned int getQuantumNumberForSOrbital(Utils::ElementType e);

unsigned int getQuantumNumberForPOrbital(Utils::ElementType e);

unsigned int getQuantumNumberForDOrbital(Utils::ElementType e);

unsigned int getNumberOfAOs(Utils::ElementType e, BasisFunctions basisFunctions);

unsigned int getNumberOneCenterTwoElectronIntegrals(Utils::ElementType e, BasisFunctions basisFunctions);

double getCoreCharge(Utils::ElementType elementType);

} // namespace PM6Elements

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_PM6PRINCIPALQUANTUMNUMBERS_H
