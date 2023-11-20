/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_STO_6G_PARAMETERS_H
#define INCLUDE_SPARROW_STO_6G_PARAMETERS_H

#include "Utils/DataStructures/AtomicGtos.h"
#include <unordered_map>

namespace Scine {
namespace Sparrow {
namespace Sto6g {

//! Default parameters for RM1
std::unordered_map<int, Utils::AtomicGtos> rm1();
//! Default parameters for AM1
std::unordered_map<int, Utils::AtomicGtos> am1();
//! Default parameters for PM3
std::unordered_map<int, Utils::AtomicGtos> pm3();
//! Default parameters for PM6
std::unordered_map<int, Utils::AtomicGtos> pm6();
//! Default parameters for MNDO
std::unordered_map<int, Utils::AtomicGtos> mndo();
//! Default parameters for DFTB
std::unordered_map<int, Utils::AtomicGtos> dftb();

} // namespace Sto6g
} // namespace Sparrow
} // namespace Scine

#endif
