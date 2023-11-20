/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ORDERTAG_H
#define SPARROW_ORDERTAG_H

namespace Scine {
namespace Sparrow {
/**
 * @brief Tag to define the order in the preconditioner evaluators.
 * This is in a single file to prevent needing to include a huge header for
 * this.
 */
struct OrderTag {};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_ORDERTAG_H
