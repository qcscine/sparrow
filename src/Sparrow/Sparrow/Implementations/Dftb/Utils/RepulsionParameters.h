/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_REPULSIONPARAMETERS_H
#define SPARROW_DFTB_REPULSIONPARAMETERS_H

#include <vector>

namespace Scine {
namespace Sparrow {

namespace dftb {

/*!
 * DFTB parameters for the repulsion of an element pair.
 */
struct RepulsionParameters {
  int nSplineInts;
  double cutoff;
  double a1, a2, a3; // Coefficients for exponential part of repulsion
  std::vector<double> splineStart, splineEnd;
  std::vector<double> c0, c1, c2, c3;
  double c4, c5; // For last spline
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_REPULSIONPARAMETERS_H