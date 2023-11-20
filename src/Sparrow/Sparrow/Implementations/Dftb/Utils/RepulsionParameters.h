/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
  struct Spline {
    double start;
    double end;
    double c0;
    double c1;
    double c2;
    double c3;
  };

  int nSplineInts;
  double cutoff;
  double a1, a2, a3; // Coefficients for exponential part of repulsion
  std::vector<Spline> splines;
  double c4, c5; // For last spline
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_REPULSIONPARAMETERS_H
