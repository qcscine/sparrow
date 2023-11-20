/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SKPair.h"
#include "SKAtom.h"
#include <Utils/IO/Regex.h>
#include <Utils/Scf/MethodExceptions.h>
#include <Utils/Technical/ScopedLocale.h>
#include <cmath>
#include <fstream>
#include <regex>

namespace Scine {
namespace Sparrow {

using std::exp;
using namespace Utils::AutomaticDifferentiation;

namespace dftb {

// Declare space for values in binary at runtime
constexpr decltype(SKPair::integralIndexes) SKPair::integralIndexes;

SKPair::SKPair(SKAtom* atomicParameters1, SKAtom* atomicParameters2, SkfData data)
  : atomType1(atomicParameters1),
    atomType2(atomicParameters2),
    gridDist(data.gridDistance),
    nGridPoints(data.integralTable.front().size()),
    rMax(gridDist * nGridPoints + distFudge),
    integralTable(std::move(data.integralTable)),
    repulsion_(std::move(data.repulsion)),
    extrC3(28),
    extrC4(28),
    extrC5(28) {
  if (integralTable.back().empty()) {
    throw std::runtime_error("Back columns of integral table are empty instead of default-initialized!");
  }

  if (data.atomicParameters) {
    const auto& p = data.atomicParameters.value();
    atomicParameters1->setEnergies(p.Es, p.Ep, p.Ed);
    atomicParameters1->setHubbardParameter(p.Us, p.Up, p.Ud);
    atomicParameters1->setOccupations(p.fs, p.fp, p.fd);
  }

  // Set the number of integrals for the atom pair
  const int nAOsA = atomicParameters1->getnAOs();
  const int nAOsB = atomicParameters2->getnAOs();
  if (nAOsA + nAOsB == 2) { // Only s orbitals
    nIntegrals = 2;
  }
  else if (nAOsA + nAOsB == 5) { // one atom only has s, the other one s and p
    nIntegrals = 4;
  }
  else if (nAOsA + nAOsB == 8) { // Both atoms have s and p
    nIntegrals = 10;
  }
  else if (nAOsA + nAOsB == 10) { // one has s, one has d
    nIntegrals = 14;
  }
  else if (nAOsA + nAOsB == 13) { // one has p, one has d
    nIntegrals = 22;
  }
  else if (nAOsA + nAOsB == 18) { // both have d
    nIntegrals = 28;
  }
}

void SKPair::precalculateGammaTerms() {
  // Precondition: not same atom type. In this case the gamma terms are not needed
  if (atomType1 == atomType2)
    return;

  double ta = atomType1->getHubbardParameter() * 3.2;
  double tb = atomType2->getHubbardParameter() * 3.2;
  double ta2 = ta * ta;
  double tb2 = tb * tb;

  gamma.g1a = (tb2 * tb2 * ta) / (2.0 * (ta2 - tb2) * (ta2 - tb2));
  gamma.g1b = (ta2 * ta2 * tb) / (2.0 * (tb2 - ta2) * (tb2 - ta2));
  gamma.g2a = tb2 * tb2 * (tb2 - 3 * ta2) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  gamma.g2b = ta2 * ta2 * (ta2 - 3 * tb2) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));

  gammaDerivative.dgab1a = (tb2 * tb2 * tb2 + 3.0 * ta2 * tb2 * tb2) / (2.0 * (ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  gammaDerivative.dgab2a = (12.0 * ta2 * ta * tb2 * tb2) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  gammaDerivative.dgab1b = (2.0 * ta2 * ta * tb2 * tb) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  gammaDerivative.dgab2b = (12.0 * ta2 * ta2 * tb2 * tb) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  gammaDerivative.dgba1a = (ta2 * ta2 * ta2 + 3.0 * tb2 * ta2 * ta2) / (2.0 * (tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
  gammaDerivative.dgba2a = (12.0 * tb2 * tb * ta2 * ta2) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
  gammaDerivative.dgba1b = (2.0 * tb2 * tb * ta2 * ta) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
  gammaDerivative.dgba2b = (12.0 * tb2 * tb2 * ta2 * ta) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));

  gammaDerivative.dgadr = (tb2 * tb2 * tb2 - 3.0 * ta2 * tb2 * tb2) / ((ta2 - tb2) * (ta2 - tb2) * (ta2 - tb2));
  gammaDerivative.dgbdr = (ta2 * ta2 * ta2 - 3.0 * tb2 * ta2 * ta2) / ((tb2 - ta2) * (tb2 - ta2) * (tb2 - ta2));
}

const SKPair::GammaTerms& SKPair::getGammaTerms() const {
  return gamma;
}

const SKPair::GammaDerivativeTerms& SKPair::getGammaDerTerms() const {
  return gammaDerivative;
}

void SKPair::complete(SKPair* p) {
  assert(p != nullptr);
  constexpr double factor = -1.0;
  for (int i = 0; i < nGridPoints; i++) {
    integralTable[20][i] = factor * p->integralTable[3][i];
    integralTable[21][i] = factor * p->integralTable[4][i];
    integralTable[22][i] = p->integralTable[7][i];
    integralTable[23][i] = factor * p->integralTable[8][i];
    integralTable[24][i] = factor * p->integralTable[13][i];
    integralTable[25][i] = factor * p->integralTable[14][i];
    integralTable[26][i] = p->integralTable[17][i];
    integralTable[27][i] = factor * p->integralTable[18][i];
  }
  // Precompute the coefficients for the extrapolation
  precompute5Extrapolation();
}

void SKPair::precompute5Extrapolation() {
  int indStart = nGridPoints - nInter;
  InterpolationValues<Utils::DerivativeOrder::One> t1;
  InterpolationValues<Utils::DerivativeOrder::One> t2;

  // Calculate interpolated values right after and right before the last given
  // point in parameter table; these values are used to compute the numerical
  // derivatives at the last point
  interpolate(t1, static_cast<double>(nGridPoints - 1) - deltaR, indStart); // Get values right before last given point
  interpolate(t2, static_cast<double>(nGridPoints - 1) + deltaR, indStart); // Get values right after last given point

  double yd0, yd1, yd2; // zero, first and second order derivatives
  double dx1, dx2;      // Temporary variables needed for the extrapolation

  for (int L = 0; L < nIntegrals; L++) {
    assert(integralIndexes[L] < static_cast<int>(integralTable.size()));
    yd0 = integralTable[integralIndexes[L]][nGridPoints - 1];
    yd1 = t1.derivIntegral[L].derivative();
    yd2 = (t2.derivIntegral[L].derivative() - t1.derivIntegral[L].derivative()) / (deltaR);

    // 5th order extrapolation (impose derivatives at last gridpoint, and zero value and derivatives at rMax)
    dx1 = -yd1 * (distFudge / gridDist);
    dx2 = yd2 * (distFudge / gridDist) * (distFudge / gridDist);
    extrC3[L] = 10 * yd0 - 4 * dx1 + 0.5 * dx2;
    extrC4[L] = -15 * yd0 + 7 * dx1 - dx2;
    extrC5[L] = 6 * yd0 - 3 * dx1 + 0.5 * dx2;
  }
}

template<Utils::DerivativeOrder O>
Value1DType<O> SKPair::getRepulsion(double const& r) const {
  auto R = variableWithUnitDerivative<O>(r);
  if (r > repulsion_.cutoff)
    return constant1D<O>(0);
  if (r < repulsion_.splines[0].start)
    return exp(-repulsion_.a1 * R + repulsion_.a2) + repulsion_.a3;

  auto i = static_cast<int>((r - repulsion_.splines[0].start) / (repulsion_.cutoff - repulsion_.splines[0].start) *
                            repulsion_.nSplineInts);

  // If not the right bin, find the right one:
  if (repulsion_.splines[i].start > r) {
    while (repulsion_.splines[--i].start > r) {
    }
  }
  else if (repulsion_.splines[i].end < r) {
    while (repulsion_.splines[++i].end < r) {
    }
  }

  const auto& spline = repulsion_.splines[i];

  auto dr = R - spline.start;
  auto conditional = (i == repulsion_.nSplineInts - 1 ? dr * (repulsion_.c4 + dr * repulsion_.c5) : constant1D<O>(0.0));

  // clang-format off
  return (
    spline.c0 + dr * (
      spline.c1 + dr * (
        spline.c2 + dr * (
          spline.c3 + conditional
        )
      )
    )
  );
  // clang-format on
}

template<Utils::DerivativeOrder O>
int SKPair::getHS(double dist, InterpolationValues<O>& val) const {
  // If 'dist' too short or too large, return 0
  if (dist < gridDist || dist > rMax) {
    for (int L = 0; L < nIntegrals; L++)
      val.derivIntegral[L] = constant1D<O>(0.0);
    return 0;
  }

  // If 'dist' is in the extrapolation zone, take simple formula
  if (dist > nGridPoints * gridDist) {
    double dr = dist - rMax;
    auto xr = variableWithUnitDerivative<O>(-dr / distFudge);

    // 5th order extrapolation (impose derivatives at last gridpoint, and zero value and derivatives at rMax)
    for (int L = 0; L < nIntegrals; L++)
      val.derivIntegral[L] = ((extrC5[L] * xr + extrC4[L]) * xr + extrC3[L]) * xr * xr * xr;

    return 1;
  }

  // Else: do interpolation
  getHSIntegral(val, dist);

  return 1;
}

template<Utils::DerivativeOrder O>
int SKPair::getHSIntegral(InterpolationValues<O>& val, double dist) const {
  double position = dist / gridDist - 1.0;
  int ind; // Current index during interpolation
  ind = static_cast<int>(position);

  int indStart; // first index for interpolation

  indStart = (ind >= nGridPoints - nInterRight) ? nGridPoints - nInter : ((ind < nInterLeft) ? 0 : ind - nInterLeft + 1);
  // -1 because loop afterward goes from i = 0 to nInter - 1 ->
  // first item accessed is indStart + 0, last element
  // accessed is indStart + nInter - 1
  assert(indStart + nInter - 1 < nGridPoints);
  interpolate(val, position, indStart);

  return 1;
}

template<Utils::DerivativeOrder O>
void SKPair::interpolate(InterpolationValues<O>& val, double x, int start) const {
  /*
   * Taken from Numerical recipes
   * Adapted to treat simultaneously the different orbitals
   */
  double dift, dif = x - start;
  int ns = 0; // index of closest table entry
  for (int L = 0; L < nIntegrals; L++) {
    for (int i = 0; i < nInter; i++) {
      val.ya[L][i] = integralTable[integralIndexes[L]][start + i];
      val.C[L][i] = constant1D<O>(integralTable[integralIndexes[L]][start + i]);
      val.D[L][i] = constant1D<O>(integralTable[integralIndexes[L]][start + i]);
    }
  }
  for (int i = 0; i < nInter; i++) {
    val.xa[i] = static_cast<double>(start + i);
    if ((dift = std::abs(x - val.xa[i])) < dif) {
      dif = dift;
      ns = i;
    }
  }
  for (int L = 0; L < nIntegrals; L++) {
    val.derivIntegral[L] = constant1D<O>(val.ya[L][ns]);
  }
  ns--;

  Value1DType<O> w, den;
  double ho, hp;
  for (int m = 1; m < nInter; m++) {
    for (int i = 0; i <= nInter - 1 - m; i++) {
      ho = val.xa[i] - x;
      hp = val.xa[i + m] - x;
      for (int L = 0; L < nIntegrals; L++) {
        w = val.C[L][i + 1] - val.D[L][i];
        den = w / (ho - hp);
        val.D[L][i] = getFromFull<O>(hp, -1.0 / gridDist, 0) * den;
        val.C[L][i] = getFromFull<O>(ho, -1.0 / gridDist, 0) * den;
      }
    }

    for (int L = 0; L < nIntegrals; L++) {
      val.derivIntegral[L] += (2 * (ns + 1) < (nInter - m)) ? val.C[L][ns + 1] : val.D[L][ns];
    }
    if (!(2 * (ns + 1) < (nInter - m)))
      ns--;
  }
}

template int SKPair::getHS<Utils::DerivativeOrder::Zero>(double, InterpolationValues<Utils::DerivativeOrder::Zero>&) const;
template int SKPair::getHS<Utils::DerivativeOrder::One>(double, InterpolationValues<Utils::DerivativeOrder::One>&) const;
template int SKPair::getHS<Utils::DerivativeOrder::Two>(double, InterpolationValues<Utils::DerivativeOrder::Two>&) const;
template int SKPair::getHSIntegral<Utils::DerivativeOrder::Zero>(InterpolationValues<Utils::DerivativeOrder::Zero>&, double) const;
template int SKPair::getHSIntegral<Utils::DerivativeOrder::One>(InterpolationValues<Utils::DerivativeOrder::One>&, double) const;
template int SKPair::getHSIntegral<Utils::DerivativeOrder::Two>(InterpolationValues<Utils::DerivativeOrder::Two>&, double) const;
template void SKPair::interpolate<Utils::DerivativeOrder::Zero>(InterpolationValues<Utils::DerivativeOrder::Zero>&,
                                                                double, int) const;
template void SKPair::interpolate<Utils::DerivativeOrder::One>(InterpolationValues<Utils::DerivativeOrder::One>&,
                                                               double, int) const;
template void SKPair::interpolate<Utils::DerivativeOrder::Two>(InterpolationValues<Utils::DerivativeOrder::Two>&,
                                                               double, int) const;
template double SKPair::getRepulsion<Utils::DerivativeOrder::Zero>(double const& r) const;
template First1D SKPair::getRepulsion<Utils::DerivativeOrder::One>(double const& r) const;
template Second1D SKPair::getRepulsion<Utils::DerivativeOrder::Two>(double const& r) const;

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
