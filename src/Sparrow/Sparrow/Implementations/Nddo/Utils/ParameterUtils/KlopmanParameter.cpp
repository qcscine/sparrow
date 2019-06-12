/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "KlopmanParameter.h"
#include <cmath>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

KlopmanParameter::KlopmanParameter() {
  reset();
}

void KlopmanParameter::reset() {
  for (double& i : K_)
    i = 0;
}

void KlopmanParameter::generateUpToS(double gss) {
  reset();
  double v = 1.0 / (2.0 * gss);
  set(multipolePair_t::ss0, v);
}

void KlopmanParameter::generateUpToP(double gss, double hsp, double D1sp, double hpp, double D2pp) {
  generateUpToS(gss);
  set(multipolePair_t::pp0, get(multipolePair_t::ss0));
  // In Husch, Vaucger, Reiher, 2018 the factor is 4/3 * hsp, but it is an uncorrected mistake.
  set(sp1, findRootVar1(D1sp, 4 * hsp));
  set(pp2, findRootVar2(D2pp, 8 * hpp));
}

void KlopmanParameter::generateUpToD(double gss, double hsp, double D1sp, double hpp, double D2pp, double F0dd,
                                     double G1pd, double D1pd, double G2sd, double D2sd, double F2dd, double D2dd) {
  generateUpToP(gss, hsp, D1sp, hpp, D2pp);
  set(dd0, 1.0 / (2.0 * F0dd));
  set(pd1, findRootVar1(D1pd, 16.0 / 15.0 * G1pd));
  set(sd2, findRootVar2(D2sd, 8.0 / 5.0 * G2sd));
  set(dd2, findRootVar2(D2dd, 24.0 / 49.0 * F2dd));
}

double KlopmanParameter::findRootVar1(double D, double shift) {
  // Newton's method
  double guess = 0.001;
  double previous = 0;

  while (std::abs(guess - previous) > 1e-8) {
    double r = guess;
    double r2 = r * r;
    double d2 = D * D;
    double sq = std::sqrt(r2 + d2);
    // Formula
    double f = 1.0 / r - 1.0 / sq - shift;
    // Derivative of the formula
    double der = -1.0 / r2 + r / (sq * sq * sq);

    previous = guess;
    guess = r - f / der;
  }

  return guess;
}

double KlopmanParameter::findRootVar2(double D, double shift) {
  // Newton's method
  double guess = 0.001;
  double previous = 0;

  while (std::abs(guess - previous) > 1e-8) {
    double r = guess;
    double r2 = r * r;
    double d2 = D * D;
    double sq1 = std::sqrt(r2 + d2);
    double sq2 = std::sqrt(r2 + 2 * d2);
    // Formula
    double f = 1.0 / r - 2.0 / sq1 + 1.0 / sq2 - shift;
    // Derivative of the formula
    double der = -1.0 / r2 + 2 * r / (sq1 * sq1 * sq1) - r / (sq2 * sq2 * sq2);

    previous = guess;
    guess = r - f / der;
  }

  return guess;
}

double KlopmanParameter::findRootVar2Wrong(double D, double shift) {
  // Newton's method
  double guess = 0.001;
  double previous = 0;

  while (std::abs(guess - previous) > 1e-8) {
    double r = guess;
    double r2 = r * r;
    double d2 = D * D;
    double sq1 = std::sqrt(r2 + d2 / 2);
    double sq2 = std::sqrt(r2 + d2);
    // Formula
    double f = 1.0 / r - 2.0 / sq1 + 1.0 / sq2 - shift;
    // Derivative of the formula
    double der = -1.0 / r2 + 2 * r / (sq1 * sq1 * sq1) - r / (sq2 * sq2 * sq2);

    previous = guess;
    guess = r - f / der;
  }

  return guess;
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
