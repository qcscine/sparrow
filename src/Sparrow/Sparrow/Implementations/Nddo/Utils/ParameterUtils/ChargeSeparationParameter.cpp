/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "ChargeSeparationParameter.h"
#include <cmath>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

ChargeSeparationParameter::ChargeSeparationParameter() {
  for (double& i : D_) {
    i = 0;
  }
}

void ChargeSeparationParameter::computeFromExponents(unsigned int ns, unsigned int np, double zs, double zp) {
  // set(sp1, A(ns,np,zs,zp,1) / std::sqrt(3)); THIS SHOULD BE THE CORRECT FORMULA!
  unsigned int n = (ns <= np) ? ns : np;
  set(sp1, A(n, n, zs, zp, 1) / std::sqrt(3));
  // set(pp2, std::sqrt(A(np,np,zp,zp,2)) / std::sqrt(5));THIS SHOULD BE THE CORRECT FORMULA!
  set(pp2, std::sqrt(A(n, n, zp, zp, 2)) / std::sqrt(5));
}

void ChargeSeparationParameter::computeFromExponents(unsigned int ns, unsigned int np, unsigned int nd, double zs,
                                                     double zp, double zd) {
  computeFromExponents(ns, np, zs, zp);
  set(pd1, A(np, nd, zp, zd, 1) / std::sqrt(5));
  set(sd2, std::sqrt(A(ns, nd, zs, zd, 2)) / std::pow(15.0, 0.25));
  set(dd2, std::sqrt(A(nd, nd, zd, zd, 2)) / std::sqrt(7));
}

double ChargeSeparationParameter::A(unsigned int n1, unsigned int n2, double z1, double z2, int L) const {
  double f1 = std::pow(2 * z1, n1 + 0.5);
  double f2 = std::pow(2 * z2, n2 + 0.5);
  double f3 = std::pow(z1 + z2, (-1.0) * (n1 + n2 + L + 1));
  double f4 = 1.0 / std::sqrt(factorial(2 * n1) * factorial(2 * n2));
  auto f5 = static_cast<double>(factorial(n1 + n2 + L));
  return f1 * f2 * f3 * f4 * f5;
}

long long ChargeSeparationParameter::factorial(long long n) const {
  return (n == 0 || n == 1) ? 1 : n * factorial(n - 1);
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
