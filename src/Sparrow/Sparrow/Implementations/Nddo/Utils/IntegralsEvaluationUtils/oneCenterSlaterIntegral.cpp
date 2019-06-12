/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "oneCenterSlaterIntegral.h"
#include <cassert>
#include <cmath>

namespace Scine {
namespace Sparrow {

namespace nddo {

void OneCenterSlaterIntegral::setPrincipal(int na, int nb, int nc, int nd) {
  na_ = na;
  nb_ = nb;
  nc_ = nc;
  nd_ = nd;
}

void OneCenterSlaterIntegral::setAngular(int la, int lb, int lc, int ld) {
  la_ = la;
  lb_ = lb;
  lc_ = lc;
  ld_ = ld;
}

void OneCenterSlaterIntegral::setExponents(double ea, double eb, double ec, double ed) {
  ea_ = ea;
  eb_ = eb;
  ec_ = ec;
  ed_ = ed;
}

double OneCenterSlaterIntegral::calculate(int l) {
  l_ = l;
  if (!validLValue())
    return 0.0;

  computeHelperVariables();

  double c = calculateFactor();
  double s1 = calculateFirstTerm();
  double s2 = calculateSecondTerm();
  double s3 = calculateThirdTerm();
  return c * (s1 - s2 + s3);
}

bool OneCenterSlaterIntegral::validLValue() const {
  /*
  if(l_ < (la_>lc_ ? la_-lc_ : lc_-la_) || l_ > la_+lc_)
    return false;
  else
    return true;
    */
  // NB: in kumar1987, they say it should be zero for the above condition,
  // but this would make F0pd zero, which it is not... TODO: understand why
  return true;
}

void OneCenterSlaterIntegral::computeHelperVariables() {
  n1_ = na_ + nb_;
  n2_ = nc_ + nd_;
  e1_ = ea_ + eb_;
  e2_ = ec_ + ed_;
}

double OneCenterSlaterIntegral::calculateFactor() {
  double x = std::pow(2 * ec_, nc_ + 0.5) * std::pow(2 * ed_, nd_ + 0.5) / std::sqrt(factorial(2 * nc_) * factorial(2 * nd_));
  double y = std::pow(2 * ea_, na_ + 0.5) * std::pow(2 * eb_, nb_ + 0.5) / std::sqrt(factorial(2 * na_) * factorial(2 * nb_));
  return x * y * factorial(n2_ + l_) / std::pow(e2_, n2_ + l_ + 1);
}

double OneCenterSlaterIntegral::calculateFirstTerm() {
  return factorial(n1_ - l_ - 1) / std::pow(e1_, n1_ - l_);
}

double OneCenterSlaterIntegral::calculateSecondTerm() {
  double s2 = 0.0;
  for (int ll = 1; ll <= n2_ + l_ + 1; ll++)
    s2 += calculateSecondSumTerm(ll);
  return s2;
}

double OneCenterSlaterIntegral::calculateSecondSumTerm(int ll) {
  double firstFraction = std::pow(e2_, n2_ + l_ - ll + 1) / factorial(n2_ + l_ - ll + 1);
  double secondFraction = factorial(n1_ + n2_ - ll) / std::pow(e1_ + e2_, n1_ + n2_ - ll + 1);
  return firstFraction * secondFraction;
}

double OneCenterSlaterIntegral::calculateThirdTerm() {
  double s3 = 0.0;
  for (int ll = 1; ll <= n2_ - l_; ll++)
    s3 += calculateThirdSumTerm(ll);
  return s3;
}

double OneCenterSlaterIntegral::calculateThirdSumTerm(int ll) {
  double numerator = std::pow(e2_, n2_ + l_ - ll + 1) * factorial(n2_ - l_ - 1) * factorial(n1_ + n2_ - ll);
  double denominator = factorial(n2_ + l_) * factorial(n2_ - l_ - ll) * std::pow(e1_ + e2_, n1_ + n2_ - ll + 1);
  return numerator / denominator;
}

long long OneCenterSlaterIntegral::factorial(int n) {
  assert(n <= 20);
  static auto factorialArray = createFactorialArrayUpTo20();
  return factorialArray[n];
}

std::array<long long int, 20> OneCenterSlaterIntegral::createFactorialArrayUpTo20() {
  std::array<long long int, 20> array{};
  array[0] = 1;
  for (long long int i = 1; i < 20; ++i) {
    array[i] = array[i - 1] * i;
  }
  return array;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
