/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ONECENTERSLATERINTEGRAL_H
#define SPARROW_ONECENTERSLATERINTEGRAL_H

#include <array>

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * Calculation of the one-centre integrals (related to Slater-Condon parameters)
 * according to Kumar, Mishra, J. Phys., 1987, 29, 385-390.
 * The calculation is not optimized for performance as it does not
 * need to be fast. NB: this returns
 * U^l(a,b,c,d) = R^k(a,c,b,d) as defined in other papers.
 */

class OneCenterSlaterIntegral {
 public:
  OneCenterSlaterIntegral() = default;
  void setPrincipal(int na, int nb, int nc, int nd);
  void setAngular(int la, int lb, int lc, int ld);
  void setExponents(double ea, double eb, double ec, double ed);

  double calculate(int l);

 private:
  bool validLValue() const;
  void computeHelperVariables();
  double calculateFactor();
  double calculateFirstTerm();
  double calculateSecondTerm();
  double calculateThirdTerm();
  double calculateSecondSumTerm(int ll);
  double calculateThirdSumTerm(int ll);
  /*! Factorial, works up to n = 20. */
  long long int factorial(int n);
  static std::array<long long int, 20> createFactorialArrayUpTo20();

  int l_;
  int na_, nb_, nc_, nd_;
  int la_, lb_, lc_, ld_;
  double ea_, eb_, ec_, ed_;
  int n1_, n2_;    // NB: double of value in Ref.
  double e1_, e2_; // NB: double of value in Ref.
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ONECENTERSLATERINTEGRAL_H
