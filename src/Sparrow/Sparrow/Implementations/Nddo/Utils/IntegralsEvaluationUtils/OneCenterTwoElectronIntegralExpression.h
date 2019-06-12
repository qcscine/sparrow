/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ONECENTERTWOELECTRONINTEGRALEXPRESSION_H
#define SPARROW_ONECENTERTWOELECTRONINTEGRALEXPRESSION_H

#include "Sparrow/Implementations/Nddo/Utils/ParameterUtils/SlaterCondonParameters.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * This class allows the creation of instances containing, so to say,
 * the analytical expression for the calculation of a one-center
 * two-electron integral based on Slater-Type parameters.
 */
class OneCenterTwoElectronIntegralExpression {
 public:
  OneCenterTwoElectronIntegralExpression(double F1, sc_t E1, double F2, sc_t E2, double F3, sc_t E3)
    : f1(F1), f2(F2), f3(F3), e1(E1), e2(E2), e3(E3) {
  }

  OneCenterTwoElectronIntegralExpression(double F1, sc_t E1, double F2, sc_t E2)
    : OneCenterTwoElectronIntegralExpression(F1, E1, F2, E2, 0, F0ss) {
  }

  OneCenterTwoElectronIntegralExpression(double F1, sc_t E1) : OneCenterTwoElectronIntegralExpression(F1, E1, 0, F0ss) {
  }

  OneCenterTwoElectronIntegralExpression() : OneCenterTwoElectronIntegralExpression(0, F0ss) {
  }

  double result(const SlaterCondonParameters* p) {
    return f1 * p->get(e1) + f2 * p->get(e2) + f3 * p->get(e3);
  }

 private:
  double f1, f2, f3;
  sc_t e1, e2, e3;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ONECENTERTWOELECTRONINTEGRALEXPRESSION_H
