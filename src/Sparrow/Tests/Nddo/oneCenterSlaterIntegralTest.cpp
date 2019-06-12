/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterSlaterIntegral.h>
#include <gmock/gmock.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class OneCenterSlaterIntegralCalculation : public Test {
 public:
  OneCenterSlaterIntegral integral;
};

TEST_F(OneCenterSlaterIntegralCalculation, ReturnsCorrectValueForF01s1s) {
  double sExp = 1.0;
  integral.setPrincipal(1, 1, 1, 1);
  integral.setAngular(0, 0, 0, 0);
  integral.setExponents(sExp, sExp, sExp, sExp);

  ASSERT_THAT(integral.calculate(0), DoubleNear(0.625 * sExp, 1e-4));
}

TEST_F(OneCenterSlaterIntegralCalculation, ReturnsCorrectValueForG21s3d) {
  double sExp = 0.22;
  double dExp = 14.21;
  integral.setPrincipal(1, 3, 3, 1);
  integral.setAngular(0, 2, 2, 0);
  integral.setExponents(sExp, dExp, dExp, sExp);

  double expected = 36.000026 * std::pow(sExp, 3) * std::pow(dExp, 7) / std::pow(sExp + dExp, 9);
  ASSERT_THAT(integral.calculate(2), DoubleNear(expected, 1e-4));
}
} // namespace Sparrow
} // namespace Scine
