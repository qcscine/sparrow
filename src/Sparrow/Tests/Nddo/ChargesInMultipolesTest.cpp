/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/ChargesInMultipoles.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace multipole;

class AChargesInMultipoles : public Test {
 public:
  MultipoleCharge c1, c2;
  void SetUp() override {
    c1.x = ChargeDistance::d0;
    c1.y = ChargeDistance::d0;
    c1.z = ChargeDistance::d0;
    c1.q = 1;
    c2.x = ChargeDistance::dp1;
    c2.y = ChargeDistance::d0;
    c2.z = ChargeDistance::d0;
    c2.q = 1;
  }
};

TEST_F(AChargesInMultipoles, HasCorrectNumberOfChargesInMultipoles) {
  using CM = ChargesInMultipoles;

  ASSERT_THAT(CM::getChargeConfiguration(M00).size(), Eq(1u));
  ASSERT_THAT(CM::getChargeConfiguration(Qxx).size(), Eq(3u));
  ASSERT_THAT(CM::getChargeConfiguration(Qyy).size(), Eq(3u));
  ASSERT_THAT(CM::getChargeConfiguration(M1m1).size(), Eq(2u));
  ASSERT_THAT(CM::getChargeConfiguration(M10).size(), Eq(2u));
  ASSERT_THAT(CM::getChargeConfiguration(M11).size(), Eq(2u));
  ASSERT_THAT(CM::getChargeConfiguration(M2m2).size(), Eq(4u));
  ASSERT_THAT(CM::getChargeConfiguration(M2m1).size(), Eq(4u));
  ASSERT_THAT(CM::getChargeConfiguration(M20).size(), Eq(6u));
  ASSERT_THAT(CM::getChargeConfiguration(M21).size(), Eq(4u));
  ASSERT_THAT(CM::getChargeConfiguration(M22).size(), Eq(4u));
  ASSERT_THAT(CM::getChargeConfiguration(Qzx).size(), Eq(4u));
}
} // namespace Sparrow
} // namespace Scine
