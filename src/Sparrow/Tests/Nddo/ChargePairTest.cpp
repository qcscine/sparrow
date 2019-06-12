/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/MultipoleChargePair.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace multipole;

class AMultipoleChargePair : public Test {
 public:
  MultipoleCharge c1, c2;
  void SetUp() override {
    c1.x = ChargeDistance::d0;
    c1.y = ChargeDistance::d0;
    c1.z = ChargeDistance::d0;
    c1.q = 1;
    c2.x = ChargeDistance::d0;
    c2.y = ChargeDistance::d0;
    c2.z = ChargeDistance::d0;
    c2.q = 1;
  }
};

TEST_F(AMultipoleChargePair, CanBeSet) {
  MultipoleChargePair cp(c1, c2);
}

TEST_F(AMultipoleChargePair, GetsCorrectDistanceForChargesInCenter) {
  MultipoleChargePair cp(c1, c2);
  ASSERT_THAT(cp.getXDistance(), Eq(ChargeDistanceSeparation::d00));
}

TEST_F(AMultipoleChargePair, GetsCorrectDistancesForOtherCharges) {
  MultipoleCharge c3;
  c3.x = ChargeDistance::dp1;
  c3.y = ChargeDistance::dm1;
  c3.z = ChargeDistance::d0;
  c3.q = 0.5;
  MultipoleCharge c4;
  c4.x = ChargeDistance::dms2;
  c4.y = ChargeDistance::d0;
  c4.z = ChargeDistance::dms2;
  c4.q = -0.25;

  MultipoleChargePair cp(c1, c3);
  ASSERT_THAT(cp.getXDistance(), Eq(ChargeDistanceSeparation::d01));
  ASSERT_THAT(cp.getYDistance(), Eq(ChargeDistanceSeparation::d01));
  ASSERT_THAT(cp.getChargeProduct(), Eq(0.5));

  MultipoleChargePair cp1(c3, c3);
  ASSERT_THAT(cp1.getXDistance(), Eq(ChargeDistanceSeparation::m11));
  ASSERT_THAT(cp1.getYDistance(), Eq(ChargeDistanceSeparation::m11));
  ASSERT_THAT(cp1.getChargeProduct(), Eq(0.5 * 0.5));

  MultipoleChargePair cp2(c3, c4);
  ASSERT_THAT(cp2.getXDistance(), Eq(ChargeDistanceSeparation::p1s2));
  ASSERT_THAT(cp2.getYDistance(), Eq(ChargeDistanceSeparation::d10));
  ASSERT_THAT(cp2.getChargeProduct(), Eq(-0.5 * 0.25));
}

TEST_F(AMultipoleChargePair, AllowsToRetrieveTheCharges) {
  MultipoleChargePair cp(c1, c2);

  MultipoleCharge retrieved1 = cp.firstCharge();
  MultipoleCharge retrieved2 = cp.secondCharge();

  ASSERT_THAT(retrieved1.x, Eq(c1.x));
  ASSERT_THAT(retrieved1.y, Eq(c1.y));
  ASSERT_THAT(retrieved1.z, Eq(c1.z));
  ASSERT_THAT(retrieved1.q, Eq(c1.q));
  ASSERT_THAT(retrieved2.x, Eq(c2.x));
  ASSERT_THAT(retrieved2.y, Eq(c2.y));
  ASSERT_THAT(retrieved2.z, Eq(c2.z));
  ASSERT_THAT(retrieved2.q, Eq(c2.q));
}
} // namespace Sparrow
} // namespace Scine
