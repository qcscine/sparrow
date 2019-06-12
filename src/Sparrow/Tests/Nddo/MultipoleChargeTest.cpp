/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/MultipoleCharge.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace multipole;

class AMultipoleCharge : public Test {
 public:
  MultipoleCharge c;
  void SetUp() override {
  }
};

TEST_F(AMultipoleCharge, CanBeSet) {
  c.x = ChargeDistance::d0;
  c.y = ChargeDistance::d0;
  c.z = ChargeDistance::d0;
  c.q = 1;
}

TEST_F(AMultipoleCharge, CanBeSetInConstructor) {
  c.x = ChargeDistance::d0;
  c.y = ChargeDistance::d0;
  c.z = ChargeDistance::d0;
  c.q = 1;

  MultipoleCharge d(ChargeDistance::d0, ChargeDistance::d0, ChargeDistance::d0, 1);
}
} // namespace Sparrow
} // namespace Scine
