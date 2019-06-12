/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/SlaterCondonParameters.h>
#include <Utils/Constants.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class SlaterCondonParameterCalculation : public Test {
 public:
  const double arbitraryValue{1.23456};
  const sc_t arbitraryType{F0dd};
  SlaterCondonParameters par;

  void SetUp() override {
    par.setElement(Utils::ElementType::V);
    par.setExponents(1.094426, 0.755378, 1.099367);
  }
};

TEST_F(SlaterCondonParameterCalculation, DoesNotRecalculateGivenParameters) {
  par.set(arbitraryType, arbitraryValue);

  par.calculate();

  ASSERT_THAT(par.get(arbitraryType), Eq(arbitraryValue));
}

TEST_F(SlaterCondonParameterCalculation, GivesCorrectResultsForVanadium) {
  par.set(F0sd, 6.810021 / Utils::Constants::ev_per_hartree);
  par.set(G2sd, 1.831407 / Utils::Constants::ev_per_hartree);

  par.calculate();

  ASSERT_THAT(par.get(F0dd) * Utils::Constants::ev_per_hartree, DoubleNear(7.722275, 1e-5));
  ASSERT_THAT(par.get(F2dd) * Utils::Constants::ev_per_hartree, DoubleNear(4.076348, 1e-5));
  ASSERT_THAT(par.get(F4dd) * Utils::Constants::ev_per_hartree, DoubleNear(2.658488, 1e-5));
  ASSERT_THAT(par.get(F0sd) * Utils::Constants::ev_per_hartree, DoubleNear(6.810021, 1e-5));
  ASSERT_THAT(par.get(G2sd) * Utils::Constants::ev_per_hartree, DoubleNear(1.831407, 1e-5));
  ASSERT_THAT(par.get(F0pd) * Utils::Constants::ev_per_hartree, DoubleNear(4.947693, 1e-5));
  ASSERT_THAT(par.get(F2pd) * Utils::Constants::ev_per_hartree, DoubleNear(1.937683, 1e-5));
  ASSERT_THAT(par.get(G1pd) * Utils::Constants::ev_per_hartree, DoubleNear(1.851863, 1e-5));
  ASSERT_THAT(par.get(G3pd) * Utils::Constants::ev_per_hartree, DoubleNear(1.127754, 1e-5));
}
} // namespace Sparrow
} // namespace Scine
