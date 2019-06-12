/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Utils/Geometry/ElementTypes.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class APairwiseParameters : public Test {
 public:
  PM6DiatomicParameters p;
  const double arbitraryAlpha{3.45678};
  const double arbitraryX{8.76543};

  void SetUp() override {
  }
};

TEST_F(APairwiseParameters, SetsDefaultElementsToNone) {
  ASSERT_THAT(p.firstElement(), Eq(Utils::ElementType::none));
  ASSERT_THAT(p.secondElement(), Eq(Utils::ElementType::none));
}

TEST_F(APairwiseParameters, CanSetAtomTypes) {
  p.setFirstElement(Utils::ElementType::Ga);
  p.setSecondElement(Utils::ElementType::He);

  ASSERT_THAT(p.firstElement(), Eq(Utils::ElementType::Ga));
  ASSERT_THAT(p.secondElement(), Eq(Utils::ElementType::He));
}

TEST_F(APairwiseParameters, IsNotValidIfNoElementsSet) {
  ASSERT_FALSE(p.isValid());

  p.setFirstElement(Utils::ElementType::Ga);
  ASSERT_FALSE(p.isValid());

  p.setSecondElement(Utils::ElementType::He);
  ASSERT_TRUE(p.isValid());

  p.setFirstElement(Utils::ElementType::none);
  ASSERT_FALSE(p.isValid());
}

TEST_F(APairwiseParameters, CanSetAlpha) {
  p.setAlpha(arbitraryAlpha);

  ASSERT_THAT(p.alpha(), Eq(arbitraryAlpha));
}

TEST_F(APairwiseParameters, CanSetX) {
  p.setX(arbitraryX);

  ASSERT_THAT(p.x(), Eq(arbitraryX));
}
} // namespace Sparrow
} // namespace Scine
