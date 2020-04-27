/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class AAtomicParameters : public Test {
 public:
  AtomicParameters p;
  const double arbitraryCharge{4.0};

  void SetUp() override {
  }
};

TEST_F(AAtomicParameters, SetsDefaultElementToNone) {
  ASSERT_THAT(p.element(), Eq(Utils::ElementType::none));
}

TEST_F(AAtomicParameters, CanSetAtomType) {
  p.setElement(Utils::ElementType::Ga);

  ASSERT_THAT(p.element(), Eq(Utils::ElementType::Ga));
}

TEST_F(AAtomicParameters, SetsDefaultNumberOfAOsToZero) {
  ASSERT_THAT(p.nAOs(), Eq(0));
}

TEST_F(AAtomicParameters, CanSetNumberOfAOs) {
  p.setNAOs(4);
  ASSERT_THAT(p.nAOs(), Eq(4));
}

TEST_F(AAtomicParameters, SetsNumberOfAOsWhenElementTypeIsSet) {
  p.setElement(Utils::ElementType::V);
  ASSERT_THAT(p.nAOs(), Eq(9));
  p.setElement(Utils::ElementType::H);
  ASSERT_THAT(p.nAOs(), Eq(1));
  ASSERT_THAT(AtomicParameters(Utils::ElementType::V).nAOs(), Eq(9));
}

TEST_F(AAtomicParameters, MonoisotopicElementsWork) {
  p.setElement(Utils::ElementType::Al);
  // Al has Z = 13, 10 core electrons, so 3 core charge
  ASSERT_THAT(p.coreCharge(), Eq(3.0));
  p.setElement(Utils::ElementType::F);
  // F has Z = 9, 2 core electrons, so 7 core charge
  ASSERT_THAT(p.coreCharge(), Eq(7.0));
}

TEST_F(AAtomicParameters, CanSetCoreCharge) {
  p.setCoreCharge(arbitraryCharge);

  ASSERT_THAT(p.coreCharge(), Eq(arbitraryCharge));
}

TEST_F(AAtomicParameters, IsNotValidIfNoElementSet) {
  ASSERT_FALSE(p.isValid());

  p.setElement(Utils::ElementType::Ga);
  ASSERT_TRUE(p.isValid());

  p.setElement(Utils::ElementType::none);
  ASSERT_FALSE(p.isValid());
}

TEST_F(AAtomicParameters, CalculatesThirdRootOfCoreCharge) {
  p.setCoreCharge(arbitraryCharge);

  ASSERT_THAT(p.cubicRootOfCoreCharge(), Eq(std::pow(arbitraryCharge, 1.0 / 3.0)));
}

TEST_F(AAtomicParameters, HasNoGaussianRepulsionByDefault) {
  ASSERT_FALSE(p.hasGaussianRepulsionParameters());
}

TEST_F(AAtomicParameters, GaussianRepulsionParametersCanBeSet) {
  double a = 4.234234, b = 92.94324, c = 0.231121;
  p.clearGaussianRepulsionParameters();
  p.addGaussianRepulsionParameters(a, b, c);

  ASSERT_TRUE(p.hasGaussianRepulsionParameters());
  ASSERT_THAT(std::get<0>(p.getGaussianRepulsionParameters()[0]), Eq(a));
  ASSERT_THAT(std::get<1>(p.getGaussianRepulsionParameters()[0]), Eq(b));
  ASSERT_THAT(std::get<2>(p.getGaussianRepulsionParameters()[0]), Eq(c));
}

TEST_F(AAtomicParameters, CorrectNumberOfAOsForDifferentBasisSets) {
  p.setElement(Utils::ElementType::S, BasisFunctions::sp);
  ASSERT_EQ(p.nAOs(), 4);
  p.setElement(Utils::ElementType::S, BasisFunctions::spd);
  ASSERT_EQ(p.nAOs(), 9);
}

TEST_F(AAtomicParameters, CorrectCoreChargeForDifferentBasisSets) {
  p.setElement(Utils::ElementType::S, BasisFunctions::sp);
  ASSERT_EQ(p.coreCharge(), 6);
  p.setElement(Utils::ElementType::S, BasisFunctions::spd);
  ASSERT_EQ(p.coreCharge(), 6);
}

} // namespace Sparrow
} // namespace Scine
