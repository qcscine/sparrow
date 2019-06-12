/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Am1/AM1PairwiseRepulsion.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementTypes.h>
#include <gmock/gmock.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace Utils::AutomaticDifferentiation;
using namespace Utils::Constants;
using std::exp;
using std::sqrt;

class AAM1PairwiseRepulsion : public Test {
 public:
  AAM1PairwiseRepulsion() : ap1(arbitraryElement1), ap2(arbitraryElement2), rep(ap1, ap2) {
  }
  const Utils::ElementType arbitraryElement1{Utils::ElementType::P};
  const Utils::ElementType arbitraryElement2{Utils::ElementType::He};
  double arbitraryCharge1, arbitraryCharge2;
  double arbitraryAlpha1, arbitraryAlpha2;
  double arbitraryK1, arbitraryL1, arbitraryM1;
  double arbitraryK2, arbitraryL2, arbitraryM2;
  double arbitraryK3, arbitraryL3, arbitraryM3;
  double arbitraryIntegral;
  double arbitraryRadius;
  double aCore, bCore;
  AtomicParameters ap1;
  AtomicParameters ap2;
  AM1PairwiseRepulsion rep;

  void SetUp() override {
    arbitraryCharge1 = 4.50000;
    arbitraryCharge2 = 18.20000;
    arbitraryAlpha1 = 2.32323;
    arbitraryAlpha2 = 5.32543;
    arbitraryK1 = 0.235243, arbitraryL1 = 4.24523, arbitraryM1 = 1.772342;
    arbitraryK2 = -1.21423, arbitraryL2 = 3.25123, arbitraryM2 = 3.512213;
    arbitraryK3 = 0.51356, arbitraryL3 = 1.869943, arbitraryM3 = 7.534561;
    arbitraryIntegral = 0.123456789;
    arbitraryRadius = 1.321;

    aCore = 0.122;
    bCore = 0.922;
    double pSum = aCore + bCore;
    arbitraryIntegral = 1.0 / sqrt(arbitraryRadius * arbitraryRadius + pSum * pSum);
    ap1.setCoreCharge(arbitraryCharge1);
    ap1.setPCore(aCore);
    ap1.setAlpha(arbitraryAlpha1);
    ap2.setCoreCharge(arbitraryCharge2);
    ap2.setPCore(bCore);
    ap2.setAlpha(arbitraryAlpha2);
  }
};

TEST_F(AAM1PairwiseRepulsion, ReturnsTheCorrectValueForArbitraryParameters) {
  double expected =
      arbitraryCharge1 * arbitraryCharge2 * arbitraryIntegral *
          (1.0 + std::exp(-arbitraryAlpha1 * arbitraryRadius) + std::exp(-arbitraryAlpha2 * arbitraryRadius)) +
      arbitraryCharge1 * arbitraryCharge2 / arbitraryRadius *
          (arbitraryK1 * exp(-arbitraryL1 * (arbitraryRadius - arbitraryM1) * (arbitraryRadius - arbitraryM1)) +
           arbitraryK2 * exp(-arbitraryL2 * (arbitraryRadius - arbitraryM2) * (arbitraryRadius - arbitraryM2)) +
           arbitraryK3 * exp(-arbitraryL3 * (arbitraryRadius - arbitraryM3) * (arbitraryRadius - arbitraryM3))) /
          ev_per_hartree;

  ap1.clearGaussianRepulsionParameters();
  ap1.addGaussianRepulsionParameters(arbitraryK1, arbitraryL1, arbitraryM1);
  ap1.addGaussianRepulsionParameters(arbitraryK2, arbitraryL2, arbitraryM2);
  ap2.clearGaussianRepulsionParameters();
  ap2.addGaussianRepulsionParameters(arbitraryK3, arbitraryL3, arbitraryM3);
  ASSERT_THAT(rep.calculateRepulsion<Utils::derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}

TEST_F(AAM1PairwiseRepulsion, ReturnsTheCorrectDerivativeForArbitraryParameters) {
  double s = (aCore + bCore) * (aCore + bCore);
  First1D R2(arbitraryRadius * arbitraryRadius, 2 * arbitraryRadius);
  First1D integral = First1D(1.0, 0) / sqrt(R2 + First1D(s, 0));
  double expected = arbitraryCharge1 * arbitraryCharge2 * integral.value() *
                    (1.0 + std::exp(-arbitraryAlpha1 * arbitraryRadius) + std::exp(-arbitraryAlpha2 * arbitraryRadius));

  double derExp = arbitraryCharge1 * arbitraryCharge2 * integral.derivative() *
                      (1.0 + std::exp(-arbitraryAlpha1 * arbitraryRadius) + std::exp(-arbitraryAlpha2 * arbitraryRadius)) +
                  arbitraryCharge1 * arbitraryCharge2 * integral.value() *
                      (-arbitraryAlpha1 * std::exp(-arbitraryAlpha1 * arbitraryRadius) -
                       arbitraryAlpha2 * std::exp(-arbitraryAlpha2 * arbitraryRadius));

  First1D ed = rep.calculateRepulsion<Utils::derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(ed.value(), DoubleEq(expected));
  ASSERT_THAT(ed.derivative(), DoubleEq(derExp));
}
/*
TEST_F(AAM1PairwiseRepulsion, CorrectDerivativeForGaussianRepulsionIfParametersGiven) {
  double expected = arbitraryCharge1 * arbitraryCharge2 / arbitraryRadius *
                    (arbitraryK1 * exp(-arbitraryL1 * (arbitraryRadius - arbitraryM1) * (arbitraryRadius -
arbitraryM1))+ arbitraryK2 * exp(-arbitraryL2 * (arbitraryRadius - arbitraryM2) * (arbitraryRadius - arbitraryM2))+
                     arbitraryK3 * exp(-arbitraryL3 * (arbitraryRadius - arbitraryM3) * (arbitraryRadius -
arbitraryM3))) / ev_per_hartree; double DerExpected = -arbitraryCharge1 * arbitraryCharge2 / (arbitraryRadius *
arbitraryRadius) * (arbitraryK1 * exp(-arbitraryL1 * (arbitraryRadius - arbitraryM1) * (arbitraryRadius - arbitraryM1))+
                     arbitraryK2 * exp(-arbitraryL2 * (arbitraryRadius - arbitraryM2) * (arbitraryRadius -
arbitraryM2))+ arbitraryK3 * exp(-arbitraryL3 * (arbitraryRadius - arbitraryM3) * (arbitraryRadius - arbitraryM3)))/
          ev_per_hartree +
      arbitraryCharge1 * arbitraryCharge2 / arbitraryRadius *
          (arbitraryK1 * exp(-arbitraryL1 * (arbitraryRadius - arbitraryM1) * (arbitraryRadius - arbitraryM1))+
                     arbitraryK2 * exp(-arbitraryL2 * (arbitraryRadius - arbitraryM2) * (arbitraryRadius -
arbitraryM2))+ arbitraryK3 * exp(-arbitraryL3 * (arbitraryRadius - arbitraryM3) * (arbitraryRadius - arbitraryM3))) *
           (-2 * arbitraryL1 * (arbitraryRadius - arbitraryM1) -2 * arbitraryL2 * (arbitraryRadius - arbitraryM2) -2 *
arbitraryL3 * (arbitraryRadius - arbitraryM3)) / ev_per_hartree;

  ap1.clearGaussianRepulsionParameters();
  ap1.addGaussianRepulsionParameters(arbitraryK1, arbitraryL1, arbitraryM1);
  ap1.addGaussianRepulsionParameters(arbitraryK2, arbitraryL2, arbitraryM2);
  ap2.clearGaussianRepulsionParameters();
  ap2.addGaussianRepulsionParameters(arbitraryK3, arbitraryL3, arbitraryM3);
  First1D GT = rep.gaussianRepulsionTerm<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(GT.value(), DoubleEq(expected));
  ASSERT_THAT(GT.derivative(), DoubleEq(DerExpected));

}

TEST_F(APM6PairwiseRepulsion, CorrectDerivativeForStandardParenthesis) {
  double FACTOR = 2;
  double expected = (1.0 + FACTOR * arbitraryX *
                               std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 *
                                                                                 std::pow(arbitraryRadius, 6))));
  double DerExpected =
      FACTOR * arbitraryX *
      std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 * std::pow(arbitraryRadius, 6))) *
      (-arbitraryAlpha * (1 + 0.00001244878365801758178693335 * 6 * std::pow(arbitraryRadius, 5)));

  First1D sT = rep.standardParenthesis<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(sT.value(), DoubleEq(expected));
  ASSERT_THAT(sT.derivative(), DoubleEq(DerExpected));
} */
} // namespace Sparrow
} // namespace Scine
