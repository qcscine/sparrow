/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Pm6/PM6PairwiseRepulsion.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
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
using Utils::derivOrder;

class APM6PairwiseRepulsion : public Test {
 public:
  APM6PairwiseRepulsion()
    : ap1(arbitraryElement1), ap2(arbitraryElement2), pp(arbitraryElement1, arbitraryElement2), rep(ap1, ap2, pp) {
  }
  const Utils::ElementType arbitraryElement1{Utils::ElementType::Ga};
  const Utils::ElementType arbitraryElement2{Utils::ElementType::He};
  double arbitraryCharge1, arbitraryCharge2;
  double arbitraryAlpha;
  double arbitraryX;
  double arbitraryIntegral;
  double arbitraryRadius;
  double aCore, bCore;
  AtomicParameters ap1;
  AtomicParameters ap2;
  PM6DiatomicParameters pp;
  double additionalTerm;
  PM6PairwiseRepulsion rep;

  void SetUp() override {
    arbitraryCharge1 = 4.50000;
    arbitraryCharge2 = 18.20000;
    arbitraryAlpha = 2.32323;
    arbitraryX = 1.010101010;
    arbitraryIntegral = 0.123456789;
    arbitraryRadius = 1.321;

    aCore = 0.122;
    bCore = 0.922;
    double pSum = aCore + bCore;
    arbitraryIntegral = 1.0 / sqrt(arbitraryRadius * arbitraryRadius + pSum * pSum);
    ap1.setCoreCharge(arbitraryCharge1);
    ap1.setPCore(aCore);
    ap2.setCoreCharge(arbitraryCharge2);
    ap2.setPCore(bCore);
    pp.setAlpha(arbitraryAlpha);
    pp.setX(arbitraryX);
  }
};

TEST_F(APM6PairwiseRepulsion, ReturnsTheCorrectValueForArbitraryParameters) {
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(31, 1.0 / 3.0) + std::pow(2, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  double expected = arbitraryCharge1 * arbitraryCharge2 * arbitraryIntegral *
                        (1.0 + FACTOR * arbitraryX *
                                   std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 *
                                                                                     std::pow(arbitraryRadius, 6)))) +
                    additionalTerm;
  ASSERT_THAT(rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}

TEST_F(APM6PairwiseRepulsion, UsesOtherFormulaForOH) {
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(8, 1.0 / 3.0) + std::pow(1, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  ap1.setElement(Utils::ElementType::H);
  ap2.setElement(Utils::ElementType::O);

  double expected = 1 * 6 * arbitraryIntegral *
                        (1.0 + FACTOR * arbitraryX * std::exp(-arbitraryAlpha * (arbitraryRadius * arbitraryRadius))) +
                    additionalTerm;
  ASSERT_THAT(rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}

TEST_F(APM6PairwiseRepulsion, UsesOtherFormulaForCC) {
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(6, 1.0 / 3.0) + std::pow(6, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  ap1.setElement(Utils::ElementType::C);
  ap2.setElement(Utils::ElementType::C);

  double expected =
      4 * 4 * arbitraryIntegral *
          (1.0 +
           FACTOR * arbitraryX *
               std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 * std::pow(arbitraryRadius, 6))) +
           9.28 * std::exp(-3.1644797213016 * arbitraryRadius)) +
      additionalTerm;
  ASSERT_THAT(rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}

TEST_F(APM6PairwiseRepulsion, UsesOtherFormulaForSiO) {
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(8, 1.0 / 3.0) + std::pow(14, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  ap1.setElement(Utils::ElementType::Si);
  ap2.setElement(Utils::ElementType::O);

  double expected =
      4 * 6 * arbitraryIntegral *
          (1.0 +
           FACTOR * arbitraryX *
               std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 * std::pow(arbitraryRadius, 6))) -
           0.0007 * std::exp(-(arbitraryRadius - 5.480205761238679760340424) * (arbitraryRadius - 5.480205761238679760340424))) +
      additionalTerm;
  ASSERT_THAT((rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius)), DoubleEq(expected));
}

TEST_F(APM6PairwiseRepulsion, IncludesGaussianRepulsionIfParametersGiven) {
  double arbitraryA = 2.1231, arbitraryB = 0.23552, arbitraryC = 5.001;
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(2, 1.0 / 3.0) + std::pow(31, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  double expected =
      arbitraryCharge1 * arbitraryCharge2 * arbitraryIntegral *
          (1.0 + FACTOR * arbitraryX *
                     std::exp(-arbitraryAlpha *
                              (arbitraryRadius + 0.00001244878365801758178693335 * std::pow(arbitraryRadius, 6)))) +
      additionalTerm +
      arbitraryCharge1 * arbitraryCharge2 / arbitraryRadius *
          (arbitraryA * exp(-arbitraryB * (arbitraryRadius - arbitraryC) * (arbitraryRadius - arbitraryC))) / ev_per_hartree;
  ap1.clearGaussianRepulsionParameters();
  ap1.addGaussianRepulsionParameters(arbitraryA, arbitraryB, arbitraryC);

  ASSERT_THAT(rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}

TEST_F(APM6PairwiseRepulsion, ReturnsTheCorrectDerivativeForArbitraryParameters) {
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  double s = (aCore + bCore) * (aCore + bCore);
  First1D R2(arbitraryRadius * arbitraryRadius, 2 * arbitraryRadius);
  First1D integral = First1D(1.0, 0) / sqrt(R2 + First1D(s, 0));
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(31, 1.0 / 3.0) + std::pow(2, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  double expected = arbitraryCharge1 * arbitraryCharge2 * integral.value() *
                        (1.0 + FACTOR * arbitraryX *
                                   std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 *
                                                                                     std::pow(arbitraryRadius, 6)))) +
                    additionalTerm;

  double derAddTerm = (-12.0 / arbitraryRadius) * 0.0000207389327971913490822259 *
                      std::pow((std::pow(31, 1.0 / 3.0) + std::pow(2, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  double derExp =
      arbitraryCharge1 * arbitraryCharge2 * integral.derivative() *
          (1.0 + FACTOR * arbitraryX *
                     std::exp(-arbitraryAlpha *
                              (arbitraryRadius + 0.00001244878365801758178693335 * std::pow(arbitraryRadius, 6)))) +
      arbitraryCharge1 * arbitraryCharge2 * integral.value() *
          (FACTOR * arbitraryX *
           std::exp(-arbitraryAlpha * (arbitraryRadius + 0.00001244878365801758178693335 * std::pow(arbitraryRadius, 6)))) *
          (-arbitraryAlpha * (1 + 0.00001244878365801758178693335 * 6 * std::pow(arbitraryRadius, 5))) +
      derAddTerm;

  First1D ed = rep.calculateRepulsion<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(ed.value(), DoubleEq(expected));
  ASSERT_THAT(ed.derivative(), DoubleEq(derExp));
}

TEST_F(APM6PairwiseRepulsion, CorrectDerivativeForAdditionalTerm) {
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(31, 1.0 / 3.0) + std::pow(2, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;
  double derAddTerm = (-12.0 / arbitraryRadius) * 0.0000207389327971913490822259 *
                      std::pow((std::pow(31, 1.0 / 3.0) + std::pow(2, 1.0 / 3.0)) / arbitraryRadius, 12) / ev_per_hartree;

  First1D addT = rep.additionalTerm<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(addT.value(), DoubleEq(additionalTerm));
  ASSERT_THAT(addT.derivative(), DoubleEq(derAddTerm));
}

TEST_F(APM6PairwiseRepulsion, CorrectDerivativeForGaussianRepulsionIfParametersGiven) {
  double arbitraryA = 2.1231, arbitraryB = 0.23552, arbitraryC = 5.001;

  double expected = arbitraryCharge1 * arbitraryCharge2 / arbitraryRadius *
                    (arbitraryA * exp(-arbitraryB * (arbitraryRadius - arbitraryC) * (arbitraryRadius - arbitraryC))) /
                    ev_per_hartree;
  double DerExpected =
      -arbitraryCharge1 * arbitraryCharge2 / (arbitraryRadius * arbitraryRadius) *
          (arbitraryA * exp(-arbitraryB * (arbitraryRadius - arbitraryC) * (arbitraryRadius - arbitraryC))) / ev_per_hartree +
      arbitraryCharge1 * arbitraryCharge2 / arbitraryRadius *
          (arbitraryA * exp(-arbitraryB * (arbitraryRadius - arbitraryC) * (arbitraryRadius - arbitraryC)) *
           (-2 * arbitraryB * (arbitraryRadius - arbitraryC))) /
          ev_per_hartree;
  ap1.clearGaussianRepulsionParameters();
  ap1.addGaussianRepulsionParameters(arbitraryA, arbitraryB, arbitraryC);

  First1D GT = rep.gaussianRepulsionTerm<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(GT.value(), DoubleEq(expected));
  ASSERT_THAT(GT.derivative(), DoubleEq(DerExpected));

  // Check that ok if two times
  ap2.clearGaussianRepulsionParameters();
  ap2.addGaussianRepulsionParameters(arbitraryA, arbitraryB, arbitraryC);
  First1D GT2 = rep.gaussianRepulsionTerm<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(GT2.value(), DoubleEq(2 * expected));
  ASSERT_THAT(GT2.derivative(), DoubleEq(2 * DerExpected));
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
}
} // namespace Sparrow
} // namespace Scine
