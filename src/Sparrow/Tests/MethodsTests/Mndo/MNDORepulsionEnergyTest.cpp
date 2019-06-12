/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Utils/Constants.h"
#include "Utils/Geometry/ElementTypes.h"
#include <Sparrow/Implementations/Nddo/Mndo/MNDOPairwiseRepulsion.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <gmock/gmock.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;
using namespace Utils::AutomaticDifferentiation;
using std::sqrt;
using Utils::derivOrder;

class AMNDOPairwiseRepulsion : public Test {
 public:
  AMNDOPairwiseRepulsion() : ap1(arbitraryElement1), ap2(arbitraryElement2), rep(ap1, ap2) {
  }
  const Utils::ElementType arbitraryElement1{Utils::ElementType::C};
  const Utils::ElementType arbitraryElement2{Utils::ElementType::N};
  double arbitraryCharge1, arbitraryCharge2;
  double arbitraryAlpha1, arbitraryAlpha2;
  double arbitraryIntegral;
  double arbitraryRadius;
  double aCore, bCore;
  AtomicParameters ap1;
  AtomicParameters ap2;
  MNDOPairwiseRepulsion rep;

  void SetUp() override {
    arbitraryCharge1 = 4.50000;
    arbitraryCharge2 = 18.20000;
    arbitraryAlpha1 = 2.32323;
    arbitraryAlpha2 = 5.34624;
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
    ap1.setAlpha(arbitraryAlpha1);
    ap2.setAlpha(arbitraryAlpha2);
  }
};

TEST_F(AMNDOPairwiseRepulsion, ReturnsTheCorrectValueForArbitraryParameters) {
  double expected = arbitraryCharge1 * arbitraryCharge2 * arbitraryIntegral *
                    (1.0 + std::exp(-arbitraryAlpha1 * arbitraryRadius) + std::exp(-arbitraryAlpha2 * arbitraryRadius));
  ASSERT_THAT(rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}

/*TEST_F(APM6PairwiseRepulsion, UsesOtherFormulaForOH) {
  double FACTOR = 2; // factor that has to be taken to get the same results as MOPAC, but that is not in the theory...
  additionalTerm = 0.0000207389327971913490822259 *
                   std::pow((std::pow(8, 1.0 / 3.0) + std::pow(1, 1.0 / 3.0)) / arbitraryRadius, 12) /
Utils::ev_per_hartree; ap1.setElement(Utils::ElementType::H); ap2.setElement(Utils::ElementType::O);

  double expected = 1 * 6 * arbitraryIntegral *
                        (1.0 + FACTOR * arbitraryX * std::exp(-arbitraryAlpha * (arbitraryRadius * arbitraryRadius))) +
                    additionalTerm;
  ASSERT_THAT(rep.calculateRepulsion<derivOrder::zero>(arbitraryRadius), DoubleEq(expected));
}*/

TEST_F(AMNDOPairwiseRepulsion, ReturnsTheCorrectDerivativeForArbitraryParameters) {
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

  First1D ed = rep.calculateRepulsion<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(ed.value(), DoubleEq(expected));
  ASSERT_THAT(ed.derivative(), DoubleEq(derExp));
}

TEST_F(AMNDOPairwiseRepulsion, CorrectDerivativeForStandardParenthesis) {
  double expected = (1.0 + std::exp(-arbitraryAlpha1 * arbitraryRadius) + std::exp(-arbitraryAlpha2 * arbitraryRadius));
  double DerExpected = -arbitraryAlpha1 * std::exp(-arbitraryAlpha1 * arbitraryRadius) -
                       arbitraryAlpha2 * std::exp(-arbitraryAlpha2 * arbitraryRadius);

  First1D sT = rep.standardParenthesis<derivOrder::one>(arbitraryRadius);
  ASSERT_THAT(sT.value(), DoubleEq(expected));
  ASSERT_THAT(sT.derivative(), DoubleEq(DerExpected));
}
} // namespace Sparrow
} // namespace Scine
