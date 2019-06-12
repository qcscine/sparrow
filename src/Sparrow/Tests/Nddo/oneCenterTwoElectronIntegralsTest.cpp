/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementTypes.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class OneCenterTwoElectronIntegralCalculation : public Test {
 public:
  const Utils::ElementType arbitraryElement{Utils::ElementType::V};
  const double arbitraryIntegral{0.12345678};
  const double arbExp1{0.33445566}, arbExp2{0.989898}, arbExp3{1.111111};
  OneCenterTwoElectronIntegrals integral;
};

TEST_F(OneCenterTwoElectronIntegralCalculation, GivesCorrectValuesForH) {
  integral.setElement(Utils::ElementType::H);
  integral.set(0, 0, 0, 0, 14.448686 / Utils::Constants::ev_per_hartree);

  integral.calculateIntegrals();

  ASSERT_THAT(integral.getNumberIntegrals(), Eq(1));
  ASSERT_THAT(integral.get(0, 0, 0, 0), DoubleNear(14.448686 / Utils::Constants::ev_per_hartree, 1e-5));
}

TEST_F(OneCenterTwoElectronIntegralCalculation, GivesCorrectValuesForGa) {
  integral.setElement(Utils::ElementType::Ga);
  integral.set(0, 0, 0, 0, 10.354885 / Utils::Constants::ev_per_hartree);
  integral.set(0, 0, 1, 1, 7.993674 / Utils::Constants::ev_per_hartree);
  integral.set(1, 1, 1, 1, 6.090184 / Utils::Constants::ev_per_hartree);
  integral.set(1, 1, 2, 2, 6.299226 / Utils::Constants::ev_per_hartree);
  integral.set(0, 1, 0, 1, 1.295974 / Utils::Constants::ev_per_hartree);

  integral.calculateIntegrals();

  ASSERT_THAT(integral.getNumberIntegrals(), Eq(6));
  ASSERT_THAT(integral.get(0), DoubleNear(10.3549 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(1), DoubleNear(7.99367 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(2), DoubleNear(1.29597 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(3), DoubleNear(6.09018 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(4), DoubleNear(6.29923 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(5), DoubleNear(-0.104521 / Utils::Constants::ev_per_hartree, 1e-5));
}

TEST_F(OneCenterTwoElectronIntegralCalculation, GivesCorrectValuesForVanadium) {
  SlaterCondonParameters sc;
  sc.setElement(Utils::ElementType::V);
  sc.setExponents(1.094426, 0.755378, 1.099367);
  sc.set(F0sd, 6.810021 / Utils::Constants::ev_per_hartree);
  sc.set(G2sd, 1.831407 / Utils::Constants::ev_per_hartree);
  sc.calculate();

  integral.setSlaterCondonParameters(&sc);
  integral.setElement(Utils::ElementType::V);
  integral.set(0, 0, 0, 0, 5.983116 / Utils::Constants::ev_per_hartree);
  integral.set(0, 0, 1, 1, 4.736769 / Utils::Constants::ev_per_hartree);
  integral.set(1, 1, 1, 1, 4.499763 / Utils::Constants::ev_per_hartree);
  integral.set(1, 1, 2, 2, 3.944481 / Utils::Constants::ev_per_hartree);
  integral.set(0, 1, 0, 1, 0.901105 / Utils::Constants::ev_per_hartree);

  integral.calculateIntegrals();

  ASSERT_THAT(integral.getNumberIntegrals(), Eq(58));
  ASSERT_THAT(integral.get(0), DoubleNear(5.98312 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(1), DoubleNear(4.73677 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(2), DoubleNear(0.901105 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(3), DoubleNear(4.49976 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(4), DoubleNear(3.94448 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(5), DoubleNear(0.277641 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(6), DoubleNear(0.206313 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(7), DoubleNear(-0.237753 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(8), DoubleNear(0.273708 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(9), DoubleNear(-0.273708 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(10), DoubleNear(4.83697 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(11), DoubleNear(-0.191781 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(12), DoubleNear(0.237753 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(13), DoubleNear(0.618114 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(14), DoubleNear(0.0690462 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(15), DoubleNear(0.480847 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(16), DoubleNear(0.355915 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(17), DoubleNear(0.342754 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(18), DoubleNear(-0.0690462 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(19), DoubleNear(-0.118162 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(20), DoubleNear(-0.342754 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(21), DoubleNear(5.16914 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(22), DoubleNear(0.0958905 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(23), DoubleNear(4.72624 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(24), DoubleNear(0.166087 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(25), DoubleNear(-0.166087 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(26), DoubleNear(5.05842 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(27), DoubleNear(0.191781 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(28), DoubleNear(0.366281 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(29), DoubleNear(6.81002 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(30), DoubleNear(8.27206 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(31), DoubleNear(0.423188 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(32), DoubleNear(0.26404 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(33), DoubleNear(7.42569 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(34), DoubleNear(0.091884 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(35), DoubleNear(-0.091884 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(36), DoubleNear(7.74398 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(37), DoubleNear(-0.183768 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(38), DoubleNear(0.183768 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(39), DoubleNear(0.210991 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(40), DoubleNear(0.370139 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(41), DoubleNear(7.85008 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(42), DoubleNear(0.159148 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(43), DoubleNear(-0.159148 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(44), DoubleNear(7.53178 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(45), DoubleNear(0.450563 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(46), DoubleNear(-0.450563 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(47), DoubleNear(0.225281 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(48), DoubleNear(0.390199 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(49), DoubleNear(-0.390199 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(50), DoubleNear(-0.182068 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(51), DoubleNear(-0.569124 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(52), DoubleNear(-0.315351 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(53), DoubleNear(0.657168 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(54), DoubleNear(0.569124 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(55), DoubleNear(0.364136 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(56), DoubleNear(0.315351 / Utils::Constants::ev_per_hartree, 1e-5));
  ASSERT_THAT(integral.get(57), DoubleNear(-0.328584 / Utils::Constants::ev_per_hartree, 1e-5));
}
} // namespace Sparrow
} // namespace Scine
