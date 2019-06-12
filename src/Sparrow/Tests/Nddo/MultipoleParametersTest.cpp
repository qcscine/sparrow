/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <Utils/Constants.h>
#include <gmock/gmock.h>
#include <cmath>

namespace Scine {
namespace Sparrow {

#define SQRT2 1.41421356237

using namespace testing;
using namespace nddo;
using namespace multipole;

class AChargeSeparationParameter : public Test {
 public:
  ChargeSeparationParameter d;
  void SetUp() override {
  }
};

TEST_F(AChargeSeparationParameter, HasZeroAsDefaultValue) {
  ASSERT_THAT(d.get(sp1), Eq(0.0));
}

TEST_F(AChargeSeparationParameter, CanBeSet) {
  d.set(sp1, 1.0);
  ASSERT_THAT(d.get(sp1), Eq(1.0));
}

TEST_F(AChargeSeparationParameter, CanBeCalculatedFromExponentsForVanadium) {
  d.computeFromExponents(4, 4, 3, 1.974330, 1.063106, 1.394806);

  ASSERT_THAT(d.get(sp1), DoubleNear(1.11908712, 1e-8));
  ASSERT_THAT(d.get(pp2), DoubleNear(2.82191992 / SQRT2, 1e-8));
  ASSERT_THAT(d.get(sd2), DoubleNear(1.79796529 / SQRT2, 1e-8));
  ASSERT_THAT(d.get(pd1), DoubleNear(1.10444703, 1e-8));
  ASSERT_THAT(d.get(dd2), DoubleNear(1.43389116 / SQRT2, 1e-8));
}

TEST_F(AChargeSeparationParameter, CanBeCalculatedFromExponentsForHelium) {
  d.computeFromExponents(1, 2, 3.31320400, 3.65713300);

  ASSERT_THAT(d.get(sp1), DoubleNear(0.24758191, 1e-8));
  ASSERT_THAT(d.get(pp2), DoubleNear(0.29953658 / SQRT2, 1e-8));
}

TEST_F(AChargeSeparationParameter, CanBeCalculatedFromExponentsForXenon) {
  d.computeFromExponents(6, 5, 2.759787, 1.977446);

  ASSERT_THAT(d.get(sp1), DoubleNear(1.15147932, 1e-8));
  ASSERT_THAT(d.get(pp2), DoubleNear(1.83730955 / SQRT2, 1e-8));
}

class AKlopmanParameter : public Test {
 public:
  KlopmanParameter k;
  void SetUp() override {
  }
};

TEST_F(AKlopmanParameter, HasZeroAsDefaultValue) {
  ASSERT_THAT(k.get(sp1), Eq(0.0));
}

TEST_F(AKlopmanParameter, CanBeSet) {
  k.set(sp1, 1.0);
  ASSERT_THAT(k.get(sp1), Eq(1.0));
}

TEST_F(AKlopmanParameter, CanBeGotFromOrbitalPair) {
  multipolePair_t type = pairType(0, 1, 1);
  ASSERT_THAT(type, Eq(sp1));
}

TEST_F(AKlopmanParameter, CanBeGeneratedForH) {
  double gss = 14.44868600;
  k.generateUpToS(gss * Utils::Constants::hartree_per_ev);
  ASSERT_THAT(k.get(ss0), DoubleNear(0.94165599, 1e-6));
}

TEST_F(AKlopmanParameter, CanBeGeneratedForC) {
  double gss = 13.33551900;
  double hsp = 0.71732200;
  // hpp = 0.5 * (gpp - gp2)
  double hpp = 0.5 * (10.77832600 - 9.48621200);
  double D1sp = 0.75356425;
  double D2pp = 1.01715357 / sqrt(2);
  k.generateUpToP(gss * Utils::Constants::hartree_per_ev, hsp * Utils::Constants::hartree_per_ev, D1sp,
                  hpp * Utils::Constants::hartree_per_ev, D2pp);
  ASSERT_THAT(k.get(ss0), DoubleNear(1.02025963, 1e-6));
  ASSERT_THAT(k.get(pp0), DoubleNear(1.02025963, 1e-6));
  ASSERT_THAT(k.get(sp1), DoubleNear(1.29184422, 1e-6));
  ASSERT_THAT(k.get(pp2), DoubleNear(0.76267645, 1e-6));
}

TEST_F(AKlopmanParameter, CannotBeGeneratedForXe) {
  double gss = 20.000252000;
  double hsp = 4.418843000;
  // hpp = 0.5 * (gpp - gp2)
  double hpp = 0.5 * (2.30578700 - 4.06322000);
  double D1sp = 1.15147932;
  double D2pp = 1.837309557 / sqrt(2);
  k.generateUpToP(gss * Utils::Constants::hartree_per_ev, hsp * Utils::Constants::hartree_per_ev, D1sp,
                  hpp * Utils::Constants::hartree_per_ev, D2pp);
  ASSERT_THAT(k.get(ss0), DoubleNear(0.68027601, 1e-6));
  ASSERT_THAT(k.get(pp0), DoubleNear(0.68027601, 1e-6));
  ASSERT_THAT(k.get(sp1), DoubleNear(0.72182200, 1e-6));
  ASSERT_TRUE(std::isnan(k.get(pp2)) || std::isinf(k.get(pp2)));
}

TEST_F(AKlopmanParameter, CanBeGeneratedForCl) {
  double gss = 11.14265400;
  double hsp = 5.00426700;
  // hpp = 0.5 * (gpp - gp2)
  double hpp = 0.5 * (9.55188600 - 8.12843600);
  double D1sp = 0.81500430;
  double D2pp = 1.11721851 / sqrt(2);
  double f0dd = 45.02799881;
  double g1pd = 4.67265712;
  double D1pd = 0.75101197;
  double G2sd = 0.05252178;
  double D2sd = 1.10741013 / sqrt(2);
  double F2dd = 23.76887806;
  double D2dd = 1.51053637 / sqrt(2);
  k.generateUpToD(gss * Utils::Constants::hartree_per_ev, hsp * Utils::Constants::hartree_per_ev, D1sp,
                  hpp * Utils::Constants::hartree_per_ev, D2pp, f0dd * Utils::Constants::hartree_per_ev,
                  g1pd * Utils::Constants::hartree_per_ev, D1pd, G2sd * Utils::Constants::hartree_per_ev, D2sd,
                  F2dd * Utils::Constants::hartree_per_ev, D2dd);
  ASSERT_THAT(k.get(ss0), DoubleNear(1.22104587, 1e-6));
  ASSERT_THAT(k.get(pp0), DoubleNear(1.22104587, 1e-6));
  ASSERT_THAT(k.get(sp1), DoubleNear(0.57538280, 1e-6));
  ASSERT_THAT(k.get(pp2), DoubleNear(0.78649683, 1e-6));
  ASSERT_THAT(k.get(dd0), DoubleNear(0.30216070, 1e-6));
  ASSERT_THAT(k.get(pd1), DoubleNear(1.03732318, 1e-6));
  ASSERT_THAT(k.get(sd2), DoubleNear(2.34289326, 1e-6));
  ASSERT_THAT(k.get(dd2), DoubleNear(0.72430125, 1e-6));
}
} // namespace Sparrow
} // namespace Scine
