/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/Local2c2eIntegralCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <Utils/Constants.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {
using namespace testing;
using namespace nddo;
using namespace multipole;
using namespace GeneralTypes;

class ALocal2c2eIntegralCalculator : public Test {
 public:
  using Calc = Local2c2eIntegralCalculator;
  Calc l;
  ChargeSeparationParameter DC, DO, DH, DV;
  KlopmanParameter KC, KO, KH, KV;

  void SetUp() override {
    KC.set(ss0, 1.02025963);
    KC.set(sp1, 1.29184422);
    KC.set(pp0, 1.02025963);
    KC.set(pp2, 0.76267645);
    KO.set(ss0, 1.20361298);
    KO.set(sp1, 0.29685467);
    KO.set(pp0, 1.20361298);
    KO.set(pp2, 0.45722970);
    KH.set(ss0, 0.94165599);
    KH.set(pp0, 0.94165599);
    KV.set(ss0, 2.27401429);
    KV.set(sp1, 1.49936943);
    KV.set(pp2, 1.99647866);
    KV.set(sd2, 1.34719346);
    KV.set(pd1, 1.88688065);
    KV.set(dd2, 1.29691630);
    KV.set(pp0, 2.27401429);
    KV.set(dd0, 1.76187597);
    DC.computeFromExponents(2, 2, 2.047558, 1.702841);
    DO.computeFromExponents(2, 2, 5.421751, 2.270960);
    DV.computeFromExponents(4, 4, 3, 1.974330, 1.063106, 1.394806);
  }
};

TEST_F(ALocal2c2eIntegralCalculator, GivesCorrectResultsForCOPair) {
  double R = 2.267671186;

  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(10, 10);

  mat(0, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(0, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(1, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(1, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_x, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(2, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(2, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_x, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(3, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(3, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_y, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(4, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(4, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_y, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(5, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(5, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_y, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(6, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(6, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_z, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(7, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(7, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::x_z, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(8, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(8, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::y_z, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(9, 0) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::s_s, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 1) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::s_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 2) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::x_x, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 3) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::s_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 4) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::x_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 5) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::y_y, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 6) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::s_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 7) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::x_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 8) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::y_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;
  mat(9, 9) =
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::z_z, twoElIntegral_t::z_z, R, DC, DO, KC, KO).value()) *
      Utils::Constants::ev_per_hartree;

  mat(4, 4) = 0.5 * (mat(2, 2) - mat(2, 5));

  Eigen::MatrixXd expectedMatrix(10, 10);
  expectedMatrix << 8.5674, 0, 8.2102, 0, 0, 8.2102, -0.8137, 0, 0, 8.9872, 0, 0.2039, 0, 0, 0, 0, 0, -0.1705, 0, 0,
      8.1222, 0, 7.9166, 0, 0, 7.8299, -0.649, 0, 0, 8.3368, 0, 0, 0, 0.2039, 0, 0, 0, 0, -0.1705, 0, 0, 0, 0, 0,
      0.0434, 0, 0, 0, 0, 0, 8.1222, 0, 7.8299, 0, 0, 7.9166, -0.649, 0, 0, 8.3368, 1.1892, 0, 1.0006, 0, 0, 1.0006,
      -0.2148, 0, 0, 1.1075, 0, 0.2142, 0, 0, 0, 0, 0, -0.2334, 0, 0, 0, 0, 0, 0.2142, 0, 0, 0, 0, -0.2334, 0, 8.8437,
      0, 8.2156, 0, 0, 8.2156, -0.9892, 0, 0, 8.8243;

  // std::cout << "The expected matrix is\n" << expectedMatrix <<
  //   "\nand the calculated matrix is\n" << mat << std::endl;

  for (long long i = 0; i < 10; i++)
    for (long long j = 0; j < 10; j++) {
      SCOPED_TRACE("... for the pair (" + std::to_string(i) + "," + std::to_string(j) + ")");
      EXPECT_THAT(mat(i, j), DoubleNear(expectedMatrix(i, j), 1e-3));
    }
}

TEST_F(ALocal2c2eIntegralCalculator, GivesCorrectResultsForHVPair) {
  double R = 2.267671186;

  EXPECT_THAT((Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_s, R, DH, DV, KH, KV)).value() *
                  Utils::Constants::ev_per_hartree,
              DoubleNear(6.9155, 1e-3));
  EXPECT_THAT((Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_x, R, DH, DV, KH, KV)).value() *
                  Utils::Constants::ev_per_hartree,
              DoubleNear(5.7462, 1e-3));
  EXPECT_THAT((Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_y, R, DH, DV, KH, KV)).value() *
                  Utils::Constants::ev_per_hartree,
              DoubleNear(5.7462, 1e-3));
  EXPECT_THAT((Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_z, R, DH, DV, KH, KV)).value() *
                  Utils::Constants::ev_per_hartree,
              DoubleNear(-1.7843, 1e-3));
  EXPECT_THAT((Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z_z, R, DH, DV, KH, KV)).value() *
                  Utils::Constants::ev_per_hartree,
              DoubleNear(6.2308, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_z2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.7876, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_xz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_yz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::s_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_z2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_xz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(-1.3574, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_z2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_yz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(-1.3574, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::y_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z_z2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(-1.5674, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z_xz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z_yz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z2_z2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(8.4137, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z2_xz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z2_yz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z2_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::z2_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::xz_xz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(8.0626, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::xz_yz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::xz_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::xz_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::yz_yz, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(8.0626, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::yz_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::yz_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(0.0000, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::x2y2_x2y2, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(7.0093, 1e-3));
  EXPECT_THAT(
      (Calc::getIntegral<Utils::derivOrder::one>(twoElIntegral_t::s_s, twoElIntegral_t::xy_xy, R, DH, DV, KH, KV)).value() *
          Utils::Constants::ev_per_hartree,
      DoubleNear(7.0093, 1e-3));
}
} // namespace Sparrow
} // namespace Scine
