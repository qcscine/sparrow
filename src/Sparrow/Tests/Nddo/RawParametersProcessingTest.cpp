/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParameterProcessor.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParametersContainer.h>
#include <Utils/Constants.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

#define SQRT2 1.41421356237309504880

using namespace testing;
using namespace nddo;

class ARawParameterProcessing : public Test {
 public:
  ARawParameterProcessing() : processor(cont) {
  }

  RawAtomicParameters rawAtomicPar;
  RawDiatomicParameters rawDiatomicPar;
  std::unique_ptr<AtomicParameters> atomicPar;
  std::shared_ptr<PM6DiatomicParameters> diatomicPar;
  std::shared_ptr<OneCenterTwoElectronIntegrals> integral;
  RawParametersContainer cont;
  RawParameterProcessor processor;

  void SetUp() override {
    // Those are the parameters for iodine
    rawAtomicPar.uss = -59.9732;
    rawAtomicPar.upp = -56.4598;
    rawAtomicPar.udd = -28.8226;
    rawAtomicPar.bs = -30.5225;
    rawAtomicPar.bp = -5.94212;
    rawAtomicPar.bd = -7.67611;
    rawAtomicPar.gss = 7.23476;
    rawAtomicPar.gpp = 9.87747;
    rawAtomicPar.gsp = 9.15441;
    rawAtomicPar.gp2 = 8.03592;
    rawAtomicPar.hsp = 5.00422;
    rawAtomicPar.zs = 4.49865;
    rawAtomicPar.zp = 1.91707;
    rawAtomicPar.zd = 1.87518;
    rawAtomicPar.zsn = 9.13524;
    rawAtomicPar.zpn = 6.88819;
    rawAtomicPar.zdn = 3.79152;
    rawAtomicPar.pcore = 1.8806;
    rawAtomicPar.f0sd = 20.495;
    rawAtomicPar.g2sd = 2.16014;
    // rawAtomicPar.grepa = -0.035519;
    // rawAtomicPar.grepb = 1.74439;
    // rawAtomicPar.grepc = 1.22384;

    // Those are parameters for the Cd-Cl pair:
    rawDiatomicPar.exponent = 0.943547;
    rawDiatomicPar.factor = 0.140424;

    cont.setAtomicParameters(Utils::ElementType::I, rawAtomicPar);
    cont.setDiatomicParameters(Utils::ElementType::Cd, Utils::ElementType::Cl, rawDiatomicPar);
    auto par = processor.processAtomicParameters(Utils::ElementType::I);
    atomicPar = std::move(par.first);
    integral = std::move(par.second);
    diatomicPar = processor.runtimeDiatomicParameters(Utils::ElementType::Cd, Utils::ElementType::Cl);
  }
};

TEST_F(ARawParameterProcessing, SetsBetaValuesCorrectly) {
  ASSERT_THAT(atomicPar->betaS(), DoubleEq(rawAtomicPar.bs * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->betaP(), DoubleEq(rawAtomicPar.bp * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->betaD(), DoubleEq(rawAtomicPar.bd * Utils::Constants::hartree_per_ev));
}

TEST_F(ARawParameterProcessing, SetsUValuesCorrectly) {
  ASSERT_THAT(atomicPar->Uss(), DoubleEq(rawAtomicPar.uss * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->Upp(), DoubleEq(rawAtomicPar.upp * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->Udd(), DoubleEq(rawAtomicPar.udd * Utils::Constants::hartree_per_ev));
}

TEST_F(ARawParameterProcessing, SetsChargeSeparationsCorrectly) {
  auto d = atomicPar->chargeSeparations();

  ASSERT_THAT(d.get(multipole::sp1), DoubleNear(0.3746525, 1e-4));
  ASSERT_THAT(d.get(multipole::pp2), DoubleNear(1.89517161 / SQRT2, 1e-5));
  ASSERT_THAT(d.get(multipole::sd2), DoubleNear(0.77747018 / SQRT2, 1e-5));
  ASSERT_THAT(d.get(multipole::pd1), DoubleNear(1.29634171, 1e-5));
  ASSERT_THAT(d.get(multipole::dd2), DoubleNear(1.63749938 / SQRT2, 1e-5));
}

TEST_F(ARawParameterProcessing, GtoExpansionHasBeenPerformed) {
  ASSERT_TRUE(atomicPar->GTOs().hasS());
  ASSERT_TRUE(atomicPar->GTOs().hasP());
  ASSERT_TRUE(atomicPar->GTOs().hasD());
}

/*TEST_F(ARawParameterProcessing, SetsGaussianRepulsionParameters) {
  auto g = atomicPar->getGaussianRepulsionParameters();

  ASSERT_TRUE(atomicPar->hasGaussianRepulsionParameters());
  ASSERT_THAT(std::get<0>(g), DoubleEq(rawAtomicPar.grepa * Utils::Constants::bohr_per_angstrom));
  ASSERT_THAT(std::get<1>(g), DoubleEq(rawAtomicPar.grepb * Utils::Constants::angstrom_per_bohr *
Utils::Constants::angstrom_per_bohr)); ASSERT_THAT(std::get<2>(g), DoubleEq(rawAtomicPar.grepc *
Utils::Constants::bohr_per_angstrom));
}*/

TEST_F(ARawParameterProcessing, SetsDiatomicExponentCorrectly) {
  ASSERT_THAT(diatomicPar->alpha(), Eq(rawDiatomicPar.exponent * Utils::Constants::angstrom_per_bohr));
}

TEST_F(ARawParameterProcessing, SetsDiatomicFactorCorrectly) {
  ASSERT_THAT(diatomicPar->x(), Eq(rawDiatomicPar.factor));
}

TEST_F(ARawParameterProcessing, SetsPCoreCorrectly) {
  ASSERT_THAT(atomicPar->pCore(), DoubleEq(rawAtomicPar.pcore));
}

TEST_F(ARawParameterProcessing, Returns1c2eIntegrals) {
  ASSERT_THAT(integral->getNumberIntegrals(), Eq(58));
  ASSERT_THAT(integral->get(0) * Utils::Constants::ev_per_hartree, DoubleNear(7.23480, 1e-4));
  ASSERT_THAT(integral->get(1) * Utils::Constants::ev_per_hartree, DoubleNear(9.15440, 1e-4));
  ASSERT_THAT(integral->get(2) * Utils::Constants::ev_per_hartree, DoubleNear(5.004200, 1e-4));
  ASSERT_THAT(integral->get(3), DoubleNear(9.87750 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(4), DoubleNear(8.03590 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(5), DoubleNear(0.920800 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(11), DoubleNear(-0.815800 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(24), DoubleNear(0.7065 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(25), DoubleNear(-0.7065 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(27), DoubleNear(0.815800 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(30), DoubleNear(18.3654 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(41), DoubleNear(17.2884 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(52), DoubleNear(-0.962700 / Utils::Constants::ev_per_hartree, 1e-4));
  ASSERT_THAT(integral->get(56), DoubleNear(0.962700 / Utils::Constants::ev_per_hartree, 1e-4));
}
} // namespace Sparrow
} // namespace Scine
