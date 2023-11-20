/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParameterProcessor.h>
#include <Utils/Constants.h>
#include <Utils/Geometry/ElementInfo.h>
#include <gmock/gmock.h>

namespace Scine {
namespace Sparrow {

#define SQRT2 1.41421356237309504880

using namespace testing;
using namespace nddo;

class ARawParameterProcessing : public Test {
 public:
  ARawParameterProcessing() : processor(parameters) {
  }

  Parameters parameters;
  Parameters::Atomic iodine;
  Parameters::Diatomic cadmiumChlorine;
  std::unique_ptr<AtomicParameters> atomicPar;
  std::shared_ptr<PM6DiatomicParameters> diatomicPar;
  std::shared_ptr<OneCenterTwoElectronIntegrals> integral;
  RawParameterProcessor processor;

  void SetUp() override {
    // Those are the parameters for iodine
    iodine.pack.oneCenterEnergy = {-59.9732, -56.4598, -28.8226};
    iodine.pack.beta = {-30.5225, -5.94212, -7.67611};
    iodine.pack.gss = 7.23476;
    iodine.pack.gpp = 9.87747;
    iodine.pack.gsp = 9.15441;
    iodine.pack.gp2 = 8.03592;
    iodine.pack.hsp = 5.00422;
    iodine.pack.orbitalExponent = {4.49865, 1.91707, 1.87518};
    iodine.pack.internalExponent = {9.13524, 6.88819, 3.79152};
    iodine.pack.pcore = 1.8806;
    iodine.pack.f0sd = 20.495;
    iodine.pack.g2sd = 2.16014;
    // iodine.grepa = -0.035519;
    // iodine.grepb = 1.74439;
    // iodine.grepc = 1.22384;

    // Those are parameters for the Cd-Cl pair:
    cadmiumChlorine.exponent = 0.943547;
    cadmiumChlorine.factor = 0.140424;

    parameters.atomic.emplace(Utils::ElementInfo::Z(Utils::ElementType::I), iodine);
    parameters.diatomic.emplace(Parameters::key(Utils::ElementType::Cd, Utils::ElementType::Cl), cadmiumChlorine);
    auto par = processor.processAtomicParameters(Utils::ElementType::I);
    atomicPar = std::move(par.first);
    integral = std::move(par.second);
    diatomicPar = processor.runtimeDiatomicParameters(Utils::ElementType::Cd, Utils::ElementType::Cl);
  }
};

TEST_F(ARawParameterProcessing, SetsBetaValuesCorrectly) {
  ASSERT_THAT(atomicPar->betaS(), DoubleEq(iodine.pack.beta.s * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->betaP(), DoubleEq(iodine.pack.beta.p * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->betaD(), DoubleEq(iodine.pack.beta.d * Utils::Constants::hartree_per_ev));
}

TEST_F(ARawParameterProcessing, SetsUValuesCorrectly) {
  ASSERT_THAT(atomicPar->Uss(), DoubleEq(iodine.pack.oneCenterEnergy.s * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->Upp(), DoubleEq(iodine.pack.oneCenterEnergy.p * Utils::Constants::hartree_per_ev));
  ASSERT_THAT(atomicPar->Udd(), DoubleEq(iodine.pack.oneCenterEnergy.d * Utils::Constants::hartree_per_ev));
}

TEST_F(ARawParameterProcessing, SetsChargeSeparationsCorrectly) {
  auto d = atomicPar->chargeSeparations();

  ASSERT_THAT(d.get(multipole::MultipolePair::sp1), DoubleNear(0.3746525, 1e-4));
  ASSERT_THAT(d.get(multipole::MultipolePair::pp2), DoubleNear(1.89517161 / SQRT2, 1e-5));
  ASSERT_THAT(d.get(multipole::MultipolePair::sd2), DoubleNear(0.77747018 / SQRT2, 1e-5));
  ASSERT_THAT(d.get(multipole::MultipolePair::pd1), DoubleNear(1.29634171, 1e-5));
  ASSERT_THAT(d.get(multipole::MultipolePair::dd2), DoubleNear(1.63749938 / SQRT2, 1e-5));
}

TEST_F(ARawParameterProcessing, GtoExpansionHasBeenPerformed) {
  ASSERT_TRUE(atomicPar->GTOs().s);
  ASSERT_TRUE(atomicPar->GTOs().p);
  ASSERT_TRUE(atomicPar->GTOs().d);
}

/*TEST_F(ARawParameterProcessing, SetsGaussianRepulsionParameters) {
  auto g = atomicPar->getGaussianRepulsionParameters();

  ASSERT_TRUE(atomicPar->hasGaussianRepulsionParameters());
  ASSERT_THAT(std::get<0>(g), DoubleEq(iodine.grepa * Utils::Constants::bohr_per_angstrom));
  ASSERT_THAT(std::get<1>(g), DoubleEq(iodine.grepb * Utils::Constants::angstrom_per_bohr *
Utils::Constants::angstrom_per_bohr)); ASSERT_THAT(std::get<2>(g), DoubleEq(iodine.grepc *
Utils::Constants::bohr_per_angstrom));
}*/

TEST_F(ARawParameterProcessing, SetsDiatomicExponentCorrectly) {
  ASSERT_THAT(diatomicPar->alpha(), Eq(cadmiumChlorine.exponent * Utils::Constants::angstrom_per_bohr));
}

TEST_F(ARawParameterProcessing, SetsDiatomicFactorCorrectly) {
  ASSERT_THAT(diatomicPar->x(), Eq(cadmiumChlorine.factor));
}

TEST_F(ARawParameterProcessing, SetsPCoreCorrectly) {
  ASSERT_THAT(atomicPar->pCore(), DoubleEq(iodine.pack.pcore));
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
