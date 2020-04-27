/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RawParameterProcessor.h"
#include "AtomicParameters.h"
#include "PM6DiatomicParameters.h"
#include "PrincipalQuantumNumbers.h"
#include "RawParameters.h"
#include "RawParametersContainer.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GeneralTypes.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/multipoleTypes.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <Utils/Constants.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/SlaterToGaussian.h>
#include <algorithm>

namespace Scine {
namespace Sparrow {

namespace nddo {
RawParameterProcessor::RawParameterProcessor(const RawParametersContainer& rawParameters, BasisFunctions basisFunctions)
  : rawParameters_(rawParameters), basisFunctions_(basisFunctions) {
}

std::unique_ptr<PM6DiatomicParameters> RawParameterProcessor::runtimeDiatomicParameters(Utils::ElementType e1,
                                                                                        Utils::ElementType e2) {
  const RawDiatomicParameters& p = rawParameters_.getDiatomicParameters(e1, e2);
  auto par = std::make_unique<PM6DiatomicParameters>(e1, e2);

  par->setX(p.factor);
  setDiatomicExponent(*par, e1, e2, p);

  return par;
}

std::unique_ptr<OneCenterTwoElectronIntegrals> RawParameterProcessor::get1c2eIntegrals(Utils::ElementType e,
                                                                                       const RawAtomicParameters& p) const {
  auto integrals = std::make_unique<OneCenterTwoElectronIntegrals>();

  integrals->setElement(e, basisFunctions_);
  integrals->setSlaterCondonParameters(&scParameters_);
  auto s = static_cast<int>(GeneralTypes::orb_t::s);
  auto p1 = static_cast<int>(GeneralTypes::orb_t::x);
  auto p2 = static_cast<int>(GeneralTypes::orb_t::y);
  if (p.gss != 0)
    integrals->set(s, s, s, s, p.gss * Utils::Constants::hartree_per_ev);
  if (p.gsp != 0)
    integrals->set(s, s, p1, p1, p.gsp * Utils::Constants::hartree_per_ev);
  if (p.gpp != 0)
    integrals->set(p1, p1, p1, p1, p.gpp * Utils::Constants::hartree_per_ev);
  if (p.gp2 != 0)
    integrals->set(p1, p1, p2, p2, p.gp2 * Utils::Constants::hartree_per_ev);
  if (p.hsp != 0)
    integrals->set(s, p1, s, p1, p.hsp * Utils::Constants::hartree_per_ev);

  integrals->calculateIntegrals();

  return integrals;
}

std::pair<std::unique_ptr<AtomicParameters>, std::unique_ptr<OneCenterTwoElectronIntegrals>>
RawParameterProcessor::processAtomicParameters(Utils::ElementType e) {
  const RawAtomicParameters& p = rawParameters_.getAtomicParameters(e);
  auto runtimeAtomicPar = std::make_unique<AtomicParameters>(e, basisFunctions_);

  // If the element has d orbitals, the slater-condon parameters are needed
  if (runtimeAtomicPar->nAOs() == 9)
    computeSlaterCondonParameters(*runtimeAtomicPar, p);

  runtimeAtomicPar->setBetaS(p.bs * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setBetaP(p.bp * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setBetaD(p.bd * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setUss(p.uss * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setUpp(p.upp * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setUdd(p.udd * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setAlpha(p.alpha * Utils::Constants::angstrom_per_bohr);
  setChargeSeparations(e, *runtimeAtomicPar, p);
  setKlopman(*runtimeAtomicPar, p);
  setGtoExpansion(e, *runtimeAtomicPar, p);

  if (p.gaussianRepulsionParameters.size() > 0) {
    for (int i = 0; i < p.gaussianRepulsionParameters.size(); i += 1) {
      runtimeAtomicPar->addGaussianRepulsionParameters(
          p.gaussianRepulsionParameters[i].grepa * Utils::Constants::bohr_per_angstrom,
          p.gaussianRepulsionParameters[i].grepb * Utils::Constants::angstrom_per_bohr * Utils::Constants::angstrom_per_bohr,
          p.gaussianRepulsionParameters[i].grepc * Utils::Constants::bohr_per_angstrom);
    };
  };
  /*if (p.grepa != 0)
    runtimeAtomicPar->setGaussianRepulsionParameters(p.grepa * Utils::bohr_per_angstrom,
                                                     p.grepb * Utils::angstrom_per_bohr * Utils::angstrom_per_bohr,
                                                     p.grepc * Utils::bohr_per_angstrom);*/

  auto atomicOneCenterIntegrals = get1c2eIntegrals(e, p);

  return std::make_pair(std::move(runtimeAtomicPar), std::move(atomicOneCenterIntegrals));
}

void RawParameterProcessor::setDiatomicExponent(PM6DiatomicParameters& par, Utils::ElementType e1,
                                                Utils::ElementType e2, const RawDiatomicParameters& p) {
  // For the pairs N-H and O-H, the given alpha has the unit [A^-2] instead of [A^-1]
  if ((e1 == Utils::ElementType::H &&
       (e2 == Utils::ElementType::O || e2 == Utils::ElementType::N || e2 == Utils::ElementType::C)) ||
      (e2 == Utils::ElementType::H &&
       (e1 == Utils::ElementType::O || e1 == Utils::ElementType::N || e1 == Utils::ElementType::C)))
    par.setAlpha(p.exponent * Utils::Constants::angstrom_per_bohr * Utils::Constants::angstrom_per_bohr);
  else
    par.setAlpha(p.exponent * Utils::Constants::angstrom_per_bohr);
}

void RawParameterProcessor::computeSlaterCondonParameters(AtomicParameters& runtimeAtomicPar, const RawAtomicParameters& p) {
  scParameters_ = SlaterCondonParameters();
  scParameters_.setElement(runtimeAtomicPar.element());
  scParameters_.setExponents(p.zsn, p.zpn, p.zdn);
  if (p.f0sd != 0)
    scParameters_.set(F0sd, p.f0sd * Utils::Constants::hartree_per_ev);
  if (p.g2sd != 0)
    scParameters_.set(G2sd, p.g2sd * Utils::Constants::hartree_per_ev);
  scParameters_.calculate();
}

void RawParameterProcessor::setKlopman(AtomicParameters& par, const RawAtomicParameters& p) const {
  multipole::KlopmanParameter k;
  auto D = par.chargeSeparations();

  if (par.nAOs() == 1) {
    k.generateUpToS(p.gss * Utils::Constants::hartree_per_ev);
  }
  else if (par.nAOs() == 4) {
    double hpp = 0.5 * (p.gpp - p.gp2) * Utils::Constants::hartree_per_ev;
#ifdef REPRODUCE_PM6_ERRORS
    hpp = std::max(0.1 * Utils::Constants::hartree_per_ev, hpp);
#endif
    k.generateUpToP(p.gss * Utils::Constants::hartree_per_ev, p.hsp * Utils::Constants::hartree_per_ev,
                    D.get(multipole::sp1), hpp, D.get(multipole::pp2));
  }
  else {
    double hpp = 0.5 * (p.gpp - p.gp2) * Utils::Constants::hartree_per_ev;
    //#ifdef REPRODUCE_PM6_ERRORS // Also for d elements?
    //    hpp = std::max(0.1 * Utils::Constants::hartree_per_ev, hpp);
    //#endif
    k.generateUpToD(p.gss * Utils::Constants::hartree_per_ev, p.hsp * Utils::Constants::hartree_per_ev, D.get(multipole::sp1),
                    hpp, D.get(multipole::pp2), scParameters_.get(F0dd), scParameters_.get(G1pd), D.get(multipole::pd1),
                    scParameters_.get(G2sd), D.get(multipole::sd2), scParameters_.get(F2dd), D.get(multipole::dd2));
  }
  par.setKlopmanParameters(k);
  if (p.pcore == 0)
    par.setPCore(k.get(multipole::ss0));
  else
    par.setPCore(p.pcore);
  par.setPCoreSpecified(p.pcore != 0);
}

void RawParameterProcessor::setChargeSeparations(Utils::ElementType e, AtomicParameters& par, const RawAtomicParameters& p) const {
  multipole::ChargeSeparationParameter d;

  unsigned int ns = PM6Elements::getQuantumNumberForSOrbital(e);
  unsigned int np = PM6Elements::getQuantumNumberForPOrbital(e);
  unsigned int nd = PM6Elements::getQuantumNumberForDOrbital(e);
  if (PM6Elements::getNumberOfAOs(e, basisFunctions_) == 4)
    d.computeFromExponents(ns, np, p.zs, p.zp);
  else if (PM6Elements::getNumberOfAOs(e, basisFunctions_) == 9)
    d.computeFromExponents(ns, np, nd, p.zs, p.zp, p.zd);
  par.setChargeSeparations(d);
}

void RawParameterProcessor::setGtoExpansion(Utils::ElementType e, AtomicParameters& par, const RawAtomicParameters& p) const {
  Utils::AtomicGtos gto;
  unsigned int N = 6; // STO-6G
  unsigned int ns = PM6Elements::getQuantumNumberForSOrbital(e);
  unsigned int np = PM6Elements::getQuantumNumberForPOrbital(e);
  unsigned int nd = PM6Elements::getQuantumNumberForDOrbital(e);
  auto gtoS = Utils::SlaterToGaussian::getGTOExpansion(N, ns, 0, p.zs);
  gto.setS(gtoS);
  if (PM6Elements::getNumberOfAOs(e, basisFunctions_) >= 4) {
    auto gtoP = Utils::SlaterToGaussian::getGTOExpansion(N, np, 1, p.zp);
    gto.setP(gtoP);
  }
  if (PM6Elements::getNumberOfAOs(e, basisFunctions_) == 9) {
    auto gtoD = Utils::SlaterToGaussian::getGTOExpansion(N, nd, 2, p.zd);
    gto.setD(gtoD);
  }
  par.setGTOs(gto);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
