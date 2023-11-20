/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "RawParameterProcessor.h"
#include "AtomicParameters.h"
#include "PM6DiatomicParameters.h"
#include "PrincipalQuantumNumbers.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GeneralTypes.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/multipoleTypes.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ChargeSeparationParameter.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/KlopmanParameter.h>
#include <Utils/Constants.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/SlaterToGaussian.h>
#include <Utils/Geometry/ElementInfo.h>
#include <algorithm>

namespace Scine {
namespace Sparrow {

namespace nddo {
RawParameterProcessor::RawParameterProcessor(const Parameters& rawParameters, BasisFunctions basisFunctions)
  : rawParameters_(rawParameters), basisFunctions_(basisFunctions) {
}

std::unique_ptr<PM6DiatomicParameters> RawParameterProcessor::runtimeDiatomicParameters(Utils::ElementType e1,
                                                                                        Utils::ElementType e2) {
  const Parameters::Diatomic& p = rawParameters_.diatomic.at(Parameters::key(e1, e2));
  auto par = std::make_unique<PM6DiatomicParameters>(e1, e2);

  par->setX(p.factor);
  setDiatomicExponent(*par, e1, e2, p);

  return par;
}

std::unique_ptr<OneCenterTwoElectronIntegrals> RawParameterProcessor::get1c2eIntegrals(Utils::ElementType e,
                                                                                       const Parameters::Atomic& p) const {
  auto integrals = std::make_unique<OneCenterTwoElectronIntegrals>();

  integrals->setElement(e, basisFunctions_);
  integrals->setSlaterCondonParameters(&scParameters_);
  auto s = static_cast<int>(GeneralTypes::orb_t::s);
  auto p1 = static_cast<int>(GeneralTypes::orb_t::x);
  auto p2 = static_cast<int>(GeneralTypes::orb_t::y);
  if (p.pack.gss != 0)
    integrals->set(s, s, s, s, p.pack.gss * Utils::Constants::hartree_per_ev);
  if (p.pack.gsp != 0)
    integrals->set(s, s, p1, p1, p.pack.gsp * Utils::Constants::hartree_per_ev);
  if (p.pack.gpp != 0)
    integrals->set(p1, p1, p1, p1, p.pack.gpp * Utils::Constants::hartree_per_ev);
  if (p.pack.gp2 != 0)
    integrals->set(p1, p1, p2, p2, p.pack.gp2 * Utils::Constants::hartree_per_ev);
  if (p.pack.hsp != 0)
    integrals->set(s, p1, s, p1, p.pack.hsp * Utils::Constants::hartree_per_ev);

  integrals->calculateIntegrals();

  return integrals;
}

std::pair<std::unique_ptr<AtomicParameters>, std::unique_ptr<OneCenterTwoElectronIntegrals>>
RawParameterProcessor::processAtomicParameters(Utils::ElementType e) {
  const Parameters::Atomic& p = rawParameters_.atomic.at(Utils::ElementInfo::Z(e));
  auto runtimeAtomicPar = std::make_unique<AtomicParameters>(e, basisFunctions_);

  // If the element has d orbitals, the slater-condon parameters are needed
  if (runtimeAtomicPar->nAOs() == 9) {
    computeSlaterCondonParameters(*runtimeAtomicPar, p);
  }

  // TODO can refactor to single-arg setBeta, setOneCenterEnergy
  runtimeAtomicPar->setBetaS(p.pack.beta.s * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setBetaP(p.pack.beta.p * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setBetaD(p.pack.beta.d * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setUss(p.pack.oneCenterEnergy.s * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setUpp(p.pack.oneCenterEnergy.p * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setUdd(p.pack.oneCenterEnergy.d * Utils::Constants::hartree_per_ev);
  runtimeAtomicPar->setAlpha(p.pack.alpha * Utils::Constants::angstrom_per_bohr);
  setChargeSeparations(e, *runtimeAtomicPar, p);
  setKlopman(*runtimeAtomicPar, p);
  setGtoExpansion(e, *runtimeAtomicPar, p);

  for (const auto& gaussianRepulsion : p.gaussianRepulsion) {
    runtimeAtomicPar->addGaussianRepulsionParameters(gaussianRepulsion.a * Utils::Constants::bohr_per_angstrom,
                                                     gaussianRepulsion.b * Utils::Constants::angstrom_per_bohr *
                                                         Utils::Constants::angstrom_per_bohr,
                                                     gaussianRepulsion.c * Utils::Constants::bohr_per_angstrom);
  }

  auto atomicOneCenterIntegrals = get1c2eIntegrals(e, p);

  return std::make_pair(std::move(runtimeAtomicPar), std::move(atomicOneCenterIntegrals));
}

void RawParameterProcessor::setDiatomicExponent(PM6DiatomicParameters& par, Utils::ElementType e1,
                                                Utils::ElementType e2, const Parameters::Diatomic& p) {
  // For the pairs N-H and O-H, the given alpha has the unit [A^-2] instead of [A^-1]
  if ((e1 == Utils::ElementType::H &&
       (e2 == Utils::ElementType::O || e2 == Utils::ElementType::N || e2 == Utils::ElementType::C)) ||
      (e2 == Utils::ElementType::H &&
       (e1 == Utils::ElementType::O || e1 == Utils::ElementType::N || e1 == Utils::ElementType::C)))
    par.setAlpha(p.exponent * Utils::Constants::angstrom_per_bohr * Utils::Constants::angstrom_per_bohr);
  else
    par.setAlpha(p.exponent * Utils::Constants::angstrom_per_bohr);
}

// TODO can make static
void RawParameterProcessor::computeSlaterCondonParameters(AtomicParameters& runtimeAtomicPar, const Parameters::Atomic& p) {
  scParameters_ = SlaterCondonParameters();
  scParameters_.setElement(runtimeAtomicPar.element());
  // TODO can refactor to setInternalExponent
  scParameters_.setExponents(p.pack.internalExponent.s, p.pack.internalExponent.p, p.pack.internalExponent.d);
  if (p.pack.f0sd != 0)
    scParameters_.set(F0sd, p.pack.f0sd * Utils::Constants::hartree_per_ev);
  if (p.pack.g2sd != 0)
    scParameters_.set(G2sd, p.pack.g2sd * Utils::Constants::hartree_per_ev);
  scParameters_.calculate();
}

void RawParameterProcessor::setKlopman(AtomicParameters& par, const Parameters::Atomic& p) const {
  multipole::KlopmanParameter k;
  auto D = par.chargeSeparations();

  if (par.nAOs() == 1) {
    k.generateUpToS(p.pack.gss * Utils::Constants::hartree_per_ev);
  }
  else if (par.nAOs() == 4) {
    double hpp = 0.5 * (p.pack.gpp - p.pack.gp2) * Utils::Constants::hartree_per_ev;
#ifdef REPRODUCE_PM6_ERRORS
    hpp = std::max(0.1 * Utils::Constants::hartree_per_ev, hpp);
#endif
    k.generateUpToP(p.pack.gss * Utils::Constants::hartree_per_ev, p.pack.hsp * Utils::Constants::hartree_per_ev,
                    D.get(multipole::MultipolePair::sp1), hpp, D.get(multipole::MultipolePair::pp2));
  }
  else {
    double hpp = 0.5 * (p.pack.gpp - p.pack.gp2) * Utils::Constants::hartree_per_ev;
    //#ifdef REPRODUCE_PM6_ERRORS // Also for d elements?
    //    hpp = std::max(0.1 * Utils::Constants::hartree_per_ev, hpp);
    //#endif
    k.generateUpToD(p.pack.gss * Utils::Constants::hartree_per_ev, p.pack.hsp * Utils::Constants::hartree_per_ev,
                    D.get(multipole::MultipolePair::sp1), hpp, D.get(multipole::MultipolePair::pp2), scParameters_.get(F0dd),
                    scParameters_.get(G1pd), D.get(multipole::MultipolePair::pd1), scParameters_.get(G2sd),
                    D.get(multipole::MultipolePair::sd2), scParameters_.get(F2dd), D.get(multipole::MultipolePair::dd2));
  }
  par.setKlopmanParameters(k);
  if (p.pack.pcore == 0) {
    par.setPCore(k.get(multipole::MultipolePair::ss0));
  }
  else {
    par.setPCore(p.pack.pcore);
  }
  par.setPCoreSpecified(p.pack.pcore != 0);
}

void RawParameterProcessor::setChargeSeparations(Utils::ElementType e, AtomicParameters& par, const Parameters::Atomic& p) const {
  multipole::ChargeSeparationParameter d;

  unsigned int ns = PM6Elements::getQuantumNumberForSOrbital(e);
  unsigned int np = PM6Elements::getQuantumNumberForPOrbital(e);
  unsigned int nd = PM6Elements::getQuantumNumberForDOrbital(e);
  if (PM6Elements::getNumberOfAOs(e, basisFunctions_) == 4) {
    d.computeFromExponents(ns, np, p.pack.orbitalExponent.s, p.pack.orbitalExponent.p);
  }
  else if (PM6Elements::getNumberOfAOs(e, basisFunctions_) == 9) {
    d.computeFromExponents(ns, np, nd, p.pack.orbitalExponent.s, p.pack.orbitalExponent.p, p.pack.orbitalExponent.d);
  }
  par.setChargeSeparations(d);
}

void RawParameterProcessor::setGtoExpansion(Utils::ElementType e, AtomicParameters& par, const Parameters::Atomic& p) const {
  Utils::AtomicGtos gto;
  unsigned int N = 6; // STO-6G
  unsigned int ns = PM6Elements::getQuantumNumberForSOrbital(e);
  unsigned int np = PM6Elements::getQuantumNumberForPOrbital(e);
  unsigned int nd = PM6Elements::getQuantumNumberForDOrbital(e);
  gto.s = Utils::SlaterToGaussian::getGTOExpansion(N, ns, 0, p.pack.orbitalExponent.s);
  if (PM6Elements::getNumberOfAOs(e, basisFunctions_) >= 4) {
    gto.p = Utils::SlaterToGaussian::getGTOExpansion(N, np, 1, p.pack.orbitalExponent.p);
  }
  if (PM6Elements::getNumberOfAOs(e, basisFunctions_) == 9) {
    gto.d = Utils::SlaterToGaussian::getGTOExpansion(N, nd, 2, p.pack.orbitalExponent.d);
  }
  par.setGTOs(std::move(gto));
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
