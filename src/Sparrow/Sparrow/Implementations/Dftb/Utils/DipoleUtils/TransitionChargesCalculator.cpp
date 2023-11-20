/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "TransitionChargesCalculator.h"
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>

namespace Scine {
namespace Sparrow {

TransitionChargesCalculator::TransitionChargesCalculator(const Utils::MolecularOrbitals& molecularOrbitals,
                                                         const Eigen::MatrixXd& overlapMatrix,
                                                         const Utils::AtomsOrbitalsIndexes& aoIndex)
  : mos_(molecularOrbitals), overlapMatrix_(overlapMatrix), aoIndex_(aoIndex) {
  fillOverlapProductMatrix();
}

std::vector<Eigen::MatrixXd> TransitionChargesCalculator::calculateMORestrictedAtomicChargeMatrices() const {
  std::vector<Eigen::MatrixXd> atomicChargeMatrices(aoIndex_.getNAtoms());
  assert(mos_.numberOrbitals() != 0);
  assert(overlapProductMatrix_.restrictedMatrix().size() != 0);
  // Get correct coefficient matrix
  const Eigen::MatrixXd& coefMatrix = mos_.restrictedMatrix();
  for (int atom = 0; atom < aoIndex_.getNAtoms(); ++atom) {
    auto firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    auto nAOsOnAtom = aoIndex_.getNOrbitals(atom);
    Eigen::MatrixXd frontTransition =
        overlapProductMatrix_.restrictedMatrix().middleRows(firstIndex, nAOsOnAtom).transpose() *
        coefMatrix.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransition = coefMatrix.middleRows(firstIndex, nAOsOnAtom).transpose() *
                                     overlapProductMatrix_.restrictedMatrix().middleRows(firstIndex, nAOsOnAtom);

    atomicChargeMatrices[atom] = 0.5 * (frontTransition + backTransition);
  }
  return atomicChargeMatrices;
}

std::vector<Eigen::MatrixXd> TransitionChargesCalculator::calculateMOUnrestrictedAtomicChargeMatrices() const {
  std::vector<Eigen::MatrixXd> atomicChargeMatrices(aoIndex_.getNAtoms());
  assert(mos_.numberOrbitals() != 0);
  assert(overlapProductMatrix_.alphaMatrix().size() + overlapProductMatrix_.betaMatrix().size() != 0);
  // Get correct coefficient matrix
  auto alphaMatrix = mos_.alphaMatrix();
  auto betaMatrix = mos_.betaMatrix();
  auto const& SCAlpha = overlapProductMatrix_.alphaMatrix();
  auto const& SCBeta = overlapProductMatrix_.betaMatrix();
  for (int atom = 0; atom < static_cast<int>(atomicChargeMatrices.size()); ++atom) {
    auto firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    auto nAOsOnAtom = aoIndex_.getNOrbitals(atom);
    Eigen::MatrixXd frontTransitionAlpha =
        alphaMatrix.middleRows(firstIndex, nAOsOnAtom).transpose() * SCAlpha.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransitionAlpha =
        SCAlpha.middleRows(firstIndex, nAOsOnAtom).transpose() * alphaMatrix.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd frontTransitionBeta =
        betaMatrix.middleRows(firstIndex, nAOsOnAtom).transpose() * SCBeta.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransitionBeta =
        SCBeta.middleRows(firstIndex, nAOsOnAtom).transpose() * betaMatrix.middleRows(firstIndex, nAOsOnAtom);

    atomicChargeMatrices[atom] = 0.5 * (frontTransitionAlpha + backTransitionAlpha + frontTransitionBeta + backTransitionBeta);
  }
  return atomicChargeMatrices;
}

Eigen::MatrixXd TransitionChargesCalculator::calculateRestrictedTransitionChargeMatrices(
    const Utils::LcaoUtils::ElectronicOccupation& occupation) const {
  if (!occupation.isFilledUpFromTheBottom()) {
    throw InvalidOccupationException();
  }
  auto nOccupied = occupation.numberOccupiedRestrictedOrbitals();
  auto nVirtual = mos_.numberOrbitals() - nOccupied;
  Eigen::MatrixXd transitionAOCharges = Eigen::MatrixXd::Zero(nOccupied * nVirtual, mos_.restrictedMatrix().rows());

  // Get correct coefficient matrix
  const Eigen::MatrixXd& occupiedBlock = mos_.restrictedMatrix().leftCols(nOccupied);
  const Eigen::MatrixXd& virtualBlock = mos_.restrictedMatrix().rightCols(nVirtual);
  const Eigen::MatrixXd& SCvir = overlapProductMatrix_.restrictedMatrix().rightCols(nVirtual);
  const Eigen::MatrixXd& SCocc = overlapProductMatrix_.restrictedMatrix().leftCols(nOccupied);

  for (int AO = 0; AO < mos_.restrictedMatrix().rows(); ++AO) {
    Eigen::MatrixXd frontTransition = occupiedBlock.row(AO).transpose() * SCvir.row(AO);
    Eigen::MatrixXd backTransition = SCocc.transpose().col(AO) * virtualBlock.row(AO);

    transitionAOCharges.col(AO) = 0.5 * (Eigen::Map<Eigen::VectorXd>(frontTransition.data(), frontTransition.size()) +
                                         Eigen::Map<Eigen::VectorXd>(backTransition.data(), backTransition.size()));
  }
  return transitionAOCharges;
}

Eigen::MatrixXd TransitionChargesCalculator::calculateUnrestrictedTransitionChargeMatrices(
    const Utils::LcaoUtils::ElectronicOccupation& occupation) const {
  throw SpinPolarizedTransitionChargesNotImplementedException();
  if (!occupation.isFilledUpFromTheBottom()) {
    throw InvalidOccupationException();
  }
  // get number of shells
  int nShells = 0;
  for (int atom = 0; atom < aoIndex_.getNAtoms(); ++atom) {
    if (aoIndex_.getNOrbitals(atom) == 1)
      nShells += 1;
    else if (aoIndex_.getNOrbitals(atom) == 4)
      nShells += 2;
    else if (aoIndex_.getNOrbitals(atom) == 9)
      nShells += 3;
    else
      throw std::runtime_error("Only orbitals up to D supported.");
  }

  /*
   * Implement according to
   * A. Dominguez, B. Aradi, T. Frauenheim, V. Lutsker, T.A. Niehaus,
   * Extensions of the Time-Dependent Density Functional Based Tight-Binding Approach,
   * 2013
   *
   * q_{A,l}^\sigma = 1/2 * (c_\mu^\sigma * c'_\mu^{\sigma, T} + c'_\mu^\sigma * c_\mu^{\sigma, T})
   * A: Atomic index.
   * l: Orbital quantum number.
   * q_{A,l}^\sigma: Transition charge matrix for spin projection \sigma.
   * \mu: all atomic orbital indices in the set {A,l}.
   * c_\mu^\sigma: Matrix with the molecular orbital coefficients
   *               for atomic orbitals \mu and spin projection \sigma.
   * c' = c*S.
   * S: overlap matrix.
   *
   */
  // Assumes Aufbau principle construction
  auto nOccupiedAlpha = occupation.numberAlphaElectrons();
  auto nVirtualAlpha = mos_.alphaMatrix().cols() - nOccupiedAlpha;
  auto nOccupiedBeta = occupation.numberBetaElectrons();
  auto nVirtualBeta = mos_.betaMatrix().cols() - nOccupiedBeta;
  Eigen::MatrixXd transitionChargeMatrix(nOccupiedAlpha * nVirtualAlpha + nOccupiedBeta * nVirtualBeta, nShells);
  const Eigen::MatrixXd& occupiedBlockAlpha = mos_.alphaMatrix().leftCols(nOccupiedAlpha);
  const Eigen::MatrixXd& virtualBlockAlpha = mos_.alphaMatrix().rightCols(nVirtualAlpha);
  const Eigen::MatrixXd& occupiedBlockBeta = mos_.betaMatrix().leftCols(nOccupiedBeta);
  const Eigen::MatrixXd& virtualBlockBeta = mos_.betaMatrix().rightCols(nVirtualBeta);
  const Eigen::MatrixXd& SCAlphavir = overlapProductMatrix_.alphaMatrix().rightCols(nVirtualAlpha);
  const Eigen::MatrixXd& SCAlphaocc = overlapProductMatrix_.alphaMatrix().leftCols(nOccupiedAlpha);
  const Eigen::MatrixXd& SCBetavir = overlapProductMatrix_.betaMatrix().rightCols(nVirtualBeta);
  const Eigen::MatrixXd& SCBetaocc = overlapProductMatrix_.betaMatrix().leftCols(nOccupiedBeta);

  for (int atom = 0, shellIndex = 0; atom < aoIndex_.getNAtoms(); ++atom) {
    int nAOsOnA = aoIndex_.getNOrbitals(atom);
    int maxAngMomentum = (nAOsOnA == 1 ? 0 : (nAOsOnA == 4 ? 1 : 2));
    int firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    for (int angularMomentum = 0; angularMomentum <= maxAngMomentum; ++angularMomentum, ++shellIndex) {
      int nAOsInShell = angularMomentum * 2 + 1;
      Eigen::MatrixXd transitionChargesAlpha = 0.5 * (occupiedBlockAlpha.middleRows(firstIndex, nAOsInShell).transpose() *
                                                          SCAlphavir.middleRows(firstIndex, nAOsInShell) +
                                                      SCAlphaocc.middleRows(firstIndex, nAOsInShell).transpose() *
                                                          virtualBlockAlpha.middleRows(firstIndex, nAOsInShell));
      Eigen::MatrixXd transitionChargesBeta = 0.5 * (occupiedBlockBeta.middleRows(firstIndex, nAOsInShell).transpose() *
                                                         SCBetavir.middleRows(firstIndex, nAOsInShell) +
                                                     SCBetaocc.middleRows(firstIndex, nAOsInShell).transpose() *
                                                         virtualBlockBeta.middleRows(firstIndex, nAOsInShell));
      transitionChargeMatrix.col(shellIndex)
          << Eigen::Map<Eigen::VectorXd>(transitionChargesAlpha.data(), transitionChargesAlpha.size()),
          Eigen::Map<Eigen::VectorXd>(transitionChargesBeta.data(), transitionChargesBeta.size());
      firstIndex += nAOsInShell;
    }
  }
  return transitionChargeMatrix;
}

void TransitionChargesCalculator::fillOverlapProductMatrix() {
  if (mos_.isRestricted()) {
    overlapProductMatrix_.setRestrictedMatrix(overlapMatrix_ * mos_.restrictedMatrix());
  }
  else {
    overlapProductMatrix_.setAlphaMatrix(overlapMatrix_ * mos_.alphaMatrix());
    overlapProductMatrix_.setBetaMatrix(overlapMatrix_ * mos_.betaMatrix());
  }
}

} // namespace Sparrow
} // namespace Scine
