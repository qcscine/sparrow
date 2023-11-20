/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_TRANSITIONCHARGESCALCULATOR_H
#define SPARROW_TRANSITIONCHARGESCALCULATOR_H

#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <map>
#include <vector>

namespace Scine {
namespace Utils {
class MolecularOrbitals;
class AtomsOrbitalsIndexes;
namespace LcaoUtils {
class ElectronicOccupation;
} // namespace LcaoUtils
} // namespace Utils
namespace Sparrow {

/// @brief Enum Class to store the transition charge matrices.
enum class AngularMomentum { None, S, P, D };
/// @brief Index composed by a <AtomIndex, AngularMomentum> pair.
using AtomAndOrbitalShellIndex = std::pair<int, AngularMomentum>;
template<Utils::Reference restrictedness>
using TransitionChargesContainer =
    Utils::SpinAdaptedContainer<restrictedness, std::map<AtomAndOrbitalShellIndex, Eigen::VectorXd>>;

class InvalidOccupationException : public std::exception {
  const char* what() const noexcept final {
    return "Occupation does not respect the Aufbau principle.";
  }
};

class SpinPolarizedTransitionChargesNotImplementedException : public std::exception {
  const char* what() const noexcept final {
    return "Spin polarized version non yet implemented.";
  }
};

class SpinPolarizedOrbitalsNotAvailableException : public std::exception {
  const char* what() const noexcept final {
    return "No unrestricted orbitals present in unrestricted calculation.";
  }
};
class TransitionChargesCalculator {
 public:
  TransitionChargesCalculator(const Utils::MolecularOrbitals& molecularOrbitals, const Eigen::MatrixXd& overlapMatrix,
                              const Utils::AtomsOrbitalsIndexes& aoIndex);

  /**
   * @brief Prepares the intermediate matrix S*c
   * This function must be called if a restricted calculation was run and
   * unrestricted quantities are needed. Or whenever the intermediate needs to be
   * reupdated (changes in molecularOrbitals or overlap matrix).
   */
  void fillOverlapProductMatrix();
  /**
   * @brief Calculates the nAtom matrices with element q^A_{ij} with i,j all the molecular orbitals.
   * q^A_ij is the partitioned transition charge from orbital i to j on the atom A.
   * @return A std::vector of size N_{Atoms}, containing for each element an MOxMO matrix.
   */
  std::vector<Eigen::MatrixXd> calculateMORestrictedAtomicChargeMatrices() const;
  std::vector<Eigen::MatrixXd> calculateMOUnrestrictedAtomicChargeMatrices() const;
  /**
   * @brief Calculates the nAtom matrices with element q^A_{ij} with i,j in the occ/vir block.
   * q^A_ij is the partitioned transition charge from orbital i to j on the atom A.
   * Assumes Aufbau construction.
   * @return A N_{transitions}xN_{Atoms} Eigen::Matrix containing the relative transition charge.
   * @tparam restrictedness Whether restricted or unrestricted charges needed.
   *
   * Implemented according to
   * A. Dominguez, B. Aradi, T. Frauenheim, V. Lutsker, T.A. Niehaus,
   * Extensions of the Time-Dependent Density Functional Based Tight-Binding Approach,
   * 2013
   *
   * q_{A}^\sigma = 1/2 * (c_\mu^\sigma * c'_\mu^{\sigma, T} + c'_\mu^\sigma * c_\mu^{\sigma, T})
   * A: Atomic index.
   * l: Orbital quantum number.
   * q_{A}^\sigma: Transition charge matrix for spin projection \sigma.
   * \mu: all atomic orbital indices centered on the atom A.
   * c_\mu^\sigma: Matrix with the molecular orbital coefficients
   *               for atomic orbitals \mu and spin projection \sigma.
   * c' = c*S.
   * S: overlap matrix.
   *
   */
  template<Utils::Reference restrictedness>
  Eigen::MatrixXd calculateAtomicTransitionChargeMatrices(const Utils::LcaoUtils::ElectronicOccupation& /*occupation*/) const {
    throw std::runtime_error("Wrong specialization in calculateAtomicTransitionChargeMatrices()");
  }
  /**
   * @brief Calculates the nAtomicOrbitals matrices with element q^\mu_{ij} with i,j all the molecular orbitals.
   * q^mu_ij is the partitioned transition charge from orbital i to j due to atomic orbital mu.
   * @return A N_{transitions}xN_{AOs} Eigen::Matrix containing the relative transition charge.
   */
  Eigen::MatrixXd calculateRestrictedTransitionChargeMatrices(const Utils::LcaoUtils::ElectronicOccupation& occupation) const;
  /**
   * @brief Calculates the nAtomicOrbitals matrices with element q^\mu_{ij} with i,j all the molecular orbitals.
   * q^mu_ij is the partitioned transition charge from orbital i to j due to atomic orbital mu.
   *
   * Implemented according to
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
   * @return A ((nOcc x nVir)^\alpha + (nOcc x nVir)^\beta) x nShellsTotal
   *         dimensional matrix (shell size is 1 per s orbital, 3 per p orbital, 5 for d,...),
   *         which contain as columns the linearized transition matrices containing the q_{A,l}^\sigma elements.
   */
  Eigen::MatrixXd calculateUnrestrictedTransitionChargeMatrices(const Utils::LcaoUtils::ElectronicOccupation& occupation) const;

 private:
  const Utils::MolecularOrbitals& mos_;
  const Eigen::MatrixXd& overlapMatrix_;
  const Utils::AtomsOrbitalsIndexes& aoIndex_;
  Utils::SpinAdaptedMatrix overlapProductMatrix_;
};

template<>
inline Eigen::MatrixXd TransitionChargesCalculator::calculateAtomicTransitionChargeMatrices<Utils::Reference::Restricted>(
    const Utils::LcaoUtils::ElectronicOccupation& occupation) const {
  if (!occupation.isFilledUpFromTheBottom()) {
    throw InvalidOccupationException();
  }

  auto nOccupied = occupation.numberOccupiedRestrictedOrbitals();
  auto nVirtual = mos_.numberOrbitals() - nOccupied;
  Eigen::MatrixXd transitionAtomicCharges = Eigen::MatrixXd::Zero(nOccupied * nVirtual, aoIndex_.getNAtoms());

  // Get correct coefficient matrix
  const Eigen::MatrixXd& occupiedBlock = mos_.restrictedMatrix().leftCols(nOccupied);
  const Eigen::MatrixXd& virtualBlock = mos_.restrictedMatrix().rightCols(nVirtual);
  const Eigen::MatrixXd& SCvir = overlapProductMatrix_.restrictedMatrix().rightCols(nVirtual);
  const Eigen::MatrixXd& SCocc = overlapProductMatrix_.restrictedMatrix().leftCols(nOccupied);
  for (int atom = 0; atom < aoIndex_.getNAtoms(); ++atom) {
    auto firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    auto nAOsOnAtom = aoIndex_.getNOrbitals(atom);
    Eigen::MatrixXd frontTransition =
        SCvir.middleRows(firstIndex, nAOsOnAtom).transpose() * occupiedBlock.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransition =
        virtualBlock.middleRows(firstIndex, nAOsOnAtom).transpose() * SCocc.middleRows(firstIndex, nAOsOnAtom);

    transitionAtomicCharges.col(atom) = 0.5 * (Eigen::Map<Eigen::VectorXd>(frontTransition.data(), frontTransition.size()) +
                                               Eigen::Map<Eigen::VectorXd>(backTransition.data(), backTransition.size()));
  }
  return transitionAtomicCharges;
}

template<>
inline Eigen::MatrixXd TransitionChargesCalculator::calculateAtomicTransitionChargeMatrices<Utils::Reference::Unrestricted>(
    const Utils::LcaoUtils::ElectronicOccupation& occupation) const {
  if (!occupation.isFilledUpFromTheBottom()) {
    throw InvalidOccupationException();
  }
  if (mos_.alphaMatrix().cols() == 0 || overlapProductMatrix_.alphaMatrix().cols() == 0) {
    throw SpinPolarizedOrbitalsNotAvailableException();
  }

  // Assumes Aufbau principle construction
  auto nOccupiedAlpha = occupation.numberAlphaElectrons();
  auto nVirtualAlpha = mos_.alphaMatrix().cols() - nOccupiedAlpha;
  auto nOccupiedBeta = occupation.numberBetaElectrons();
  auto nVirtualBeta = mos_.betaMatrix().cols() - nOccupiedBeta;
  Eigen::MatrixXd transitionChargeMatrix(nOccupiedAlpha * nVirtualAlpha + nOccupiedBeta * nVirtualBeta, aoIndex_.getNAtoms());
  const Eigen::MatrixXd& occupiedBlockAlpha = mos_.alphaMatrix().leftCols(nOccupiedAlpha);
  const Eigen::MatrixXd& virtualBlockAlpha = mos_.alphaMatrix().rightCols(nVirtualAlpha);
  const Eigen::MatrixXd& occupiedBlockBeta = mos_.betaMatrix().leftCols(nOccupiedBeta);
  const Eigen::MatrixXd& virtualBlockBeta = mos_.betaMatrix().rightCols(nVirtualBeta);
  const Eigen::MatrixXd& SCAlphavir = overlapProductMatrix_.alphaMatrix().rightCols(nVirtualAlpha);
  const Eigen::MatrixXd& SCAlphaocc = overlapProductMatrix_.alphaMatrix().leftCols(nOccupiedAlpha);
  const Eigen::MatrixXd& SCBetavir = overlapProductMatrix_.betaMatrix().rightCols(nVirtualBeta);
  const Eigen::MatrixXd& SCBetaocc = overlapProductMatrix_.betaMatrix().leftCols(nOccupiedBeta);

  for (int atom = 0; atom < aoIndex_.getNAtoms(); ++atom) {
    int firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    int nAOsOnAtom = aoIndex_.getNOrbitals(atom);
    Eigen::MatrixXd transitionChargesAlpha =
        0.5 *
        (SCAlphavir.middleRows(firstIndex, nAOsOnAtom).transpose() * occupiedBlockAlpha.middleRows(firstIndex, nAOsOnAtom) +
         virtualBlockAlpha.middleRows(firstIndex, nAOsOnAtom).transpose() * SCAlphaocc.middleRows(firstIndex, nAOsOnAtom));
    Eigen::MatrixXd transitionChargesBeta =
        0.5 *
        (SCBetavir.middleRows(firstIndex, nAOsOnAtom).transpose() * occupiedBlockBeta.middleRows(firstIndex, nAOsOnAtom) +
         virtualBlockBeta.middleRows(firstIndex, nAOsOnAtom).transpose() * SCBetaocc.middleRows(firstIndex, nAOsOnAtom));
    transitionChargeMatrix.col(atom)
        << Eigen::Map<Eigen::VectorXd>(transitionChargesAlpha.data(), transitionChargesAlpha.size()),
        Eigen::Map<Eigen::VectorXd>(transitionChargesBeta.data(), transitionChargesBeta.size());
  }
  return transitionChargeMatrix;
}

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_TRANSITIONCHARGESCALCULATOR_H
