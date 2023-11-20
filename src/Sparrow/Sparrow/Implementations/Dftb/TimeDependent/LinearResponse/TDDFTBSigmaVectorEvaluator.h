/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_TDDFTBSIGMAVECTOREVALUATOR_H
#define SPARROW_TDDFTBSIGMAVECTOREVALUATOR_H

#include <Sparrow/Implementations/Dftb/TimeDependent/LinearResponse/TDDFTBData.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/Math/IterativeDiagonalizer/SigmaVectorEvaluator.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>

namespace Scine {
namespace Sparrow {

template<Utils::Reference restrictedness>
class OrderedInput;

enum class TDDFTBType { TDDFTB, TDA };

/**
 * @class TDDFTBSigmaVectorEvaluator
 * @file TDDFTBigmaVectorEvaluator
 * @brief This class evaluates the sigma vector of the TDDFTB matrix for a KS reference.
 * The sigma matrix calculated previously is cached to be reused in subsequent iterations
 * of the Davidson-Liu algorithm.
 *
 * Implementation details for the singlet closed-shell case:
 * R. Rüger, E. van Lenthe, Y. Lu, J. Frenzel, T. Heine, L. Visscher,
 * Efficient Calculation of Electronic Absorption Spectra by Means of Intensity-Selected
 * Time-Dependent Density Functional Tight Binding, JCTC, 2015,
 * https://doi.org/10.1021/ct500838h
 *
 * Adaptations to Triplet and Unrestricted cases:
 * T. A. Niehaus, S. Suhai, F. Della Sala, P. Lugli, M. Elstner, G. Seifert, and Th. Frauenheim,
 * Tight-binding approach to time-dependent density-functional response theory,
 * Phys. Rev. B 63, 085108, 2001,
 * https://doi.org/10.1103/PhysRevB.63.085108
 *
 * Note on conditional inheritance:
 * std::conditional<SomeBool, A, B> gives A if SomeBool is true, B otherwise.
 * So, the class derives from the (empty) detail::RestrictedBase struct if restrictedness
 * is Restricted, and from the struct providing the isBeta_ member if restrictedness is
 * Unrestricted.
 */
template<Utils::Reference restrictedness>
class TDDFTBSigmaVectorEvaluator final : public Utils::SigmaVectorEvaluator {
 public:
  /**
   * @brief Constructor taking a OrderedInput<restrictedness> object.
   *
   * This object contains all the relevant matrices defined by the basis in the increasing energy order.
   */
  TDDFTBSigmaVectorEvaluator(std::shared_ptr<Eigen::MatrixXd> gammaMatrix, std::shared_ptr<Eigen::VectorXd> spinConstants,
                             const OrderedInput<restrictedness>& orderedInput,
                             Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet,
                             TDDFTBType type = TDDFTBType::TDDFTB);
  ~TDDFTBSigmaVectorEvaluator() final = default;

  /**
   * @brief Calculates one sigma vector per input guess vector.
   * @param guessVectors The guess vectors in the iterative calculation.
   * @pre TDDFTBData are not empty.
   * @pre Gamma matrix and Spin Constants dimensions match number of atoms.
   */
  auto evaluate(const Eigen::MatrixXd& guessVectors) const -> const Eigen::MatrixXd& final;

 private:
  /**
   * @brief precontract the transition charges with the square of the energy differences.
   */
  void calculateAtomicEnergyWeightedTransitionCharges(const Eigen::MatrixXd& transitionCharges);
  /**
   * @brief Contracts the guess vector with the energy-weighted transition charges to form matrix XBI.
   */
  auto calculateXBI(const Eigen::VectorXd& guessVector) const -> Eigen::MatrixXd;
  /**
   * @brief Contracts matrix XBI with the gamma/spin-coupling matrix to form matrix YAI.
   */
  auto calculateYAI(const Eigen::VectorXd& XBI) const -> Eigen::MatrixXd;
  /**
   * @brief Contracts matrix YAI with the energy-weighted transition charges.
   */
  auto calculateAtomicContraction(const Eigen::MatrixXd& YAIMatrix) const -> Eigen::VectorXd;
  /**
   * @brief In the unrestricted case, add the magnetization component to the final result.
   */
  template<typename Derived, typename OtherDerived>
  void fillAdditionalSigmaMatrixTerms(const Eigen::MatrixBase<Derived>& sigmaBlock,
                                      const Eigen::MatrixBase<OtherDerived>& guessVector) const;
  /**
   * @brief Resets the sigma matrix in the event of subspace collapse
   */
  void collapsed(int newSubspaceDimension) final;
  /**
   * @brief This gives 4 in the restricted case, 2 in the unrestricted case.
   */
  auto factor() const -> double;

  /**
   * @brief checks that the quantities needed for the caculation are present.
   * @throw If in unrestricted or triplet calculation no spinCouplings are present.
   */
  void check() const;
  const OrderedInput<restrictedness>& input_;
  std::shared_ptr<Eigen::MatrixXd> gammaMatrix_;
  std::shared_ptr<Eigen::VectorXd> spinConstantsVector_;
  const bool isTDA_;
  Utils::SpinTransition spinBlock_{Utils::SpinTransition::Singlet};
  Eigen::MatrixXd energyWeightedAtomicTransitionCharges_;
  // Caching variable declared mutable in order not to influence the API design.
  mutable Eigen::MatrixXd currentSigmaMatrix_;
};

template<>
inline auto TDDFTBSigmaVectorEvaluator<Utils::Reference::Restricted>::factor() const -> double {
  return isTDA_ ? 2.0 : 4.0;
}
template<>
inline auto TDDFTBSigmaVectorEvaluator<Utils::Reference::Unrestricted>::factor() const -> double {
  return isTDA_ ? 1.0 : 2.0;
}
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_TDDFTBSIGMAVECTOREVALUATOR_H
