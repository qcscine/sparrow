/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_BASISPRUNER_H
#define SPARROW_BASISPRUNER_H
#include "OrderedInput.h"
#include <Sparrow/Implementations/TimeDependent/LinearResponseCalculator.h>
#include <Utils/Math/IterativeDiagonalizer/SigmaVectorEvaluator.h>
#include <Utils/Math/IterativeDiagonalizer/SpinAdaptedEigenContainer.h>
#include <vector>

namespace Scine {
namespace Utils {
enum class SpinTransition;
class DipoleMatrix;
class Excitation;
} // namespace Utils
namespace Sparrow {

/**
 * @brief Strongly named int representing the number of electronic configurations.
 * @struct NumberOfConfigurations
 * Used in return type to give some infos on the number.
 */
struct NumberOfConfigurations {
  int number;
};

/**
 * @brief Strongly named double representing an energy threshold.
 * @struct EnergyThreshold
 * Used in function signature to differentiate from a PerturbativeThreshold.
 */
struct EnergyThreshold {
  double threshold;
};
/**
 * @brief Strongly named double representing a perturbation threshold.
 * @struct PerturbativeThreshold
 * Used in function signature to differentiate from an EnergyThreshold.
 */
struct PerturbativeThreshold {
  double threshold;
};

/**
 * @brief This class takes care of pruning the singly excited determinant space. It is called when a pruned calculation
 * is started. It prunes the vector of energy differences of the substituted orbitals in a determinant and the
 * respective exciations (short: energy difference vector).
 *
 * The pruning can be done with three different criteria: an energy, TODO: an intensity and TODO:a combination of both.
 *
 * Pruning with the energy criterion (pruningWithEnergy) according to
 *
 * S. Grimme, A simplified Tamm-Dancoff density functional approach for the electronic
 * excitation spectra of very large molecules, J. Chem. Phys, 2013, 130, 244104.
 *
 * includes all the basis functions with a virtual-minus-occupied molecular orbital energy difference below
 * a user-defined threshold, as well as all the basis functions contributing much to this reduced space
 * as estimated by perturbation theory. The contribution of all the other basis functions is estimated by
 * perturbation theory and added to the diagonal of the matrix (added to the energy difference vector).
 *
 * TODO: Pruning with the intensity criterion.
 * Pruning with the intensity criterion(pruningWithIntensity) builds a vector of all excitations that have a
 * corresponding dipole matrix element ($f|r^2|$f) that is higher than the threshold.
 *
 * TODO: Pruning with the hybrid method.
 * The hybrid method combines both methods by calling both functions sequentially (intensity first).
 *
 *
 * @class Basispruner
 * @file BasisPruner.h
 * @tparam restrictedness
 */
template<Utils::Reference restrictedness>
class BasisPruner {
 public:
  using BoolVector = Eigen::Matrix<bool, -1, 1>;
  BasisPruner(const OrderedInput<restrictedness>& input, std::shared_ptr<Eigen::MatrixXd> gammaMatrix,
              std::shared_ptr<Eigen::VectorXd> spinConstants);

  /**
   * @brief Function making the pruning of the space.
   * @param enThresh Strongly typed double, represents the energy threshold (in hartree).
   * @param ptThresh Strongly typed double, represents the perturbation theory threshold
   *        for the inclusion of a secondary configuration. (in hartree)
   *
   *  All the configurations lower in energy than enThresh are included, called
   *  "primary configurations". The other configurations ("secondary configurations")
   *  are tested for interaction with the primary space via perturbation theory:
   *           E_u = \sum_v^PBF |A_uv|^2 / (E_u - E_v).
   *  if E_u, the estimated interaction energy of the secondary configuration u,
   *  is higher than the ptThresh, then it is included in the calculation.
   *  Otherwise, the estimated interaction is added to the energy differences.
   *
   * @return a OrderedInput<restrictedness>.
   *
   */
  auto prune(EnergyThreshold enThresh, PerturbativeThreshold ptThresh,
             Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet) -> OrderedInput<restrictedness>;

  /**
   * @brief Function making the pruning of the space.
   * @param nRoots Strongly typed int, represents the number of roots sought after.
   * @param ptThresh Strongly typed double, represents the perturbation theory threshold
   *        for the inclusion of a secondary configuration. (in hartree)
   *
   *  All the first nRoots configurations are included, called
   *  "primary configurations". The other configurations ("secondary configurations")
   *  are tested for interaction with the primary space via perturbation theory:
   *           E_u = \sum_v^PBF |A_uv|^2 / (E_u - E_v).
   *  if E_u, the estimated interaction energy of the secondary configuration u,
   *  is higher than the ptThresh, then it is included in the calculation.
   *  Otherwise, the estimated interaction is added to the energy differences.
   *
   * @return a OrderedInput<restrictedness>.
   */
  auto prune(NumberOfConfigurations nRoots, PerturbativeThreshold ptThresh,
             Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet) -> OrderedInput<restrictedness>;

  /**
   * @brief Returns the number of "primary" configurations under the energy threshold.
   */
  auto getNumberOfRootsUnderThreshold() const -> int {
    return nBasisFunctionsUnderThreshold_;
  }

  /**
   * @brief Prunes a matrix with the already calculated isIncluded_ private member.
   * The matrix will have less-equal rows.
   * @pre The energy pruning must already have taken place.
   * @pre matrixToPrune.rows() == isIncluded_.size()
   */
  auto prune(const LinearResponseCalculator::GuessSpecifier& matrixToPrune) const
      -> std::shared_ptr<LinearResponseCalculator::GuessSpecifier>;

 private:
  auto generatePruningInformation(EnergyThreshold enThresh, PerturbativeThreshold ptThresh, Utils::SpinTransition spinBlock)
      -> void;
  auto generatePruningInformation(NumberOfConfigurations nRoots, PerturbativeThreshold ptThresh,
                                  Utils::SpinTransition spinBlock) -> void;
  auto perturbativeCorrection(PerturbativeThreshold ptThresh, Utils::SpinTransition spinBlock) -> void;
  auto assembleResult(Utils::SpinTransition spinBlock) -> OrderedInput<restrictedness>;
  auto perturbationContributionVector(int nOfSecondary, Utils::SpinTransition spinBlock) -> Eigen::VectorXd;
  auto generateEnergyWeightingMatrix(Utils::SpinTransition spinBlock);

  /**
   * @brief Estimates with perturbation theory the importance of basis functions above the energy threshold.
   * The estimated contribution of basis function u is given by:
   *
   * E_u = \sum_v^PBF |A_uv|^2 / (E_u - E_v)
   *
   */
  auto generatePerturbationMatrix(const Eigen::MatrixXd& primaryCharges, const Eigen::MatrixXd& secondaryCharges,
                                  Utils::SpinTransition spinBlock) -> Eigen::MatrixXd;

  /**
   * @brief Conditionally fills the result with the "isBeta" bool vector.
   * Note: this function is inlined and does nothing in the restricted case.
   */
  auto conditionalFillData(OrderedInput<restrictedness>& prunedData) const -> void;

  auto generateDiagonalCouplings(const Eigen::VectorXd& energies) -> Eigen::VectorXd;

  auto check(Utils::SpinTransition spinBlock) const -> void;

  std::shared_ptr<Eigen::MatrixXd> gammaMatrix_;
  std::shared_ptr<Eigen::VectorXd> spinConstants_;
  OrderedInput<restrictedness> input_;
  int nBasisFunctionsAfterPruning_;
  int nBasisFunctionsUnderThreshold_;
  BoolVector isIncluded_;
};

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_BASISPRUNER_H
