/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_CISSIGMAVECTOREVALUATOR_H
#define SPARROW_CISSIGMAVECTOREVALUATOR_H

#include "CISData.h"
#include <Sparrow/Implementations/Nddo/TimeDependent/LinearResponse/CISPseudoDensityBuilder.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/DataStructures/OccupiedMolecularOrbitals.h>
#include <Utils/Math/IterativeDiagonalizer/SigmaVectorEvaluator.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

namespace Scine {
namespace Sparrow {
template<Utils::Reference restrictedness>
class CISMatrixAOFockBuilderBase;
/**
 * @brief This class evaluates the \sigma-vectors for the Davidson-Liu algorithm for both reference calculations (J. B.
 * Foresman, M. Head-Gordon, J. A. Pople and M. J. Frisch; J. Phys. Chem. 1992, 96, 1, 135-149.)
 *
 * This class is especially designed for the calculation of the sigma vectors in a pruned singly excited
 * determinant space. After the preceding calculation of the two-electron integrals in atomic orbital basis by the
 * CISAOFockBuilder, the calculation of the pseudo-Fock matrix in atomic orbital basis with the pseudo-densities from
 * the CISPseudoDensityBuilder and the energy-difference vector (TimeDependentUtils), the sigma-vectors can be
 * calculated like this
 *
 * for a restricted reference (S for singlet excitation, T for triplet excitation):
 *
 * \sigma^S_{(ia)} = b_(ia) * (e_a - e_i) + C^t ~ F^S_{ia}~ C
 * \sigma^T_{(ia)} = b_(ia) * (e_a - e_i) + F^T_{ia}
 *
 * and for an unrestricted reference (\alpha and \beta respectively)
 *
 * \sigma_{(ia)} = b_(ia) * (e_a - e_i) + F_{ia}
 *
 * \f$b\f$: guess vector
 *
 * \f$e\f$: energy level
 *
 * \f$F\f$: pseudo-Fock matrix
 *
 * @class CISTestSigmaVectorEvaluator
 * @file CISTestSigmaVectorEvaluator.h
 * @tparam restrictedness
 */
template<Utils::Reference restrictedness>
class CISSigmaVectorEvaluator final : public Utils::SigmaVectorEvaluator {
 public:
  /**
   * @brief Constructor: sets the cisData_ variable and calculates the energy difference-vector generation in the
   * TimeDependentUtils class.
   * @param cisData
   */
  explicit CISSigmaVectorEvaluator(CISData cisData, const ExcitedStatesParam& excitedStatesParam,
                                   const Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>& energyDifferenceVector,
                                   const std::vector<std::multimap<double, int, std::greater<double>>>& integralsThresholds,
                                   std::vector<int> orderMap,
                                   Utils::SpinTransition spinBlock = Utils::SpinTransition::Singlet);
  /**
   * @brief Destructor
   */
  ~CISSigmaVectorEvaluator() final;

  /**
   * @brief Evaluates the sigma vectors in a loop over all new guess vectors and inserts them in the currentSigmaMatrix_
   * variable (specified for both references).
   * @param guessVectors
   * @return
   */
  const Eigen::MatrixXd& evaluate(const Eigen::MatrixXd& guessVectors) const final;

 private:
  std::map<int, std::vector<int>>
  generateAtomPairList(const Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>& guessVector) const;
  /**
   * @brief Reset dimension after subspace collapse.
   */
  void collapsed(int newSubspaceDimension) final;
  mutable Eigen::MatrixXd currentSigmaMatrix_;
  CISData cisData_;
  std::shared_ptr<CISMatrixAOFockBuilderBase<restrictedness>> aoFockBuilder_;
  Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd> energyDifferenceVector_;
  Utils::SpinTransition spinBlock_{Utils::SpinTransition::Singlet};
  std::shared_ptr<CISPseudoDensityBuilder<restrictedness>> pseudoDensityBuilder_;
  std::shared_ptr<Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>> occupiedOrbitals_;
  std::shared_ptr<Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>> virtualOrbitals_;
  const std::vector<std::multimap<double, int, std::greater<double>>>& integralsThresholds_;
  std::vector<int> orderMap_;
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_CISSIGMAVECTOREVALUATOR_H
