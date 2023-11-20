/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "BasisPruner.h"
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/TransitionChargesCalculator.h>
#include <Sparrow/Implementations/Exceptions.h>
#include <Sparrow/Implementations/TimeDependent/TimeDependentUtils.h>
#include <Utils/Constants.h>

namespace Scine {
namespace Sparrow {

template<Utils::Reference restrictedness>
BasisPruner<restrictedness>::BasisPruner(const OrderedInput<restrictedness>& input, std::shared_ptr<Eigen::MatrixXd> gammaMatrix,
                                         std::shared_ptr<Eigen::VectorXd> spinConstants)
  : gammaMatrix_(std::move(gammaMatrix)), spinConstants_(std::move(spinConstants)), input_(input) {
  assert(input_.energyDifferences().size() == static_cast<Eigen::Index>(input_.transitionCharges().rows()));
  assert(input_.energyDifferences().size() == static_cast<Eigen::Index>(input_.excitations().size()));
  isIncluded_ = BoolVector::Constant(input_.energyDifferences().size(), false);
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::prune(EnergyThreshold enThresh, PerturbativeThreshold ptThresh,
                                        Utils::SpinTransition spinBlock) -> OrderedInput<restrictedness> {
  check(spinBlock);

  generatePruningInformation(enThresh, ptThresh, spinBlock);

  return assembleResult(spinBlock);
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::prune(NumberOfConfigurations nRoots, PerturbativeThreshold ptThresh,
                                        Utils::SpinTransition spinBlock) -> OrderedInput<restrictedness> {
  check(spinBlock);

  generatePruningInformation(nRoots, ptThresh, spinBlock);

  return assembleResult(spinBlock);
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::assembleResult(Utils::SpinTransition /*spinBlock*/) -> OrderedInput<restrictedness> {
  nBasisFunctionsAfterPruning_ = std::count(isIncluded_.data(), isIncluded_.data() + isIncluded_.size(), true);

  OrderedInput<restrictedness> prunedResults;
  prunedResults.transitionCharges().resize(nBasisFunctionsAfterPruning_, input_.transitionCharges().cols());
  prunedResults.energyDifferences().resize(nBasisFunctionsAfterPruning_);
  prunedResults.excitations().reserve(nBasisFunctionsAfterPruning_);

  int index = 0;
  for (int basisFunctionIndex = 0; basisFunctionIndex < isIncluded_.size(); ++basisFunctionIndex) {
    if (isIncluded_(basisFunctionIndex)) {
      prunedResults.energyDifferences()(index) = input_.energyDifferences()(basisFunctionIndex);
      prunedResults.transitionCharges().row(index) = input_.transitionCharges().row(basisFunctionIndex);
      prunedResults.excitations().push_back(input_.excitations()[basisFunctionIndex]);
      ++index;
    }
  }
  conditionalFillData(prunedResults);

  return prunedResults;
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::generatePruningInformation(EnergyThreshold enThresh, PerturbativeThreshold ptThresh,
                                                             Utils::SpinTransition spinBlock) -> void {
  isIncluded_ = BoolVector::Constant(input_.energyDifferences().size(), false);

  auto* it = std::lower_bound(input_.energyDifferences().data(),
                              input_.energyDifferences().data() + input_.energyDifferences().size(), enThresh.threshold);

  nBasisFunctionsUnderThreshold_ = std::distance(input_.energyDifferences().data(), it);

  perturbativeCorrection(ptThresh, spinBlock);
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::generatePruningInformation(NumberOfConfigurations nRoots, PerturbativeThreshold ptThresh,
                                                             Utils::SpinTransition spinBlock) -> void {
  isIncluded_ = BoolVector::Constant(input_.energyDifferences().size(), false);

  nBasisFunctionsUnderThreshold_ = nRoots.number;

  perturbativeCorrection(ptThresh, spinBlock);
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::perturbativeCorrection(PerturbativeThreshold ptThresh, Utils::SpinTransition spinBlock)
    -> void {
  int numberOfSecondary = input_.energyDifferences().size() - nBasisFunctionsUnderThreshold_;

  if (nBasisFunctionsUnderThreshold_ == 0) {
    throw std::runtime_error("No configurations included in pruned space! Maybe energy threshold too high?");
  }

  isIncluded_.head(nBasisFunctionsUnderThreshold_).array() = true;

  if (numberOfSecondary != 0) {
    Eigen::VectorXd perturbationContributions = perturbationContributionVector(numberOfSecondary, spinBlock);

    assert(perturbationContributions.size() == numberOfSecondary);

    for (int secondaryIndex = 0; secondaryIndex < numberOfSecondary; ++secondaryIndex) {
      if (std::abs(perturbationContributions(secondaryIndex)) >= ptThresh.threshold) {
        isIncluded_(nBasisFunctionsUnderThreshold_ + secondaryIndex) = true;
      }
    }
  }
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::perturbationContributionVector(int nOfSecondary, Utils::SpinTransition spinBlock)
    -> Eigen::VectorXd {
  // Cols: transitions Rows: atoms
  // Loop needed in stead of colwise needed because of a known bug in the MKL
  // assignment in Eigen/3.2.2, Bug 1527 in changelog.
  // Original version:
  // enWtransCharges.array().colwise() *= input_.energyDifferences().array().sqrt();
  Eigen::MatrixXd enWtransCharges = input_.transitionCharges();
  for (int col = 0; col < int(enWtransCharges.cols()); ++col) {
    enWtransCharges.col(col).array() *= input_.energyDifferences().array().sqrt();
  }
  const Eigen::MatrixXd primaryCharges = enWtransCharges.topRows(nBasisFunctionsUnderThreshold_);
  const Eigen::MatrixXd secondaryCharges = enWtransCharges.bottomRows(nOfSecondary);

  return (generatePerturbationMatrix(primaryCharges, secondaryCharges, spinBlock).array().square() *
          generateEnergyWeightingMatrix(spinBlock).array())
      .colwise()
      .sum();
}

template<Utils::Reference restrictedness>
auto BasisPruner<restrictedness>::generateEnergyWeightingMatrix(Utils::SpinTransition /*spinBlock*/) {
  int nPrimary = std::count(isIncluded_.data(), isIncluded_.data() + isIncluded_.size(), true);
  int nSecondary = isIncluded_.size() - nPrimary;
  Eigen::MatrixXd energyWeighting = Eigen::MatrixXd::Zero(nPrimary, nSecondary);

  // Calculate the diagonal elements of the coupling matrix
  // I tried to also use the complete diagonal (energy squared + coupling matrix element)
  // but it made almost no difference, except for the very last roots, but also not
  // extremely much. I had an error on the last root of 40 with 2140 transitions of 0.09 eV instead of 0.04 eV.
  // All roots except last 3 were perfectly matching.
  Eigen::VectorXd diagonalElements = input_.energyDifferences().array().square();

  int indexIncluded = 0;
  int indexExcluded = 0;
  for (int i = 0; i < isIncluded_.size(); ++i) {
    if (isIncluded_(i)) {
      energyWeighting.array().row(indexIncluded++) += diagonalElements(i);
    }
    else {
      energyWeighting.array().col(indexExcluded++) -= diagonalElements(i);
    }
  }

  energyWeighting = Eigen::inverse(energyWeighting.array()).eval();
  return energyWeighting;
}

template<>
inline auto BasisPruner<Utils::Reference::Restricted>::generatePerturbationMatrix(const Eigen::MatrixXd& primaryCharges,
                                                                                  const Eigen::MatrixXd& secondaryCharges,
                                                                                  Utils::SpinTransition spinBlock)
    -> Eigen::MatrixXd {
  return 4.0 * ((spinBlock == Utils::SpinTransition::Singlet)
                    ? Eigen::MatrixXd(primaryCharges * gammaMatrix_->selfadjointView<Eigen::Lower>() * secondaryCharges.transpose())
                    : Eigen::MatrixXd(primaryCharges * spinConstants_->asDiagonal() * secondaryCharges.transpose()));
}

template<>
inline auto BasisPruner<Utils::Reference::Unrestricted>::generatePerturbationMatrix(const Eigen::MatrixXd& primaryCharges,
                                                                                    const Eigen::MatrixXd& secondaryCharges,
                                                                                    Utils::SpinTransition /*spinBlock*/)
    -> Eigen::MatrixXd {
  assert(input_.isBeta().size() == (primaryCharges.rows() + secondaryCharges.rows()));

  Eigen::MatrixXd magnetizationComponent = primaryCharges * spinConstants_->asDiagonal() * secondaryCharges.transpose();
  for (int i = 0; i < primaryCharges.rows(); ++i) {
    for (int j = 0; j < secondaryCharges.rows(); ++j) {
      if (input_.isBeta()(i) != input_.isBeta()(primaryCharges.rows() + j)) {
        magnetizationComponent(i, j) *= -1.0;
      }
    }
  }
  return 2.0 * (primaryCharges * gammaMatrix_->selfadjointView<Eigen::Lower>() * secondaryCharges.transpose() +
                magnetizationComponent);
}

template<Utils::Reference restrictedness>
inline auto BasisPruner<restrictedness>::generateDiagonalCouplings(const Eigen::VectorXd& energies) -> Eigen::VectorXd {
  Eigen::VectorXd result = energies.array().square();
  return result;
}

template<>
inline void
BasisPruner<Utils::Reference::Restricted>::conditionalFillData(OrderedInput<Utils::Reference::Restricted>& /*prunedData*/) const {
}

template<>
inline void
BasisPruner<Utils::Reference::Unrestricted>::conditionalFillData(OrderedInput<Utils::Reference::Unrestricted>& prunedData) const {
  prunedData.isBeta().resize(nBasisFunctionsAfterPruning_);
  int index = 0;
  for (int basisFunctionIndex = 0; basisFunctionIndex < isIncluded_.size(); ++basisFunctionIndex) {
    if (isIncluded_(basisFunctionIndex)) {
      prunedData.isBeta()(index) = input_.isBeta()(basisFunctionIndex);
      ++index;
    }
  }
}

template<>
inline auto BasisPruner<Utils::Reference::Restricted>::check(Utils::SpinTransition spinBlock) const -> void {
  if (spinBlock == Utils::SpinTransition::Triplet && !spinConstants_) {
    throw SpinConstantsNotAvailableException();
  }
}
template<>
inline auto BasisPruner<Utils::Reference::Unrestricted>::check(Utils::SpinTransition /*spinBlock*/) const -> void {
  if (!spinConstants_) {
    throw SpinConstantsNotAvailableException();
  }
}

template<>
auto BasisPruner<Utils::Reference::Restricted>::prune(const LinearResponseCalculator::GuessSpecifier& matrixToPrune) const
    -> std::shared_ptr<LinearResponseCalculator::GuessSpecifier> {
  LinearResponseCalculator::GuessSpecifier result;
  if (matrixToPrune.singlet.size() != 0) {
    result.singlet.resize(isIncluded_.size(), matrixToPrune.singlet.cols());
    int index = 0;
    for (int basisFunctionIndex = 0; basisFunctionIndex < isIncluded_.size(); ++basisFunctionIndex) {
      if (isIncluded_(basisFunctionIndex)) {
        result.singlet.row(index) = matrixToPrune.singlet.row(basisFunctionIndex);
        ++index;
      }
    }
  }
  if (matrixToPrune.triplet.size() != 0) {
    result.triplet.resize(isIncluded_.size(), matrixToPrune.triplet.cols());
    int index = 0;
    for (int basisFunctionIndex = 0; basisFunctionIndex < isIncluded_.size(); ++basisFunctionIndex) {
      if (isIncluded_(basisFunctionIndex)) {
        result.triplet.row(index) = matrixToPrune.triplet.row(basisFunctionIndex);
        ++index;
      }
    }
  }
  if (result.singlet.size() + result.triplet.size() == 0) {
    return std::shared_ptr<LinearResponseCalculator::GuessSpecifier>();
  }
  return std::make_shared<LinearResponseCalculator::GuessSpecifier>(result);
}

template<>
auto BasisPruner<Utils::Reference::Unrestricted>::prune(const LinearResponseCalculator::GuessSpecifier& matrixToPrune) const
    -> std::shared_ptr<LinearResponseCalculator::GuessSpecifier> {
  LinearResponseCalculator::GuessSpecifier result;
  if (matrixToPrune.unrestricted.size() != 0) {
    result.unrestricted.resize(isIncluded_.size(), matrixToPrune.unrestricted.cols());
    int index = 0;
    for (int basisFunctionIndex = 0; basisFunctionIndex < isIncluded_.size(); ++basisFunctionIndex) {
      if (isIncluded_(basisFunctionIndex)) {
        result.unrestricted.row(index) = matrixToPrune.unrestricted.row(basisFunctionIndex);
        ++index;
      }
    }
    return std::make_shared<LinearResponseCalculator::GuessSpecifier>(result);
  }
  return std::shared_ptr<LinearResponseCalculator::GuessSpecifier>();
}

template class BasisPruner<Utils::Reference::Restricted>;
template class BasisPruner<Utils::Reference::Unrestricted>;

} // namespace Sparrow
} // namespace Scine
