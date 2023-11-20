/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CISSigmaVectorEvaluator.h"
#include "CISMatrixAOFockBuilderFactory.h"
#include <Utils/DataStructures/OccupiedMolecularOrbitals.h>
#include <iostream>

namespace Scine {
namespace Sparrow {
template<Utils::Reference restrictedness>
CISSigmaVectorEvaluator<restrictedness>::CISSigmaVectorEvaluator(
    CISData cisData, const ExcitedStatesParam& excitedStatesParam,
    const Utils::SpinAdaptedContainer<restrictedness, Eigen::VectorXd>& energyDifferenceVector,
    const std::vector<std::multimap<double, int, std::greater<double>>>& integralsThresholds, std::vector<int> orderMap,
    Utils::SpinTransition spinBlock)
  : cisData_(std::move(cisData)),
    energyDifferenceVector_(energyDifferenceVector),
    spinBlock_(spinBlock),
    integralsThresholds_(integralsThresholds),
    orderMap_(std::move(orderMap)) {
  currentSigmaMatrix_ = Eigen::MatrixXd(0, 0);
  aoFockBuilder_ = CISMatrixAOFockBuilderFactory<restrictedness>::createAOFockBuilder(spinBlock, cisData, excitedStatesParam);
  pseudoDensityBuilder_ =
      std::make_shared<CISPseudoDensityBuilder<restrictedness>>(cisData_.molecularOrbitals, cisData_.occupation);
  occupiedOrbitals_ = std::make_shared<Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>>(
      pseudoDensityBuilder_->getOccupiedOrbitals());
  virtualOrbitals_ = std::make_shared<Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>>(
      pseudoDensityBuilder_->getVirtualOrbitals());
}

namespace {
template<Utils::Reference restrictedness>
double getMaxElement(const Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>& /*pseudoDensityMatrix)*/) {
  return 0.0;
}

template<>
double getMaxElement<Utils::Reference::Restricted>(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd>& pseudoDensityMatrix) {
  return std::max(std::fabs(pseudoDensityMatrix.restricted.maxCoeff()), std::fabs(pseudoDensityMatrix.restricted.minCoeff()));
}

template<>
double getMaxElement<Utils::Reference::Unrestricted>(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd>& pseudoDensityMatrix) {
  return std::max(std::max(std::fabs(pseudoDensityMatrix.alpha.maxCoeff()), std::fabs(pseudoDensityMatrix.alpha.minCoeff())),
                  std::max(std::fabs(pseudoDensityMatrix.beta.maxCoeff()), std::fabs(pseudoDensityMatrix.beta.minCoeff())));
}
} // namespace

template<Utils::Reference restrictedness>
std::map<int, std::vector<int>> CISSigmaVectorEvaluator<restrictedness>::generateAtomPairList(
    const Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>& pseudoDensityMatrix) const {
  std::map<int, std::vector<int>> atomPairList;
  double maxDensity = getMaxElement<restrictedness>(pseudoDensityMatrix);
  for (int atomI = 0; atomI < cisData_.AOInfo.getNAtoms(); ++atomI) {
    int nAOsI = cisData_.AOInfo.getNOrbitals(atomI);
    std::vector<int> connectedAtoms;
    connectedAtoms.reserve(cisData_.AOInfo.getNAtoms() - atomI);
    for (const auto& integral : integralsThresholds_[atomI]) {
      int atomJ = integral.second;
      int nAOsJ = cisData_.AOInfo.getNOrbitals(atomJ);

      if (std::fabs(integral.first * maxDensity * (nAOsI * nAOsI + nAOsJ * nAOsJ + 2 * nAOsI * nAOsJ)) > 1e-14) {
        connectedAtoms.push_back(atomJ);
      }
      else {
        break;
      }
    }
    atomPairList.insert({atomI, std::move(connectedAtoms)});
    int totalSize = 0;
    for (const auto& atom : atomPairList) {
      totalSize += atom.second.size();
    }
  }
  return atomPairList;
}

template<>
inline const Eigen::MatrixXd&
CISSigmaVectorEvaluator<Utils::Reference::Restricted>::evaluate(const Eigen::MatrixXd& guessVectors) const {
  const int colsOldGuess = currentSigmaMatrix_.cols();
  const int dimCols = guessVectors.cols();
  const int dimRows = guessVectors.rows();
  const int colsNewSigmas = dimCols - colsOldGuess;
  const int nExcitations = energyDifferenceVector_.restricted.size();

  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd> newGuessVector;

  Eigen::MatrixXd sigmaMatrix;
  sigmaMatrix.conservativeResize(dimRows, colsNewSigmas);

  assert(dimRows == energyDifferenceVector_.restricted.rows());

#pragma omp parallel for schedule(dynamic) default(none) shared(sigmaMatrix, guessVectors) \
    firstprivate(colsOldGuess, dimCols, nExcitations) private(newGuessVector)
  for (int col = colsOldGuess; col < dimCols; col++) {
    TimeDependentUtils::transformOrder(guessVectors.col(col), newGuessVector.restricted, orderMap_,
                                       TimeDependentUtils::Direction::From);
    const auto pseudoDensityMatrix = pseudoDensityBuilder_->getPseudoDensityMatrix(newGuessVector);
    Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd> aoFockMatrix =
        aoFockBuilder_->getAOFock(pseudoDensityMatrix, generateAtomPairList(pseudoDensityMatrix));
    const Eigen::MatrixXd moFock =
        virtualOrbitals_->restricted.transpose() * aoFockMatrix.restricted * occupiedOrbitals_->restricted;
    Eigen::VectorXd sigmaVector(nExcitations);
    sigmaVector = Eigen::Map<const Eigen::VectorXd>(moFock.data(), moFock.size()) +
                  newGuessVector.restricted.cwiseProduct(energyDifferenceVector_.restricted);
    sigmaMatrix.col(col - colsOldGuess) = sigmaVector;
  }
  Eigen::MatrixXd reorderedSigmaMatrix;
  TimeDependentUtils::transformOrder(sigmaMatrix, reorderedSigmaMatrix, orderMap_, TimeDependentUtils::Direction::To);
  currentSigmaMatrix_.conservativeResize(dimRows, dimCols);
  assert(reorderedSigmaMatrix.cols() == colsNewSigmas && reorderedSigmaMatrix.rows() == dimRows);
  currentSigmaMatrix_.block(0, colsOldGuess, dimRows, colsNewSigmas) = reorderedSigmaMatrix;
  return currentSigmaMatrix_;
}

template<>
inline const Eigen::MatrixXd&
CISSigmaVectorEvaluator<Utils::Reference::Unrestricted>::evaluate(const Eigen::MatrixXd& guessVectors) const {
  const int colsOldGuess = currentSigmaMatrix_.cols();
  const int dimCols = guessVectors.cols();
  const int dimRows = guessVectors.rows();
  const int colsNewSigmas = dimCols - colsOldGuess;

  const int nExcitationsAlpha = energyDifferenceVector_.alpha.size();
  const int nExcitationsBeta = energyDifferenceVector_.beta.size();
  assert(nExcitationsAlpha + nExcitationsBeta == dimRows);

  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd> newGuessVector;
  Eigen::MatrixXd sigmaMatrix(dimRows, colsNewSigmas);

#pragma omp parallel for schedule(dynamic) default(none) shared(sigmaMatrix, guessVectors) \
    firstprivate(colsOldGuess, dimCols, nExcitationsAlpha, nExcitationsBeta) private(newGuessVector)
  for (int col = colsOldGuess; col < dimCols; col++) {
    Eigen::VectorXd reorderedVector;
    TimeDependentUtils::transformOrder(guessVectors.col(col), reorderedVector, orderMap_, TimeDependentUtils::Direction::From);

    newGuessVector.alpha = reorderedVector.head(nExcitationsAlpha);
    newGuessVector.beta = reorderedVector.tail(nExcitationsBeta);
    auto const pseudoDensityMatrix = pseudoDensityBuilder_->getPseudoDensityMatrix(newGuessVector);

    Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd> aoFockMatrix =
        aoFockBuilder_->getAOFock(pseudoDensityMatrix, generateAtomPairList(pseudoDensityMatrix));

    const Eigen::MatrixXd moAlphaFock = virtualOrbitals_->alpha.transpose() * aoFockMatrix.alpha * occupiedOrbitals_->alpha;
    const Eigen::MatrixXd moBetaFock = virtualOrbitals_->beta.transpose() * aoFockMatrix.beta * occupiedOrbitals_->beta;

    Eigen::VectorXd sigmaVectorAlpha(nExcitationsAlpha);
    Eigen::VectorXd sigmaVectorBeta(nExcitationsBeta);

    sigmaVectorAlpha = Eigen::Map<const Eigen::VectorXd>(moAlphaFock.data(), moAlphaFock.size());
    sigmaVectorBeta = Eigen::Map<const Eigen::VectorXd>(moBetaFock.data(), moBetaFock.size());

    sigmaMatrix.col(col - colsOldGuess)
        << reorderedVector.head(nExcitationsAlpha).cwiseProduct(energyDifferenceVector_.alpha) + sigmaVectorAlpha,
        reorderedVector.tail(nExcitationsBeta).cwiseProduct(energyDifferenceVector_.beta) + sigmaVectorBeta;
  }

  Eigen::MatrixXd reorderedSigmaMatrix;
  TimeDependentUtils::transformOrder(sigmaMatrix, reorderedSigmaMatrix, orderMap_, TimeDependentUtils::Direction::To);
  currentSigmaMatrix_.conservativeResize(dimRows, dimCols);
  currentSigmaMatrix_.block(0, colsOldGuess, dimRows, colsNewSigmas) = reorderedSigmaMatrix;
  return currentSigmaMatrix_;
}

template<Utils::Reference restrictedness>
void CISSigmaVectorEvaluator<restrictedness>::collapsed(int /*newSubspaceDimension*/) {
  currentSigmaMatrix_ = Eigen::MatrixXd(0, 0);
}

template<Utils::Reference restrictedness>
CISSigmaVectorEvaluator<restrictedness>::~CISSigmaVectorEvaluator() = default;
template class CISSigmaVectorEvaluator<Utils::Reference::Unrestricted>;
template class CISSigmaVectorEvaluator<Utils::Reference::Restricted>;
} // namespace Sparrow
} // namespace Scine
