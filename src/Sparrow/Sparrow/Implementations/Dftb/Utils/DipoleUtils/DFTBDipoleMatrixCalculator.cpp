/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DFTBDipoleMatrixCalculator.h"
#include <Sparrow/Implementations/Dftb/Dftb0/DFTB0.h>
#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
#include <Sparrow/Implementations/Dftb/Dftb3/DFTB3.h>

namespace Scine {
namespace Sparrow {

template<class DFTBMethod>
DFTBDipoleMatrixCalculator<DFTBMethod>::DFTBDipoleMatrixCalculator(DFTBMethod& method)
  : method_(method),
    positions_(method_.getPositions()),
    coefficientMatrix_(method_.getMolecularOrbitals()),
    overlapMatrix_(method_.getOverlapMatrix()),
    aoIndex_(method_.getAtomsOrbitalsIndexesHolder()) {
  initialize();
}

template<class DFTBMethod>
std::unique_ptr<DFTBDipoleMatrixCalculator<DFTBMethod>> DFTBDipoleMatrixCalculator<DFTBMethod>::create(DFTBMethod& method) {
  DFTBDipoleMatrixCalculator<DFTBMethod> instance(method);
  return std::make_unique<DFTBDipoleMatrixCalculator<DFTBMethod>>(std::move(instance));
}

template<class DFTBMethod>
DFTBDipoleMatrixCalculator<DFTBMethod>::~DFTBDipoleMatrixCalculator() = default;

template<class DFTBMethod>
const Utils::DipoleMatrix& DFTBDipoleMatrixCalculator<DFTBMethod>::getAODipoleMatrix() const {
  throw DipoleMatrixTypeNotAvailableException();
}

template<class DFTBMethod>
Utils::DipoleMatrix DFTBDipoleMatrixCalculator<DFTBMethod>::getMODipoleMatrix() const {
  return dipoleMatrixMO_;
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::setAODipoleMatrix(Utils::DipoleMatrix dipoleMatrix) {
  throw DipoleMatrixTypeNotAvailableException();
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::setIntegralMethod(const IntegralMethod& IntegralMethod) {
  throw DipoleMatrixTypeNotAvailableException();
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::fillDipoleMatrix(const Eigen::RowVector3d& dipoleEvaluationCoordinate) {
  // move origin
  auto positions = positions_;
  positions.rowwise() -= dipoleEvaluationCoordinate;
  if (coefficientMatrix_.isRestricted())
    calculateRestrictedTransitionChargeMatrices(positions);
  else
    calculateUnrestrictedTransitionChargeMatrices(positions);

  for (int atom = 0; atom < transitionChargeMatrices_.size(); ++atom) {
    dipoleMatrixMO_.x().template get<Utils::derivOrder::zero>() += transitionChargeMatrices_[atom] * positions.row(atom).x();
    dipoleMatrixMO_.y().template get<Utils::derivOrder::zero>() += transitionChargeMatrices_[atom] * positions.row(atom).y();
    dipoleMatrixMO_.z().template get<Utils::derivOrder::zero>() += transitionChargeMatrices_[atom] * positions.row(atom).z();
  }
  valid_ = true;
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::initialize() {
  auto nMOs = coefficientMatrix_.isRestricted() ? coefficientMatrix_.restrictedMatrix().cols()
                                                : coefficientMatrix_.alphaMatrix().cols();
  dipoleMatrixMO_.reset(nMOs);
  transitionChargeMatrices_.resize(positions_.rows());
  for (auto& transitionCharge : transitionChargeMatrices_)
    transitionCharge.resize(nMOs, nMOs);
  valid_ = false;
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::invalidate() {
  valid_ = false;
}

template<class DFTBMethod>
bool DFTBDipoleMatrixCalculator<DFTBMethod>::isValid() const {
  return valid_;
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::calculateRestrictedTransitionChargeMatrices(const Utils::PositionCollection& positions) {
  // Get correct coefficient matrix
  auto coefMatrix = coefficientMatrix_.restrictedMatrix();
  // Get the two partial matrices C^T*S and S*C
  Eigen::MatrixXd CTS = coefMatrix.transpose() * overlapMatrix_.selfadjointView<Eigen::Upper>();
  Eigen::MatrixXd SC = overlapMatrix_.selfadjointView<Eigen::Upper>() * coefMatrix;

  for (int atom = 0; atom < transitionChargeMatrices_.size(); ++atom) {
    auto firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    auto nAOsOnAtom = aoIndex_.getNOrbitals(atom);
    Eigen::MatrixXd frontTransition =
        coefMatrix.transpose().middleCols(firstIndex, nAOsOnAtom) * SC.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransition = CTS.middleCols(firstIndex, nAOsOnAtom) * coefMatrix.middleRows(firstIndex, nAOsOnAtom);

    transitionChargeMatrices_[atom] = 0.5 * (frontTransition + backTransition);
  }
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::calculateUnrestrictedTransitionChargeMatrices(const Utils::PositionCollection& positions) {
  // Get correct coefficient matrix
  auto alphaMatrix = coefficientMatrix_.alphaMatrix();
  auto betaMatrix = coefficientMatrix_.betaMatrix();
  // Get the two partial matrices C^T*S and S*C
  Eigen::MatrixXd CTSalpha = alphaMatrix.transpose() * overlapMatrix_.selfadjointView<Eigen::Upper>();
  Eigen::MatrixXd SCalpha = overlapMatrix_.selfadjointView<Eigen::Upper>() * alphaMatrix;
  Eigen::MatrixXd CTSbeta = betaMatrix.transpose() * overlapMatrix_.selfadjointView<Eigen::Upper>();
  Eigen::MatrixXd SCbeta = overlapMatrix_.selfadjointView<Eigen::Upper>() * betaMatrix;

  for (int atom = 0; atom < transitionChargeMatrices_.size(); ++atom) {
    auto firstIndex = aoIndex_.getFirstOrbitalIndex(atom);
    auto nAOsOnAtom = aoIndex_.getNOrbitals(atom);
    Eigen::MatrixXd frontTransitionAlpha =
        alphaMatrix.transpose().middleCols(firstIndex, nAOsOnAtom) * SCalpha.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransitionAlpha =
        CTSalpha.middleCols(firstIndex, nAOsOnAtom) * alphaMatrix.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd frontTransitionBeta =
        betaMatrix.transpose().middleCols(firstIndex, nAOsOnAtom) * SCbeta.middleRows(firstIndex, nAOsOnAtom);
    Eigen::MatrixXd backTransitionBeta =
        CTSbeta.middleCols(firstIndex, nAOsOnAtom) * betaMatrix.middleRows(firstIndex, nAOsOnAtom);

    transitionChargeMatrices_[atom] =
        0.5 * (frontTransitionAlpha + backTransitionAlpha + frontTransitionBeta + backTransitionBeta);
  }
}

template class DFTBDipoleMatrixCalculator<dftb::DFTB0>;
template class DFTBDipoleMatrixCalculator<dftb::DFTB2>;
template class DFTBDipoleMatrixCalculator<dftb::DFTB3>;

} // namespace Sparrow
} // namespace Scine
