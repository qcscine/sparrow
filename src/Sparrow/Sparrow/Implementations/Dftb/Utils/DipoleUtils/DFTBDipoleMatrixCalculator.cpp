/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DFTBDipoleMatrixCalculator.h"
#include <Sparrow/Implementations/Dftb/Dftb0/DFTB0.h>
#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
#include <Sparrow/Implementations/Dftb/Dftb3/DFTB3.h>
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/TransitionChargesCalculator.h>

namespace Scine {
namespace Sparrow {

template<class DFTBMethod>
DFTBDipoleMatrixCalculator<DFTBMethod>::~DFTBDipoleMatrixCalculator() = default;

template<class DFTBMethod>
DFTBDipoleMatrixCalculator<DFTBMethod>::DFTBDipoleMatrixCalculator(DFTBMethod& method)
  : method_(method),
    positions_(method_.getPositions()),
    coefficientMatrix_(method_.getMolecularOrbitals()),
    overlapMatrix_(method_.getOverlapMatrix()),
    aoIndex_(method_.getAtomsOrbitalsIndexesHolder()) {
  invalidate();
}

template<class DFTBMethod>
std::unique_ptr<DFTBDipoleMatrixCalculator<DFTBMethod>> DFTBDipoleMatrixCalculator<DFTBMethod>::create(DFTBMethod& method) {
  DFTBDipoleMatrixCalculator<DFTBMethod> instance(method);
  return std::make_unique<DFTBDipoleMatrixCalculator<DFTBMethod>>(std::move(instance));
}

template<class DFTBMethod>
const Utils::DipoleMatrix& DFTBDipoleMatrixCalculator<DFTBMethod>::getAODipoleMatrix() const {
  throw DipoleMatrixTypeNotAvailableException();
}

template<class DFTBMethod>
Utils::DipoleMatrix DFTBDipoleMatrixCalculator<DFTBMethod>::getMODipoleMatrix() const {
  return dipoleMatrixMO_;
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::setAODipoleMatrix(Utils::DipoleMatrix /*dipoleMatrix*/) {
  throw DipoleMatrixTypeNotAvailableException();
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::setIntegralMethod(const IntegralMethod& /*IntegralMethod*/) {
  throw DipoleMatrixTypeNotAvailableException();
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::fillDipoleMatrix(const Eigen::RowVector3d& dipoleEvaluationCoordinate) {
  initialize();
  auto transitionChargesCalculator = TransitionChargesCalculator(coefficientMatrix_, overlapMatrix_, aoIndex_);
  // move origin
  auto positions = method_.getPositions();
  positions.rowwise() -= dipoleEvaluationCoordinate;

  std::vector<Eigen::MatrixXd> transitionChargeMatrices;
  if (coefficientMatrix_.isRestricted()) {
    transitionChargeMatrices = transitionChargesCalculator.calculateMORestrictedAtomicChargeMatrices();
  }
  else {
    transitionChargeMatrices = transitionChargesCalculator.calculateMOUnrestrictedAtomicChargeMatrices();
  }

  for (unsigned int atom = 0; atom < transitionChargeMatrices.size(); ++atom) {
    dipoleMatrixMO_.x().template get<Utils::DerivativeOrder::Zero>() +=
        transitionChargeMatrices[atom] * positions.row(atom).x();
    dipoleMatrixMO_.y().template get<Utils::DerivativeOrder::Zero>() +=
        transitionChargeMatrices[atom] * positions.row(atom).y();
    dipoleMatrixMO_.z().template get<Utils::DerivativeOrder::Zero>() +=
        transitionChargeMatrices[atom] * positions.row(atom).z();
  }
  valid_ = true;
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::initialize() {
  dipoleMatrixMO_.reset(coefficientMatrix_.numberOrbitals());
  valid_ = false;
  assert(dipoleMatrixMO_.x().template get<Utils::DerivativeOrder::Zero>().size() == overlapMatrix_.size());
}

template<class DFTBMethod>
void DFTBDipoleMatrixCalculator<DFTBMethod>::invalidate() {
  valid_ = false;
}

template<class DFTBMethod>
bool DFTBDipoleMatrixCalculator<DFTBMethod>::isValid() const {
  return valid_;
}

template class DFTBDipoleMatrixCalculator<dftb::DFTB0>;
template class DFTBDipoleMatrixCalculator<dftb::DFTB2>;
template class DFTBDipoleMatrixCalculator<dftb::DFTB3>;

} // namespace Sparrow
} // namespace Scine
