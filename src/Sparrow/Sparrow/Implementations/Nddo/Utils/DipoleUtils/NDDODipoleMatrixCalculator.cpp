/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NDDODipoleMatrixCalculator.h"
#include "Sparrow/Implementations/Nddo/Utils/DipoleUtils/AtomPairDipole.h"
#include "Sparrow/Implementations/Nddo/Utils/DipoleUtils/GTODipoleMatrixBlock.h"
#include <Sparrow/Implementations/Nddo/Am1/AM1Method.h>
#include <Sparrow/Implementations/Nddo/Mndo/MNDOMethod.h>
#include <Sparrow/Implementations/Nddo/Pm6/PM6Method.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOInitializer.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Geometry.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

template<class NDDOMethod>
NDDODipoleMatrixCalculator<NDDOMethod>::~NDDODipoleMatrixCalculator() = default;

template<class NDDOMethod>
void NDDODipoleMatrixCalculator<NDDOMethod>::initialize() {
  nAtoms_ = aoIndexes_.getNAtoms();
  nAOs_ = aoIndexes_.getNAtomicOrbitals();
  dipoleMatrix_.reset(nAOs_);
  valid_ = false;
}
template<class NDDOMethod>
void NDDODipoleMatrixCalculator<NDDOMethod>::setIntegralMethod(const IntegralMethod& method) {
  integralMethod_ = method;
  valid_ = false;
}

template<class NDDOMethod>
void NDDODipoleMatrixCalculator<NDDOMethod>::fillDipoleMatrix(const Eigen::RowVector3d& dipoleEvaluationCoordinate) {
  initialize();
  valid_ = false;
  for (int atom_i = 0; atom_i < nAtoms_; ++atom_i) {
    auto const firstAOIndex_i = aoIndexes_.getFirstOrbitalIndex(atom_i);
    auto const& parA = elementParameters_.get(elementTypes_[atom_i]);
    auto const Ri = positions_.row(atom_i);

    for (int atom_j = atom_i; atom_j < nAtoms_; ++atom_j) {
      auto const firstAOIndex_j = aoIndexes_.getFirstOrbitalIndex(atom_j);
      auto const& parB = elementParameters_.get(elementTypes_[atom_j]);
      auto const Rj = positions_.row(atom_j);

      Eigen::RowVector3d Rij = Rj - Ri;

      AtomPairDipole::fillAtomPairDipoleBlock(dipoleMatrix_, firstAOIndex_i, firstAOIndex_j, integralMethod_,
                                              parA.GTOs(), parB.GTOs(), Ri, Rj, Rij, dipoleEvaluationCoordinate);
    }
  }
  valid_ = true;
}

template<class NDDOMethod>
const Utils::DipoleMatrix& NDDODipoleMatrixCalculator<NDDOMethod>::getAODipoleMatrix() const {
  return dipoleMatrix_;
}

template<class NDDOMethod>
void NDDODipoleMatrixCalculator<NDDOMethod>::setAODipoleMatrix(Utils::DipoleMatrix dipoleMatrix) {
  dipoleMatrix_ = std::move(dipoleMatrix);
  valid_ = true;
}

template<class NDDOMethod>
bool NDDODipoleMatrixCalculator<NDDOMethod>::isValid() const {
  return valid_;
}

template<class NDDOMethod>
NDDODipoleMatrixCalculator<NDDOMethod>::NDDODipoleMatrixCalculator(NDDOMethod& method)
  : aoIndexes_(method.getAtomsOrbitalsIndexesHolder()),
    elementTypes_(method.getElementTypes()),
    positions_(method.getPositions()),
    elementParameters_(method.getInitializer().getElementParameters()),
    overlapMatrix_(method.getOverlapMatrix()),
    molecularOrbitals_(method.getMolecularOrbitals()) {
  initialize();
}

template<class NDDOMethod>
std::unique_ptr<NDDODipoleMatrixCalculator<NDDOMethod>> NDDODipoleMatrixCalculator<NDDOMethod>::create(NDDOMethod& method) {
  NDDODipoleMatrixCalculator<NDDOMethod> instance(method);
  return std::make_unique<NDDODipoleMatrixCalculator<NDDOMethod>>(std::move(instance));
}

template<class NDDOMethod>
Utils::DipoleMatrix NDDODipoleMatrixCalculator<NDDOMethod>::getMODipoleMatrix() const {
  if (molecularOrbitals_.isRestricted())
    return calculateMODipoleMatrixRestricted();
  else
    return calculateMODipoleMatrixUnrestricted();
}

template<class NDDOMethod>
void NDDODipoleMatrixCalculator<NDDOMethod>::invalidate() {
  valid_ = false;
}

template<class NDDOMethod>
Utils::DipoleMatrix NDDODipoleMatrixCalculator<NDDOMethod>::calculateMODipoleMatrixRestricted() const {
  Utils::DipoleMatrix moDipoleMatrix;
  auto& MOMatrix = molecularOrbitals_.restrictedMatrix();
  moDipoleMatrix.reset(MOMatrix.cols());
  // Get the MO matrix as D_{MO} = C^T * D_{AO} * C
  for (int dimension = 0; dimension < 3; ++dimension) {
    moDipoleMatrix[dimension] = MOMatrix.transpose() * dipoleMatrix_[dimension] * MOMatrix;
  }
  return moDipoleMatrix;
}

template<class NDDOMethod>
Utils::DipoleMatrix NDDODipoleMatrixCalculator<NDDOMethod>::calculateMODipoleMatrixUnrestricted() const {
  Utils::DipoleMatrix moDipoleMatrix;
  auto& alphaMOMatrix = molecularOrbitals_.alphaMatrix();
  auto& betaMOMatrix = molecularOrbitals_.betaMatrix();
  moDipoleMatrix.reset(alphaMOMatrix.cols());
  // Get the MO matrix as D_{MO} = C_\alpha^T * D_{AO} * C_\alpha
  //                             + C_\beta^T  * D_{AO} * C_\beta
  for (int dimension = 0; dimension < 3; ++dimension) {
    auto alphaContribution = alphaMOMatrix.transpose() * dipoleMatrix_[dimension] * alphaMOMatrix;
    auto betaContribution = betaMOMatrix.transpose() * dipoleMatrix_[dimension] * betaMOMatrix;
    moDipoleMatrix[dimension] = alphaContribution + betaContribution;
  }
  return moDipoleMatrix;
}

template class NDDODipoleMatrixCalculator<nddo::PM6Method>;
template class NDDODipoleMatrixCalculator<nddo::AM1Method>;
template class NDDODipoleMatrixCalculator<nddo::MNDOMethod>;

} // namespace Sparrow
} // namespace Scine
