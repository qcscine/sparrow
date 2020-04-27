/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "NDDODipoleMomentCalculator.h"
#include "NDDODipoleMatrixCalculator.h"
#include <Sparrow/Implementations/Nddo/Am1/AM1Method.h>
#include <Sparrow/Implementations/Nddo/Mndo/MNDOMethod.h>
#include <Sparrow/Implementations/Nddo/Pm6/PM6Method.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/multipoleTypes.h>
#include <Sparrow/Implementations/Nddo/Utils/NDDOInitializer.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PrincipalQuantumNumbers.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Geometry.h>
#include <Eigen/Eigenvalues>

namespace Scine {
namespace Sparrow {

template<class NDDOMethod>
NDDODipoleMomentCalculator<NDDOMethod>::~NDDODipoleMomentCalculator() = default;

template<class NDDOMethod>
NDDODipoleMomentCalculator<NDDOMethod>::NDDODipoleMomentCalculator(NDDOMethod& method, DipoleMatrixCalculator& dipoleMatrixCalculator)
  : method_(method), dipoleMatrixCalculator_(dipoleMatrixCalculator) {
}

template<class NDDOMethod>
std::unique_ptr<NDDODipoleMomentCalculator<NDDOMethod>>
NDDODipoleMomentCalculator<NDDOMethod>::create(NDDOMethod& method, DipoleMatrixCalculator& dipoleMatrixCalculator) {
  NDDODipoleMomentCalculator<NDDOMethod> instance(method, dipoleMatrixCalculator);
  return std::make_unique<NDDODipoleMomentCalculator<NDDOMethod>>(std::move(instance));
}

template<class NDDOMethod>
Eigen::RowVector3d NDDODipoleMomentCalculator<NDDOMethod>::calculate() const {
  auto atomicCharges = method_.getAtomicCharges();
  auto coreCharges = method_.getInitializer().getCoreCharges();
  auto densityMatrix = method_.getDensityMatrix().restrictedMatrix();
  auto positions = method_.getPositions();

  if (useNDDOApproximation_) {
    auto elements = method_.getElementTypes();
    std::vector<int> nrAOs(positions.rows());
    std::vector<double> chargeSeparationSP(nrAOs.size());
    std::vector<double> chargeSeparationPD(nrAOs.size());

    for (int el = 0; el < nrAOs.size(); ++el) {
      nrAOs[el] = method_.getAtomsOrbitalsIndexesHolder().getNOrbitals(el);
      chargeSeparationSP[el] =
          method_.getInitializer().getElementParameters().get(elements[el]).chargeSeparations().get(nddo::multipole::sp1);
      chargeSeparationPD[el] =
          method_.getInitializer().getElementParameters().get(elements[el]).chargeSeparations().get(nddo::multipole::pd1);
    }
    return calculateWithNDDOApproximation(std::move(atomicCharges), std::move(positions), std::move(densityMatrix),
                                          std::move(elements), std::move(nrAOs), std::move(chargeSeparationSP),
                                          std::move(chargeSeparationPD));
  }
  else {
    auto overlapMatrix = method_.getOverlapMatrix();
    Eigen::RowVector3d evaluationCoordinate = Utils::Position::Zero();
    if (!dipoleMatrixCalculator_.isValid()) {
      dipoleMatrixCalculator_.fillDipoleMatrix(evaluationCoordinate);
    }
    return calculateWithDipoleMatrix(std::move(coreCharges), std::move(positions), std::move(densityMatrix),
                                     dipoleMatrixCalculator_.getAODipoleMatrix(), std::move(overlapMatrix),
                                     std::move(evaluationCoordinate));
  }
}

template<class NDDOMethod>
Eigen::RowVector3d NDDODipoleMomentCalculator<NDDOMethod>::calculateWithNDDOApproximation(
    std::vector<double> atomicCharges, Utils::PositionCollection positions, Eigen::MatrixXd densityMatrix,
    Utils::ElementTypeCollection elements, std::vector<int> nrAOs, std::vector<double> chargeSeparationSP,
    std::vector<double> chargeSeparationPD) const {
  Eigen::MatrixX3d atomicDipoles(positions.rows(), 3);

  Eigen::RowVector3d totalDipole = Utils::Dipole::Zero(3);

  int currentdensityMatrixOrbital = 0;
  // Put center of mass at the origin
  auto centerOfMass = Utils::Geometry::getCenterOfMass(positions, Utils::Geometry::getMasses(elements));
  positions.rowwise() -= centerOfMass;
  for (unsigned atom = 0; atom < atomicDipoles.rows(); ++atom) {
    // Mullikan charges dipole contribution
    atomicDipoles.row(atom) = atomicCharges[atom] * positions.row(atom);

    // SP hybridization
    if (chargeSeparationSP[atom] != 0) {
      const double HybrSP = 2.0 * chargeSeparationSP[atom];
      for (int pOrbComponent = 0; pOrbComponent < 3; ++pOrbComponent) {
        atomicDipoles.row(atom)(pOrbComponent) -=
            HybrSP * densityMatrix(currentdensityMatrixOrbital, currentdensityMatrixOrbital + 1 + pOrbComponent);
      }
    }
    // PD hybridization
    //                                  0  1   2   3   4         5      6      7      8
    // d orbitals are ordered this way: s, px, py, pz, d(x2-y2), d(xz), d(z2), d(yz), d(xy)
    // where 0 corresponds to currentdensityMatrixOrbital
    // See generalTypes.h for details.
    if (chargeSeparationPD[atom] != 0) {
      const double HybrPD = 2.0 * chargeSeparationPD[atom];
      const double OneOverSqrt3 = 1.0 / std::sqrt(3);

      const double DensMatComponentsX =
          densityMatrix(currentdensityMatrixOrbital + 1, currentdensityMatrixOrbital + 4) -
          densityMatrix(currentdensityMatrixOrbital + 1, currentdensityMatrixOrbital + 6) * OneOverSqrt3 +
          densityMatrix(currentdensityMatrixOrbital + 2, currentdensityMatrixOrbital + 8) +
          densityMatrix(currentdensityMatrixOrbital + 3, currentdensityMatrixOrbital + 5);

      atomicDipoles.row(atom)(0) -= DensMatComponentsX * HybrPD;

      const double DensMatComponentsY =
          densityMatrix(currentdensityMatrixOrbital + 1, currentdensityMatrixOrbital + 8) -
          densityMatrix(currentdensityMatrixOrbital + 2, currentdensityMatrixOrbital + 4) -
          densityMatrix(currentdensityMatrixOrbital + 2, currentdensityMatrixOrbital + 6) * OneOverSqrt3 +
          densityMatrix(currentdensityMatrixOrbital + 3, currentdensityMatrixOrbital + 7);

      atomicDipoles.row(atom)(1) -= DensMatComponentsY * HybrPD;

      const double DensMatComponentsZ =
          densityMatrix(currentdensityMatrixOrbital + 1, currentdensityMatrixOrbital + 5) +
          densityMatrix(currentdensityMatrixOrbital + 2, currentdensityMatrixOrbital + 7) +
          densityMatrix(currentdensityMatrixOrbital + 3, currentdensityMatrixOrbital + 6) * 2.0 * OneOverSqrt3;

      atomicDipoles.row(atom)(2) -= DensMatComponentsZ * HybrPD;
    }

    currentdensityMatrixOrbital += nrAOs[atom];
  }

  totalDipole += atomicDipoles.colwise().sum();
  return totalDipole;
}

template<class NDDOMethod>
void NDDODipoleMomentCalculator<NDDOMethod>::useNDDOApproximation(bool useNDDOApprox) {
  useNDDOApproximation_ = useNDDOApprox;
}

template<class NDDOMethod>
Eigen::RowVector3d NDDODipoleMomentCalculator<NDDOMethod>::calculateWithDipoleMatrix(
    std::vector<double> coreCharges, Utils::PositionCollection positions, Eigen::MatrixXd densityMatrix,
    Utils::DipoleMatrix dipoleMatrix, Eigen::MatrixXd overlapMatrix, Eigen::RowVector3d dipoleEvaluationCoordinate) const {
  auto const nAtoms = coreCharges.size();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(overlapMatrix);
  Eigen::MatrixXd sqrtS = es.operatorInverseSqrt();
  Eigen::MatrixXd nonOrthogonalDensityMatrix =
      sqrtS.selfadjointView<Eigen::Upper>() * densityMatrix * sqrtS.selfadjointView<Eigen::Upper>();
  Eigen::RowVector3d dipole;
  dipole.setZero();

  positions.rowwise() -= dipoleEvaluationCoordinate;

  const auto atomicOrbitalsNumber = dipoleMatrix.x().get<Utils::derivOrder::zero>().cols();
  // Classical Core charge component
  for (int atom = 0; atom < nAtoms; ++atom) {
    dipole += coreCharges[atom] * positions.row(atom);
  }

  // electronic component (diagonal)
  for (int AO = 0; AO < atomicOrbitalsNumber; ++AO) {
    dipole.x() -= nonOrthogonalDensityMatrix(AO, AO) * dipoleMatrix.x().get<Utils::derivOrder::zero>()(AO, AO);
    dipole.y() -= nonOrthogonalDensityMatrix(AO, AO) * dipoleMatrix.y().get<Utils::derivOrder::zero>()(AO, AO);
    dipole.z() -= nonOrthogonalDensityMatrix(AO, AO) * dipoleMatrix.z().get<Utils::derivOrder::zero>()(AO, AO);
  }
  // electronic component(off-diagonal)
  for (int nu = 0; nu < atomicOrbitalsNumber; ++nu) {
    for (int mu = nu + 1; mu < atomicOrbitalsNumber; ++mu) {
      dipole.x() -= 2 * nonOrthogonalDensityMatrix(nu, mu) * dipoleMatrix.x().get<Utils::derivOrder::zero>()(nu, mu);
      dipole.y() -= 2 * nonOrthogonalDensityMatrix(nu, mu) * dipoleMatrix.y().get<Utils::derivOrder::zero>()(nu, mu);
      dipole.z() -= 2 * nonOrthogonalDensityMatrix(nu, mu) * dipoleMatrix.z().get<Utils::derivOrder::zero>()(nu, mu);
    }
  }
  return dipole;
}

template class NDDODipoleMomentCalculator<nddo::PM6Method>;
template class NDDODipoleMomentCalculator<nddo::AM1Method>;
template class NDDODipoleMomentCalculator<nddo::MNDOMethod>;
} // namespace Sparrow
} // namespace Scine
