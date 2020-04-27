/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OneElectronMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/TwoCenterIntegralContainer.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/VuvB.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <algorithm>

namespace Scine {
namespace Sparrow {
using namespace Utils::AutomaticDifferentiation;

namespace nddo {

OneElectronMatrix::OneElectronMatrix(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                                     const Eigen::MatrixXd& densityMatrix, const TwoCenterIntegralContainer& twoCIntegrals,
                                     const ElementParameters& elementPar, const Utils::AtomsOrbitalsIndexes& aoIndexes)
  : P(densityMatrix),
    twoCenterIntegrals(twoCIntegrals),
    elementParameters(elementPar),
    aoIndexes_(aoIndexes),
    elementTypes_(elements),
    positions_(positions) {
}

void OneElectronMatrix::initialize() {
  nAOs_ = 0;
  nAtoms_ = elementTypes_.size();
  for (auto e : elementTypes_)
    nAOs_ += elementParameters.get(e).nAOs();

  H_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
}

void OneElectronMatrix::calculate(const Utils::MatrixWithDerivatives& S) {
  H_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
#pragma omp parallel
  {
    calculateSameAtomBlocks();
    calculateDifferentAtomsBlocks(S);
  }
}

void OneElectronMatrix::calculateSameAtomBlocks() {
#pragma omp for schedule(static) nowait
  for (int i = 0; i < nAtoms_; ++i) {
    auto index = aoIndexes_.getFirstOrbitalIndex(i);
    auto nAOs = aoIndexes_.getNOrbitals(i);
    calculateSameAtomBlock(i, index, nAOs);
  }
}

void OneElectronMatrix::calculateSameAtomBlock(int a, int startIndex, int nAOs) {
  const auto& ap = elementParameters.get(elementTypes_[a]);

  for (int i = 0; i < nAOs; i++) {
#pragma omp atomic write
    H_(startIndex + i, startIndex + i) = (i < 1) ? ap.Uss() : (i < 4) ? ap.Upp() : ap.Udd();
  }

  multipole::VuvB V_((nAOs == 1) ? 0 : (nAOs == 4) ? 1 : 2);
  for (int b = 0; b < nAtoms_; b++) {
    if (b != a) {
      const auto& pB = elementParameters.get(elementTypes_[b]);
      if (pB.pCoreSpecified()) {
        auto posA = positions_.row(a);
        auto posB = positions_.row(b);
        Eigen::RowVector3d Rab = posB - posA;

        V_.calculate<Utils::derivOrder::zero>(Rab, ap.chargeSeparations(), ap.klopmanParameters(), pB.pCore(),
                                              pB.coreCharge());
        for (int i = 0; i < nAOs; i++) {
          for (int j = 0; j <= i; j++) {
#pragma omp atomic
            H_(startIndex + i, startIndex + j) += V_.get(i, j);
          }
        }
      }
      else {
        if (a < b) {
          const auto& m = *twoCenterIntegrals.get(a, b);
          for (int i = 0; i < nAOs; i++) {
            for (int j = 0; j <= i; j++) {
#pragma omp atomic
              H_(startIndex + i, startIndex + j) += -pB.coreCharge() * m.get(i, j, 0, 0);
            }
          }
        }
        else {
          const auto& m = *twoCenterIntegrals.get(b, a);
          for (int i = 0; i < nAOs; i++) {
            for (int j = 0; j <= i; j++) {
#pragma omp atomic
              H_(startIndex + i, startIndex + j) += -pB.coreCharge() * m.get(0, 0, i, j);
            }
          }
        }
      }
    }
  }
}

void OneElectronMatrix::calculateDifferentAtomsBlocks(const Utils::MatrixWithDerivatives& S) {
#pragma omp for schedule(static)
  for (int i = 1; i < nAtoms_; ++i) {
    auto rowIndex = aoIndexes_.getFirstOrbitalIndex(i);
    const auto& pA = elementParameters.get(elementTypes_[i]);
    for (int j = 0; j < i; j++) {
      const auto& pB = elementParameters.get(elementTypes_[j]);

      auto colIndex = aoIndexes_.getFirstOrbitalIndex(j);
      calculateDifferentAtomsBlock(rowIndex, colIndex, pA, pB, S);
    }
  }
}

void OneElectronMatrix::calculateDifferentAtomsBlock(int startRow, int startCol, const AtomicParameters& pA,
                                                     const AtomicParameters& pB, const Utils::MatrixWithDerivatives& S) {
  for (int j = 0; j < pB.nAOs(); j++) {
    for (int i = 0; i < pA.nAOs(); i++) {
      double betaA = (i < 1) ? pA.betaS() : (i < 4) ? pA.betaP() : pA.betaD();
      double betaB = (j < 1) ? pB.betaS() : (j < 4) ? pB.betaP() : pB.betaD();
      H_(startRow + i, startCol + j) = 0.5 * (betaA + betaB) * S.getValue(startRow + i, startCol + j);
    }
  }
}

template<Utils::derivativeType O>
void OneElectronMatrix::addDerivatives(DerivativeContainerType<O>& derivativeContainer,
                                       const Utils::MatrixWithDerivatives& S) const {
  for (int i = 0; i < nAtoms_; ++i) {
    auto index = aoIndexes_.getFirstOrbitalIndex(i);
    auto nAOs = aoIndexes_.getNOrbitals(i);
    addDerivativesContribution1<O>(derivativeContainer, i, index, nAOs);
  }
  for (int a = 1; a < nAtoms_; ++a) {
    for (int b = 0; b < a; b++) {
      auto indexA = aoIndexes_.getFirstOrbitalIndex(a);
      auto nAOsA = aoIndexes_.getNOrbitals(a);
      auto indexB = aoIndexes_.getFirstOrbitalIndex(b);
      auto nAOsB = aoIndexes_.getNOrbitals(b);

      addDerivativesContribution2<O>(derivativeContainer, a, b, indexA, indexB, nAOsA, nAOsB, S);
    }
  }
}

template<Utils::derivativeType O>
void OneElectronMatrix::addDerivativesContribution1(DerivativeContainerType<O>& derivativeContainer, int a,
                                                    int startIndex, int nAOs) const {
  const auto& ap = elementParameters.get(elementTypes_[a]);

  multipole::VuvB V_((nAOs == 1) ? 0 : (nAOs == 4) ? 1 : 2);
  for (int b = 0; b < nAtoms_; b++) {
    double Pij = 0.0;
    DerivativeType<O> contrib;
    contrib.setZero();
    if (b != a) {
      const auto& pB = elementParameters.get(elementTypes_[b]);
      if (pB.pCoreSpecified()) {
        auto posA = positions_.row(a);
        auto posB = positions_.row(b);
        auto Rab = posB - posA;

        V_.calculate<UnderlyingOrder<O>>(Rab, ap.chargeSeparations(), ap.klopmanParameters(), pB.pCore(), pB.coreCharge());
        for (int i = 0; i < nAOs; i++) {
          for (int j = 0; j <= i; j++) {
            Pij = P(startIndex + i, startIndex + j);
            auto vd = V_.getDerivative<O>(i, j);
            contrib += vd * Pij * (i == j ? 1 : 2);
          }
        }
      }
      else {
        if (a < b) {
          const auto& m = *twoCenterIntegrals.get(a, b);
          for (int i = 0; i < nAOs; i++) {
            for (int j = 0; j <= i; j++) {
              Pij = P(startIndex + i, startIndex + j);
              auto vd = -pB.coreCharge() * m.getDerivative<O>(i, j, 0, 0);
              contrib += vd * Pij * (i == j ? 1 : 2);
            }
          }
        }
        else {
          const auto& m = *twoCenterIntegrals.get(b, a);
          for (int i = 0; i < nAOs; i++) {
            for (int j = 0; j <= i; j++) {
              Pij = P(startIndex + i, startIndex + j);
              auto vd = -pB.coreCharge() * m.getDerivative<O>(0, 0, i, j);
              // other way around, since swapped a-b
              contrib += getOppositeDerivative<O>(vd * Pij * (i == j ? 1 : 2));
            }
          }
        }
      }
#pragma omp critical(addDerivativesOneElectronMatrix)
      { addDerivativeToContainer<O>(derivativeContainer, a, b, contrib); }
    }
  }
}

template<Utils::derivativeType O>
void OneElectronMatrix::addDerivativesContribution2(DerivativeContainerType<O>& derivativeContainer, int a, int b,
                                                    int indexA, int indexB, int nAOsA, int nAOsB,
                                                    const Utils::MatrixWithDerivatives& S) const {
  const auto& pA = elementParameters.get(elementTypes_[a]);
  const auto& pB = elementParameters.get(elementTypes_[b]);
  DerivativeType<O> derivativeContribution;
  derivativeContribution.setZero();
  for (int i = 0; i < nAOsA; i++) {
    double betaA = (i < 1) ? pA.betaS() : (i < 4) ? pA.betaP() : pA.betaD();
    for (int j = 0; j < nAOsB; j++) {
      double betaB = (j < 1) ? pB.betaS() : (j < 4) ? pB.betaP() : pB.betaD();
      derivativeContribution +=
          getDerivativeFromValueWithDerivatives<O>(S.get<UnderlyingOrder<O>>()(indexA + i, indexB + j)) *
          ((betaA + betaB) * P(indexA + i, indexB + j));
    }
  }
  addDerivativeToContainer<O>(derivativeContainer, a, b, derivativeContribution);
}

template void
OneElectronMatrix::addDerivatives<Utils::derivativeType::first>(DerivativeContainerType<Utils::derivativeType::first>&,
                                                                const Utils::MatrixWithDerivatives&) const;
template void OneElectronMatrix::addDerivatives<Utils::derivativeType::second_atomic>(
    DerivativeContainerType<Utils::derivativeType::second_atomic>&, const Utils::MatrixWithDerivatives&) const;
template void OneElectronMatrix::addDerivatives<Utils::derivativeType::second_full>(
    DerivativeContainerType<Utils::derivativeType::second_full>&, const Utils::MatrixWithDerivatives&) const;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
