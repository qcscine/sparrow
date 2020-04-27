/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SDFTB.h"
#include "DFTBCommon.h"
#include "SKAtom.h"
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/DataStructures/SpinAdaptedMatrix.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Math/AutomaticDifferentiation/TypeDefinitions.h>
#include <utility>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace dftb {

SDFTB::SDFTB(const Utils::ElementTypeCollection& elements, const std::vector<std::unique_ptr<SKAtom>>& atomicParameters)
  : atomParameters(atomicParameters), elementTypes_(elements) {
}

void SDFTB::initialize(int nAtoms, int nAOs, Utils::AtomsOrbitalsIndexes indexes) {
  nAtoms_ = nAtoms;
  nAOs_ = nAOs;
  aoIndexes_ = std::move(indexes);

  pup = std::vector<double>(3 * nAtoms);
  pdn = std::vector<double>(3 * nAtoms);
  pdif = std::vector<double>(3 * nAtoms);
  spinContribution_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
}

double SDFTB::spinEnergyContribution() const {
  double spinEnergy = 0.0;

  for (int a = 0; a < nAtoms_; a++) {
    Utils::ElementType el = elementTypes_[a];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        spinEnergy += pdif[3 * a + i] * pdif[3 * a + j] * atomParameters[Utils::ElementInfo::Z(el)]->getSpinConstant(i, j);
  }
  spinEnergy *= 0.5;
  return spinEnergy;
}

template<Utils::derivativeType O>
void SDFTB::addDerivatives(DerivativeContainerType<O>& derivativesContainer, const Utils::MatrixWithDerivatives& overlapDeriv,
                           const Eigen::MatrixXd& pUp, const Eigen::MatrixXd& pDn) const {
#pragma omp parallel for
  for (int a = 1; a < nAtoms_; ++a) {
    int nAOsA = aoIndexes_.getNOrbitals(a);
    int indexA = aoIndexes_.getFirstOrbitalIndex(a);

    for (int b = 0; b < a; ++b) {
      int nAOsB = aoIndexes_.getNOrbitals(b);
      int indexB = aoIndexes_.getFirstOrbitalIndex(b);

      DerivativeType<O> contribution;
      contribution.setZero();
      for (int i = indexA; i < indexA + nAOsA; ++i) {
        for (int j = indexB; j < indexB + nAOsB; ++j) {
          contribution += 2 * (pUp(i, j) - pDn(i, j)) * spinContribution_(i, j) *
                          getDerivativeFromValueWithDerivatives<O>(overlapDeriv.get<UnderlyingOrder<O>>()(i, j));
        }
      }
#pragma omp critical
      { addDerivativeToContainer<O>(derivativesContainer, a, b, contribution); }
    }
  }
}

void SDFTB::spinPopulationAnalysis(const Eigen::MatrixXd& densityMatrixUp, const Eigen::MatrixXd& densityMatrixDn,
                                   const Eigen::MatrixXd& overlapMatrix) {
  spinPopulationAnalysis(densityMatrixUp, overlapMatrix, pup);
  spinPopulationAnalysis(densityMatrixDn, overlapMatrix, pdn);

  // get the spin population as the difference between up- and down- populations
  for (int i = 0; i < 3 * nAtoms_; i++)
    pdif[i] = pup[i] - pdn[i];
}

// Given a density matrix and an overlap matrix for spin-up or spin-down problem,
// performs the calculation of the spin population of the corresponding angular momentum state
void SDFTB::spinPopulationAnalysis(const Eigen::MatrixXd& densityMatrix, const Eigen::MatrixXd& overlapMatrix,
                                   std::vector<double>& population) {
  // Calculation of the Mulliken charges as described in elstner1998
  Eigen::MatrixXd D = densityMatrix.cwiseProduct(overlapMatrix); // population matrix

  // Empty vector:
  for (double& i : population)
    i = 0.0;

  int l;
  for (int a = 0; a < nAtoms_; a++) {
    int nAOsA = aoIndexes_.getNOrbitals(a);
    int index = aoIndexes_.getFirstOrbitalIndex(a);
    for (int mu = 0; mu < nAOsA; mu++) {
      if (mu == 0) // s orbital
        l = 0;
      else if (mu <= 4) // p orbital
        l = 1;
      else // d orbital
        l = 2;
      for (int nu = 0; nu < nAOs_; nu++)
        population[3 * a + l] += D(index + mu, nu);
    }
  }
}

// Add the diagonal contribution to the up/down Hamiltonians
// according to frauenheim2000
void SDFTB::calculateSpinContribution() {
  spinContribution_.setZero();
#pragma omp parallel
  {
    // Update diagonal elements
#pragma omp for nowait
    for (int a = 0; a < nAtoms_; a++) {
      Utils::ElementType el = elementTypes_[a];
      SKAtom* par = atomParameters[Utils::ElementInfo::Z(el)].get();
      int nAOsA = aoIndexes_.getNOrbitals(a);
      int indexA = aoIndexes_.getFirstOrbitalIndex(a);

      // ss
      double spinTerm = 0.0;
      for (int l = 0; l < 3; l++) {
        spinTerm += par->getSpinConstant(0, l) * pdif[3 * a + l];
      }
      // S_munu=1
#pragma omp atomic
      spinContribution_(indexA, indexA) += spinTerm;

      if (nAOsA > 1) {
        // pp
        spinTerm = 0.0;
        for (int l = 0; l < 3; l++) {
          spinTerm += par->getSpinConstant(1, l) * pdif[3 * a + l];
        }
        // S_munu=1
        for (int i = 1; i < 4; i++) {
#pragma omp atomic
          spinContribution_(indexA + i, indexA + i) += spinTerm;
        }

        if (nAOsA > 4) {
          // dd
          spinTerm = 0.0;
          for (int l = 0; l < 3; l++) {
            spinTerm += par->getSpinConstant(2, l) * pdif[3 * a + l];
          }
          // S_munu=1
          for (int i = 4; i < 9; i++) {
#pragma omp atomic
            spinContribution_(indexA + i, indexA + i) += spinTerm;
          }
        }
      }
    }

    // Update non-diagonal elements
#pragma omp for
    for (int a = 1; a < nAtoms_; ++a) {
      Utils::ElementType elA = elementTypes_[a];
      SKAtom* parA = atomParameters[Utils::ElementInfo::Z(elA)].get();
      int nAOsA = aoIndexes_.getNOrbitals(a);
      int indexA = aoIndexes_.getFirstOrbitalIndex(a);

      for (int b = 0; b < a; ++b) {
        Utils::ElementType elB = elementTypes_[b];
        SKAtom* parB = atomParameters[Utils::ElementInfo::Z(elB)].get();
        int nAOsB = aoIndexes_.getNOrbitals(b);
        int indexB = aoIndexes_.getFirstOrbitalIndex(b);

        // CONTRIBUTION FROM A
        // ss
        double spinTerm = 0.0;
        for (int l = 0; l < 3; l++) {
          spinTerm += parA->getSpinConstant(0, l) * pdif[3 * a + l];
        }
        for (auto i = 0; i < nAOsB; ++i) {
#pragma omp atomic
          spinContribution_(indexA, indexB + i) += 0.5 * spinTerm;
        }

        if (nAOsA > 1) {
          // pp
          spinTerm = 0.0;
          for (int l = 0; l < 3; l++) {
            spinTerm += parA->getSpinConstant(1, l) * pdif[3 * a + l];
          }
          // S_munu=1
          for (int i = 1; i < 4; i++) {
            for (auto j = 0; j < nAOsB; ++j)
#pragma omp atomic
              spinContribution_(indexA + i, indexB + j) += 0.5 * spinTerm;
          }

          if (nAOsA > 4) {
            // dd
            spinTerm = 0.0;
            for (int l = 0; l < 3; l++) {
              spinTerm += parA->getSpinConstant(2, l) * pdif[3 * a + l];
            }
            // S_munu=1
            for (int i = 4; i < 9; i++) {
              for (auto j = 0; j < nAOsB; ++j)
#pragma omp atomic
                spinContribution_(indexA + i, indexB + j) += 0.5 * spinTerm;
            }
          }
        }

        // CONTRIBUTION FROM B
        // ss
        spinTerm = 0.0;
        for (int l = 0; l < 3; l++) {
          spinTerm += parB->getSpinConstant(0, l) * pdif[3 * b + l];
        }
        for (auto i = 0; i < nAOsA; ++i) {
#pragma omp atomic
          spinContribution_(indexA + i, indexB) += 0.5 * spinTerm;
        }

        if (nAOsB > 1) {
          // pp
          spinTerm = 0.0;
          for (int l = 0; l < 3; l++) {
            spinTerm += parB->getSpinConstant(1, l) * pdif[3 * b + l];
          }
          // S_munu=1
          for (int i = 1; i < 4; i++) {
            for (auto j = 0; j < nAOsA; ++j)
#pragma omp atomic
              spinContribution_(indexA + j, indexB + i) += 0.5 * spinTerm;
          }

          if (nAOsB > 4) {
            // dd
            spinTerm = 0.0;
            for (int l = 0; l < 3; l++) {
              spinTerm += parB->getSpinConstant(2, l) * pdif[3 * b + l];
            }
            // S_munu=1
            for (int i = 4; i < 9; i++) {
              for (auto j = 0; j < nAOsA; ++j)
#pragma omp atomic
                spinContribution_(indexA + j, indexB + i) += 0.5 * spinTerm;
            }
          }
        }
      }
    }
  }
}

void SDFTB::constructSpinHamiltonians(Utils::SpinAdaptedMatrix& H, const Eigen::MatrixXd& overlap) const {
  Eigen::MatrixXd Hup = H.restrictedMatrix();
  Eigen::MatrixXd Hdn = H.restrictedMatrix();

#pragma omp parallel for simd
  for (int i = 0; i < nAOs_; ++i) {
    for (int j = 0; j <= i; ++j) {
      Hup(i, j) += spinContribution_(i, j) * overlap(i, j);
      Hdn(i, j) -= spinContribution_(i, j) * overlap(i, j);

      if (i != j) {
        Hup(j, i) += spinContribution_(i, j) * overlap(i, j);
        Hdn(j, i) -= spinContribution_(i, j) * overlap(i, j);
      }
    }
  }

  H.setAlphaMatrix(std::move(Hup));
  H.setBetaMatrix(std::move(Hdn));
}

template void SDFTB::addDerivatives<Utils::derivativeType::first>(DerivativeContainerType<Utils::derivativeType::first>&,
                                                                  const Utils::MatrixWithDerivatives&,
                                                                  const Eigen::MatrixXd&, const Eigen::MatrixXd&) const;
template void
SDFTB::addDerivatives<Utils::derivativeType::second_atomic>(DerivativeContainerType<Utils::derivativeType::second_atomic>&,
                                                            const Utils::MatrixWithDerivatives&, const Eigen::MatrixXd&,
                                                            const Eigen::MatrixXd&) const;
template void
SDFTB::addDerivatives<Utils::derivativeType::second_full>(DerivativeContainerType<Utils::derivativeType::second_full>&,
                                                          const Utils::MatrixWithDerivatives&, const Eigen::MatrixXd&,
                                                          const Eigen::MatrixXd&) const;

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
