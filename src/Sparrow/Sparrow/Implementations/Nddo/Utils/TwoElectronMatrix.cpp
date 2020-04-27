/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TwoElectronMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/Global2c2eMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/TwoCenterIntegralContainer.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterIntegralContainer.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <omp.h>

namespace Scine {
namespace Sparrow {
namespace nddo {

#pragma omp declare reduction (+: Eigen::MatrixXd: omp_out=omp_out+omp_in)\
initializer(omp_priv=Eigen::MatrixXd::Zero(omp_orig.rows(), omp_orig.cols()))

using namespace Utils::AutomaticDifferentiation;

TwoElectronMatrix::TwoElectronMatrix(const Utils::ElementTypeCollection& elements, const Utils::DensityMatrix& densityMatrix,
                                     const OneCenterIntegralContainer& oneCIntegrals,
                                     const TwoCenterIntegralContainer& twoCIntegrals,
                                     const ElementParameters& elementPar, const Utils::AtomsOrbitalsIndexes& aoIndexes)
  : P(densityMatrix.restrictedMatrix()),
    PAlpha_(densityMatrix.alphaMatrix()),
    PBeta_(densityMatrix.betaMatrix()),
    oneCenterIntegrals(oneCIntegrals),
    twoCenterIntegrals(twoCIntegrals),
    elementParameters(elementPar),
    aoIndexes_(aoIndexes),
    elementTypes_(elements) {
}

void TwoElectronMatrix::initialize() {
  nAOs_ = 0;
  nAtoms_ = static_cast<int>(elementTypes_.size());
  for (auto e : elementTypes_)
    nAOs_ += elementParameters.get(e).nAOs();
  G_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
}

void TwoElectronMatrix::calculate(bool spinPolarized) {
  spinPolarized_ = spinPolarized;
  if (!spinPolarized_) {
    G_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
  }
  else {
    GAlpha_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
    GBeta_ = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
  }
  calculateBlocks();
}

void TwoElectronMatrix::calculateBlocks() {
  for (int i = 0; i < nAtoms_; ++i) {
    auto index = aoIndexes_.getFirstOrbitalIndex(i);
    auto nAOs = aoIndexes_.getNOrbitals(i);
    calculateSameAtomBlock(index, nAOs, elementTypes_[i], G_, GAlpha_, GBeta_);
  }

  for (int i = 0; i < nAtoms_; ++i) {
    auto indexA = aoIndexes_.getFirstOrbitalIndex(i);
    auto nAOsA = aoIndexes_.getNOrbitals(i);

    for (int j = i + 1; j < nAtoms_; j++) {
      auto indexB = aoIndexes_.getFirstOrbitalIndex(j);
      auto nAOsB = aoIndexes_.getNOrbitals(j);
      calculateDifferentAtomsBlock(indexA, indexB, nAOsA, nAOsB, *twoCenterIntegrals.get(i, j), G_, GAlpha_, GBeta_);
    }
  }
}

void TwoElectronMatrix::calculateSameAtomBlock(int startIndex, int nAOs, Utils::ElementType el, Eigen::MatrixXd& G,
                                               Eigen::MatrixXd& GAlpha, Eigen::MatrixXd& GBeta) {
  const auto& m = oneCenterIntegrals.get(el);
  for (int i = 0; i < nAOs; i++) {
    for (int j = 0; j <= i; j++) {
      // The following is only valid when only s and p orbitals are involved:
      // G_(startIndex+i, startIndex+i) += densityMatrix_(startIndex+j,startIndex+j) *
      // (m->get(i,i,j,j)-0.5*m->get(i,j,i,j));  if(i!=j) G_(startIndex+i, startIndex+j) +=
      // 0.5*densityMatrix_(startIndex+i, startIndex+j)*(3*m->get(i,j,i,j)-m->get(i,i,j,j));
      for (int k = 0; k < nAOs; k++) {
        for (int l = 0; l < nAOs; l++) {
          // From Thiel, Perspectives on Semiempirical Molecular Orbital Theory
          if (!spinPolarized_) {
            G(startIndex + i, startIndex + j) +=
                P(startIndex + k, startIndex + l) * (m.get(i, j, k, l) - 0.5 * m.get(i, k, j, l));
          }
          else {
            GAlpha(startIndex + i, startIndex + j) += P(startIndex + k, startIndex + l) * m.get(i, j, k, l) -
                                                      PAlpha_(startIndex + k, startIndex + l) * m.get(i, k, j, l);
            GBeta(startIndex + i, startIndex + j) += P(startIndex + k, startIndex + l) * m.get(i, j, k, l) -
                                                     PBeta_(startIndex + k, startIndex + l) * m.get(i, k, j, l);
          }
        }
      }
    }
  }
}
void TwoElectronMatrix::calculateDifferentAtomsBlock(int startA, int startB, int nAOsA, int nAOsB,
                                                     const multipole::Global2c2eMatrix& m, Eigen::MatrixXd& G,
                                                     Eigen::MatrixXd& GAlpha, Eigen::MatrixXd& GBeta) {
  int mu, nu, lambda, sigma;
  for (int i = 0; i < nAOsA; i++) {
    mu = startA + i;
    for (int j = 0; j <= i; j++) {
      nu = startA + j;

      for (int k = 0; k < nAOsB; k++) {
        lambda = startB + k;

        for (int l = 0; l <= k; l++) {
          sigma = startB + l;

          double integral = m.get(i, j, k, l);
          int multiplicityA = (mu == nu) ? 1 : 2;
          int multiplicityB = (lambda == sigma) ? 1 : 2;

          if (!spinPolarized_) {
            G(mu, nu) += P(lambda, sigma) * integral * multiplicityB;
            G(lambda, sigma) += P(mu, nu) * integral * multiplicityA;
            G(lambda, mu) += -0.5 * P(sigma, nu) * integral;
            if (mu > nu) {
              G(lambda, nu) += -0.5 * P(sigma, mu) * integral;
              if (lambda > sigma) {
                G(sigma, nu) += -0.5 * P(lambda, mu) * integral;
              }
            }
            if (lambda > sigma) {
              G(sigma, mu) += -0.5 * P(lambda, nu) * integral;
            }
          }
          else {
            GAlpha(mu, nu) += P(lambda, sigma) * integral * multiplicityB;
            GBeta(mu, nu) += P(lambda, sigma) * integral * multiplicityB;
            GAlpha(lambda, sigma) += P(mu, nu) * integral * multiplicityA;
            GBeta(lambda, sigma) += P(mu, nu) * integral * multiplicityA;
            GAlpha(lambda, mu) -= PAlpha_(sigma, nu) * integral;
            GBeta(lambda, mu) -= PBeta_(sigma, nu) * integral;
            if (mu > nu) {
              GAlpha(lambda, nu) -= PAlpha_(sigma, mu) * integral;
              GBeta(lambda, nu) -= PBeta_(sigma, mu) * integral;
              if (lambda > sigma) {
                GAlpha(sigma, nu) -= PAlpha_(lambda, mu) * integral;
                GBeta(sigma, nu) -= PBeta_(lambda, mu) * integral;
              }
            }
            if (lambda > sigma) {
              GAlpha(sigma, mu) -= PAlpha_(lambda, nu) * integral;
              GBeta(sigma, mu) -= PBeta_(lambda, nu) * integral;
            }
          }
        }
      }
    }
  }
}

template<Utils::derivativeType O>
void TwoElectronMatrix::addDerivatives(DerivativeContainerType<O>& derivativeContainer) const {
  for (int i = 0; i < nAtoms_; ++i) {
    auto indexA = aoIndexes_.getFirstOrbitalIndex(i);
    auto nAOsA = aoIndexes_.getNOrbitals(i);
    for (int j = i + 1; j < nAtoms_; j++) {
      auto indexB = aoIndexes_.getFirstOrbitalIndex(j);
      auto nAOsB = aoIndexes_.getNOrbitals(j);

      addDerivativesForBlock<O>(derivativeContainer, i, j, indexA, indexB, nAOsA, nAOsB, *twoCenterIntegrals.get(i, j));
    }
  }
}

template<Utils::derivativeType O>
void TwoElectronMatrix::addDerivativesForBlock(DerivativeContainerType<O>& derivativeContainer, int a, int b, int startA,
                                               int startB, int nAOsA, int nAOsB, const multipole::Global2c2eMatrix& m) const {
  int mu, nu, lambda, sigma;
  Utils::AutomaticDifferentiation::DerivativeType<O> integralD;
  Utils::AutomaticDifferentiation::DerivativeType<O> derivativeContribution;
  integralD.setZero();
  derivativeContribution.setZero();
  for (int i = 0; i < nAOsA; i++) {
    mu = startA + i;

    for (int j = 0; j <= i; j++) {
      nu = startA + j;

      for (int k = 0; k < nAOsB; k++) {
        lambda = startB + k;

        for (int l = 0; l <= k; l++) {
          sigma = startB + l;

          integralD = m.getDerivative<O>(i, j, k, l);

          // this multiple contains, in addition to what was above, a 1/2 factor for diagonal elements (NB: but not
          // visible because of two equal terms)
          double multiple = ((mu == nu) ? 1 : 2) * ((lambda == sigma) ? 1 : 2);
          double factor;

          if (!spinPolarized_) {
            factor = P(lambda, sigma) * P(mu, nu) * multiple - 0.5 * P(sigma, nu) * P(lambda, mu);
            if (mu > nu) {
              factor += (-0.5 * P(sigma, mu) * P(lambda, nu));
              if (lambda > sigma) {
                factor += (-0.5 * P(lambda, mu) * P(sigma, nu));
              }
            }
            if (lambda > sigma) {
              factor += (-0.5 * P(lambda, nu) * P(sigma, mu));
            }
          }
          else {
            factor = P(lambda, sigma) * P(mu, nu) * multiple - PAlpha_(sigma, nu) * PAlpha_(lambda, mu) -
                     PBeta_(sigma, nu) * PBeta_(lambda, mu);
            if (mu > nu) {
              factor -= PAlpha_(sigma, mu) * PAlpha_(lambda, nu) + PBeta_(sigma, mu) * PBeta_(lambda, nu);
              if (lambda > sigma) {
                factor -= PAlpha_(lambda, mu) * PAlpha_(sigma, nu) + PBeta_(lambda, mu) * PBeta_(sigma, nu);
              }
            }
            if (lambda > sigma) {
              factor -= PAlpha_(lambda, nu) * PAlpha_(sigma, mu) + PBeta_(lambda, nu) * PBeta_(sigma, mu);
            }
          }

          derivativeContribution += integralD * factor;
        }
      }
    }
  }
  Utils::AutomaticDifferentiation::addDerivativeToContainer<O>(derivativeContainer, a, b, derivativeContribution);
}

template void
TwoElectronMatrix::addDerivatives<Utils::derivativeType::first>(DerivativeContainerType<Utils::derivativeType::first>&) const;
template void TwoElectronMatrix::addDerivatives<Utils::derivativeType::second_atomic>(
    DerivativeContainerType<Utils::derivativeType::second_atomic>&) const;
template void TwoElectronMatrix::addDerivatives<Utils::derivativeType::second_full>(
    DerivativeContainerType<Utils::derivativeType::second_full>&) const;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
