/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CISMatrixAOFockBuilder.h"
#include "Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/Global2c2eMatrix.h"
#include "Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h"
namespace Scine {
namespace Sparrow {

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock>
CISMatrixAOFockBuilder<restrictedness, spinBlock>::CISMatrixAOFockBuilder(CISData cisData,
                                                                          const ExcitedStatesParam& excitedStatesParam)
  : cisData_(std::move(cisData)) {
  c1_ = excitedStatesParam.c1;
  c2_ = excitedStatesParam.c2;
  coulombContainer_ = std::make_shared<std::vector<std::map<int, std::shared_ptr<Eigen::MatrixXd>>>>();
  exchangeContainer_ = std::make_shared<std::vector<std::map<int, std::shared_ptr<Eigen::MatrixXd>>>>();
  initialize();
  calculateMatrices();
}

template<>
void CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Singlet>::initialize() {
  coulombContainer_->resize(nAtoms_);
  exchangeContainer_->resize(nAtoms_);
}
template<>
void CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Triplet>::initialize() {
  exchangeContainer_->resize(nAtoms_);
}

template<>
void CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Singlet>::initialize() {
  coulombContainer_->resize(nAtoms_);
  exchangeContainer_->resize(nAtoms_);
}
template<>
void CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Triplet>::initialize() {
  throw std::runtime_error(" CISMatrixAOFockBuilder: Invalid combination: Calculation for an unrestricted reference "
                           "with a triplet spin-transition not possible!.");
}

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock>
void CISMatrixAOFockBuilder<restrictedness, spinBlock>::calculateMatrices() {
  for (int atomI = 0; atomI < nAtoms_; ++atomI) {
    calculate(atomI);
    for (int atomJ = atomI + 1; atomJ < nAtoms_; ++atomJ) {
      calculate(atomI, atomJ);
    }
  }
}

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock>
inline Eigen::MatrixXd CISMatrixAOFockBuilder<restrictedness, spinBlock>::getCoulombIntegrals(int atomI, int atomJ) {
  int nAOsI = cisData_.AOInfo.getNOrbitals(atomI);
  int nAOsJ = cisData_.AOInfo.getNOrbitals(atomJ);
  Eigen::MatrixXd integralMatrix = Eigen::MatrixXd::Zero(nAOsI * nAOsI, nAOsJ * nAOsJ);
  if (atomI == atomJ) {
    const auto& integrals = cisData_.oneCenterIntegrals.get(cisData_.elements[atomI]);
    for (int mu = 0; mu < nAOsI; ++mu) {
      for (int nu = 0; nu < nAOsI; ++nu) {
        for (int lambda = 0; lambda < nAOsI; ++lambda) {
          for (int sigma = 0; sigma < nAOsI; ++sigma) {
            integralMatrix(mu * nAOsI + nu, lambda * nAOsJ + sigma) = integrals.get(mu, nu, lambda, sigma);
          }
        }
      }
    }
  }
  else {
    const auto& integrals = cisData_.twoCenterIntegrals.get(atomI, atomJ);
    for (int mu = 0; mu < nAOsI; ++mu) {
      for (int nu = 0; nu < nAOsI; ++nu) {
        for (int lambda = 0; lambda < nAOsJ; ++lambda) {
          for (int sigma = 0; sigma < nAOsJ; ++sigma) {
            integralMatrix(mu * nAOsI + nu, lambda * nAOsJ + sigma) = integrals->get(mu, nu, lambda, sigma);
          }
        }
      }
    }
  }

  return integralMatrix * c1_;
}

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock>
inline Eigen::MatrixXd CISMatrixAOFockBuilder<restrictedness, spinBlock>::getExchangeIntegrals(int atomI, int atomJ) {
  int nAOsI = cisData_.AOInfo.getNOrbitals(atomI);

  int nAOsJ = cisData_.AOInfo.getNOrbitals(atomJ);
  Eigen::MatrixXd integralMatrix = Eigen::MatrixXd::Zero(nAOsI * nAOsJ, nAOsJ * nAOsI);
  if (atomI == atomJ) {
    const auto& integrals = cisData_.oneCenterIntegrals.get(cisData_.elements[atomI]);
    for (int mu = 0; mu < nAOsI; ++mu) {
      for (int nu = 0; nu < nAOsI; ++nu) {
        for (int lambda = 0; lambda < nAOsI; ++lambda) {
          for (int sigma = 0; sigma < nAOsI; ++sigma) {
            integralMatrix(mu * nAOsI + nu, lambda * nAOsI + sigma) = integrals.get(mu, sigma, lambda, nu);
          }
        }
      }
    }
  }
  else {
    const auto& integrals = cisData_.twoCenterIntegrals.get(atomI, atomJ);
    for (int mu = 0; mu < nAOsI; ++mu) {
      for (int nu = 0; nu < nAOsJ; ++nu) {
        for (int lambda = 0; lambda < nAOsJ; ++lambda) {
          for (int sigma = 0; sigma < nAOsI; ++sigma) {
            integralMatrix(mu * nAOsJ + nu, lambda * nAOsI + sigma) = integrals->get(mu, sigma, lambda, nu);
          }
        }
      }
    }
  }
  return integralMatrix * c2_;
}

template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Singlet>::calculate(int atomI,
                                                                                                            int atomJ) {
  coulombContainer_->at(atomI).insert({atomJ, std::make_unique<Eigen::MatrixXd>(2 * getCoulombIntegrals(atomI, atomJ))});
  exchangeContainer_->at(atomI).insert({atomJ, std::make_unique<Eigen::MatrixXd>(getExchangeIntegrals(atomI, atomJ))});
}
template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Triplet>::calculate(int atomI,
                                                                                                            int atomJ) {
  exchangeContainer_->at(atomI).insert({atomJ, std::make_unique<Eigen::MatrixXd>(getExchangeIntegrals(atomI, atomJ))});
}

template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Singlet>::calculate(int atomI,
                                                                                                              int atomJ) {
  coulombContainer_->at(atomI).insert({atomJ, std::make_unique<Eigen::MatrixXd>(getCoulombIntegrals(atomI, atomJ))});
  exchangeContainer_->at(atomI).insert({atomJ, std::make_unique<Eigen::MatrixXd>(getExchangeIntegrals(atomI, atomJ))});
}
template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Triplet>::calculate(int /*atomI*/,
                                                                                                              int /*atomJ*/) {
  throw std::runtime_error(" CISMatrixAOFockBuilder: Invalid combination: Calculation for an unrestricted reference "
                           "with a triplet spin-transition not possible!.");
}

template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Singlet>::calculate(int atomI) {
  coulombContainer_->at(atomI).insert({atomI, std::make_unique<Eigen::MatrixXd>(2 * getCoulombIntegrals(atomI, atomI) -
                                                                                getExchangeIntegrals(atomI, atomI))});
}

template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Triplet>::calculate(int atomI) {
  exchangeContainer_->at(atomI).insert({atomI, std::make_unique<Eigen::MatrixXd>(getExchangeIntegrals(atomI, atomI))});
}
template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Singlet>::calculate(int atomI) {
  coulombContainer_->at(atomI).insert({atomI, std::make_unique<Eigen::MatrixXd>(getCoulombIntegrals(atomI, atomI))});
  exchangeContainer_->at(atomI).insert({atomI, std::make_unique<Eigen::MatrixXd>(getExchangeIntegrals(atomI, atomI))});
}

template<>
inline void CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Triplet>::calculate(int /*atomI*/) {
  throw std::runtime_error(" CISMatrixAOFockBuilder: Invalid combination: Calculation for an unrestricted reference "
                           "with a triplet spin-transition not possible!.");
}

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock>
Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd> CISMatrixAOFockBuilder<restrictedness, spinBlock>::getAOFock(
    const Utils::SpinAdaptedContainer<restrictedness, Eigen::MatrixXd>& pseudoDensity,
    std::map<int, std::vector<int>> atomPairList) const {
  return buildFock(pseudoDensity, std::move(atomPairList));
}

template<>
inline Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd>
CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Singlet>::buildFock(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd>& pseudoDensity,
    std::map<int, std::vector<int>> atomPairList) const {
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd> fockMatrix;
  fockMatrix.restricted = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
  Eigen::MatrixXd tmpBlock(81, 81);
  Eigen::VectorXd tmpV(81 * 81);

  for (int atomI = 0; atomI < nAtoms_; atomI++) {
    int nAOsI = cisData_.AOInfo.getNOrbitals(atomI);
    int startAOI = cisData_.AOInfo.getFirstOrbitalIndex(atomI);
    tmpBlock = pseudoDensity.restricted.block(startAOI, startAOI, nAOsI, nAOsI);
    tmpV = *coulombContainer_->at(atomI)[atomI] * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsI * nAOsI);
    fockMatrix.restricted.block(startAOI, startAOI, nAOsI, nAOsI) +=
        Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

    for (int atomJ : atomPairList.at(atomI)) {
      const Eigen::MatrixXd& coulombMatrixIJ = *coulombContainer_->at(atomI).at(atomJ);
      const Eigen::MatrixXd& exchangeMatrixIJ = *exchangeContainer_->at(atomI).at(atomJ);
      int nAOsJ = cisData_.AOInfo.getNOrbitals(atomJ);
      int startAOJ = cisData_.AOInfo.getFirstOrbitalIndex(atomJ);
      // coulomb
      tmpBlock = pseudoDensity.restricted.block(startAOJ, startAOJ, nAOsJ, nAOsJ);
      tmpV = coulombMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsJ * nAOsJ);
      fockMatrix.restricted.block(startAOI, startAOI, nAOsI, nAOsI) +=
          Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

      // exchange
      tmpBlock = pseudoDensity.restricted.block(startAOI, startAOJ, nAOsI, nAOsJ);
      tmpV = exchangeMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsI * nAOsJ);
      fockMatrix.restricted.block(startAOJ, startAOI, nAOsJ, nAOsI) -=
          Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsI);

      Eigen::MatrixXd transposedMatrix = coulombMatrixIJ.transpose();
      Eigen::MatrixXd transposedMatrixX = exchangeMatrixIJ.transpose();
      // coulomb
      tmpBlock = pseudoDensity.restricted.block(startAOI, startAOI, nAOsI, nAOsI);
      tmpV = transposedMatrix * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsI * nAOsI);
      fockMatrix.restricted.block(startAOJ, startAOJ, nAOsJ, nAOsJ) +=
          Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsJ);
      // exchange
      tmpBlock = pseudoDensity.restricted.block(startAOJ, startAOI, nAOsJ, nAOsI);
      tmpV = transposedMatrixX * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsJ * nAOsI);
      fockMatrix.restricted.block(startAOI, startAOJ, nAOsI, nAOsJ) -=
          Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsJ);
    }
  }
  return fockMatrix;
}
template<>
inline Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd>
CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Triplet>::buildFock(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd>& pseudoDensity,
    std::map<int, std::vector<int>> atomPairList) const {
  Eigen::MatrixXd tmpBlock(81, 81);
  Eigen::VectorXd tmpV(81 * 81);
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd> fockMatrix;
  fockMatrix.restricted = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
  for (int atomI = 0; atomI < nAtoms_; atomI++) {
    int nAOsI = cisData_.AOInfo.getNOrbitals(atomI);
    int startAOI = cisData_.AOInfo.getFirstOrbitalIndex(atomI);
    tmpBlock = pseudoDensity.restricted.block(startAOI, startAOI, nAOsI, nAOsI);
    tmpV = *exchangeContainer_->at(atomI).at(atomI) *
           Eigen::Map<const Eigen::VectorXd>(tmpBlock.transpose().data(), nAOsI * nAOsI);
    fockMatrix.restricted.block(startAOI, startAOI, nAOsI, nAOsI) -=
        Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

    for (int atomJ : atomPairList.at(atomI)) {
      const Eigen::MatrixXd& exchangeMatrixIJ = *exchangeContainer_->at(atomI).at(atomJ);
      int nAOsJ = cisData_.AOInfo.getNOrbitals(atomJ);
      int startAOJ = cisData_.AOInfo.getFirstOrbitalIndex(atomJ);
      tmpBlock = pseudoDensity.restricted.block(startAOI, startAOJ, nAOsI, nAOsJ);
      tmpV = exchangeMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsI * nAOsJ);
      fockMatrix.restricted.block(startAOJ, startAOI, nAOsJ, nAOsI) -=
          Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsI);
      tmpBlock = pseudoDensity.restricted.block(startAOJ, startAOI, nAOsJ, nAOsI);
      tmpV = exchangeMatrixIJ.transpose() * Eigen::Map<const Eigen::VectorXd>(tmpBlock.data(), nAOsJ * nAOsI);
      fockMatrix.restricted.block(startAOI, startAOJ, nAOsI, nAOsJ) -=
          Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsJ);
    }
  }
  return fockMatrix;
}
template<>
inline Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd>
CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Singlet>::buildFock(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd>& pseudoDensity,
    std::map<int, std::vector<int>> atomPairList) const {
  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd> fockMatrix;
  fockMatrix.alpha = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
  fockMatrix.beta = Eigen::MatrixXd::Zero(nAOs_, nAOs_);
  Eigen::MatrixXd tmpBlockAlpha;
  Eigen::MatrixXd tmpBlockBeta;
  Eigen::VectorXd tmpV;
  for (int atomI = 0; atomI < nAtoms_; atomI++) {
    int nAOsI = cisData_.AOInfo.getNOrbitals(atomI);
    int startAOI = cisData_.AOInfo.getFirstOrbitalIndex(atomI);

    tmpBlockAlpha = pseudoDensity.alpha.block(startAOI, startAOI, nAOsI, nAOsI);
    tmpV = *coulombContainer_->at(atomI).at(atomI) * Eigen::Map<const Eigen::VectorXd>(tmpBlockAlpha.data(), nAOsI * nAOsI);

    tmpBlockBeta = pseudoDensity.beta.block(startAOI, startAOI, nAOsI, nAOsI);
    tmpV.noalias() +=
        *coulombContainer_->at(atomI).at(atomI) * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsI * nAOsI);
    fockMatrix.alpha.block(startAOI, startAOI, nAOsI, nAOsI) += Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);
    fockMatrix.beta.block(startAOI, startAOI, nAOsI, nAOsI) += Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

    tmpV = *exchangeContainer_->at(atomI).at(atomI) * Eigen::Map<const Eigen::VectorXd>(tmpBlockAlpha.data(), nAOsI * nAOsI);
    fockMatrix.alpha.block(startAOI, startAOI, nAOsI, nAOsI) -= Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

    tmpV = *exchangeContainer_->at(atomI).at(atomI) * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsI * nAOsI);
    fockMatrix.beta.block(startAOI, startAOI, nAOsI, nAOsI) -= Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

    for (int atomJ : atomPairList.at(atomI)) {
      const Eigen::MatrixXd& coulombMatrixIJ = *coulombContainer_->at(atomI).at(atomJ);
      const Eigen::MatrixXd& exchangeMatrixIJ = *exchangeContainer_->at(atomI).at(atomJ);
      int nAOsJ = cisData_.AOInfo.getNOrbitals(atomJ);
      int startAOJ = cisData_.AOInfo.getFirstOrbitalIndex(atomJ);
      // coulomb alpha + beta
      tmpBlockAlpha = pseudoDensity.alpha.block(startAOI, startAOI, nAOsI, nAOsI);
      tmpV = coulombMatrixIJ.transpose() * Eigen::Map<const Eigen::VectorXd>(tmpBlockAlpha.data(), nAOsI * nAOsI);

      tmpBlockBeta = pseudoDensity.beta.block(startAOI, startAOI, nAOsI, nAOsI);
      tmpV.noalias() += coulombMatrixIJ.transpose() * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsI * nAOsI);
      fockMatrix.alpha.block(startAOJ, startAOJ, nAOsJ, nAOsJ) += Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsJ);
      fockMatrix.beta.block(startAOJ, startAOJ, nAOsJ, nAOsJ) += Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsJ);

      tmpBlockAlpha = pseudoDensity.alpha.block(startAOJ, startAOJ, nAOsJ, nAOsJ);
      tmpV = coulombMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlockAlpha.data(), nAOsJ * nAOsJ);

      tmpBlockBeta = pseudoDensity.beta.block(startAOJ, startAOJ, nAOsJ, nAOsJ);
      tmpV.noalias() += coulombMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsJ * nAOsJ);
      fockMatrix.alpha.block(startAOI, startAOI, nAOsI, nAOsI) += Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);
      fockMatrix.beta.block(startAOI, startAOI, nAOsI, nAOsI) += Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsI);

      // exchange alpha
      tmpBlockAlpha = pseudoDensity.alpha.block(startAOI, startAOJ, nAOsI, nAOsJ);
      tmpV = exchangeMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlockAlpha.data(), nAOsI * nAOsJ);
      fockMatrix.alpha.block(startAOJ, startAOI, nAOsJ, nAOsI) -= Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsI);
      tmpBlockBeta = pseudoDensity.alpha.block(startAOJ, startAOI, nAOsJ, nAOsI);
      tmpV = exchangeMatrixIJ.transpose() * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsJ * nAOsI);
      fockMatrix.alpha.block(startAOI, startAOJ, nAOsI, nAOsJ) -= Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsJ);
      // exchange beta

      tmpBlockBeta = pseudoDensity.beta.block(startAOI, startAOJ, nAOsI, nAOsJ);
      tmpV = exchangeMatrixIJ * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsI * nAOsJ);
      fockMatrix.beta.block(startAOJ, startAOI, nAOsJ, nAOsI) -= Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsJ, nAOsI);
      tmpBlockBeta = pseudoDensity.beta.block(startAOJ, startAOI, nAOsJ, nAOsI);
      tmpV = exchangeMatrixIJ.transpose() * Eigen::Map<const Eigen::VectorXd>(tmpBlockBeta.data(), nAOsJ * nAOsI);
      fockMatrix.beta.block(startAOI, startAOJ, nAOsI, nAOsJ) -= Eigen::Map<const Eigen::MatrixXd>(tmpV.data(), nAOsI, nAOsJ);
    }
  }
  return fockMatrix;
}

template<>
Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd>
CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Triplet>::buildFock(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd>& /*pseudoDensity*/,
    std::map<int, std::vector<int>> /*atomPairList*/) const {
  throw std::runtime_error(" CISMatrixAOFockBuilder: Invalid combination: Calculation for an unrestricted reference "
                           "with a triplet spin-transition not possible!.");
}

template<Utils::Reference restrictedness, Utils::SpinTransition spinBlock>
CISMatrixAOFockBuilder<restrictedness, spinBlock>::~CISMatrixAOFockBuilder() = default;
template class CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Singlet>;
template class CISMatrixAOFockBuilder<Utils::Reference::Restricted, Utils::SpinTransition::Triplet>;
template class CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Singlet>;
template class CISMatrixAOFockBuilder<Utils::Reference::Unrestricted, Utils::SpinTransition::Triplet>;
} // namespace Sparrow
} // namespace Scine
