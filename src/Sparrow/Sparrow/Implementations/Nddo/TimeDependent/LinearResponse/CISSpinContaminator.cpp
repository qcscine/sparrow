/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CISSpinContaminator.h"
namespace Scine {
namespace Sparrow {
Eigen::VectorXd CISSpinContaminator::calculateSpinContaminationOpenShell(
    const Utils::MolecularOrbitals& mos, const Eigen::MatrixXd& eigenVectors, const std::vector<int>& filledAlpha,
    const std::vector<int>& filledBeta,
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<int>>& excitationIndices) {
  if (eigenVectors.rows() != static_cast<Eigen::Index>(excitationIndices.alpha.size() + excitationIndices.beta.size())) {
    throw std::runtime_error("Wrong number of excitations in the calculation of the spin contamination.");
  }
  int eigenvalues = eigenVectors.cols();
  Eigen::VectorXd spinContamination = Eigen::VectorXd::Ones(eigenVectors.cols());
  const Eigen::MatrixXd spatialOverlap = mos.alphaMatrix().transpose() * mos.betaMatrix();
  int nAlphaElectrons = filledAlpha.size();
  int nBetaElectrons = filledBeta.size();

  std::vector<int> virtualAlpha(spatialOverlap.rows() - nAlphaElectrons);
  std::vector<int> virtualBeta(spatialOverlap.cols() - nBetaElectrons);

  TimeDependentUtils::generateVirtualOrbitalIndices(filledAlpha, virtualAlpha);
  TimeDependentUtils::generateVirtualOrbitalIndices(filledBeta, virtualBeta);

  auto nAlphaExc = nAlphaElectrons * virtualAlpha.size();
  auto nBetaExc = nBetaElectrons * virtualBeta.size();
  Eigen::MatrixXd sparseEigenVectorsAlpha = Eigen::MatrixXd::Zero(nAlphaExc, eigenvalues);
  Eigen::MatrixXd sparseEigenVectorsBeta = Eigen::MatrixXd::Zero(nBetaExc, eigenvalues);
  if (excitationIndices.alpha.size() == nAlphaExc && excitationIndices.beta.size() == nBetaExc) {
    sparseEigenVectorsAlpha = eigenVectors.block(0, 0, nAlphaExc, eigenvalues);
    sparseEigenVectorsBeta = eigenVectors.block(nAlphaExc, 0, nBetaExc, eigenvalues);
  }
  else {
    int iter = 0;
    for (int exc : excitationIndices.alpha) {
      sparseEigenVectorsAlpha.row(exc) = eigenVectors.row(iter);
      iter++;
    }
    for (int exc : excitationIndices.beta) {
      sparseEigenVectorsBeta.row(exc) = eigenVectors.row(iter);
      iter++;
    }
  }
  double nDiff = static_cast<double>(nAlphaElectrons - nBetaElectrons) / 2.;
  double sUHF = (nDiff * (nDiff + 1.)) + static_cast<double>(nBetaElectrons);

  for (int occI : filledAlpha) {
    for (int occJ : filledBeta) {
      sUHF -= std::pow(spatialOverlap(occI, occJ), 2);
    }
  }
  for (int eV = 0; eV < eigenvalues; ++eV) {
    const Eigen::VectorXd& alphaCoeffs = sparseEigenVectorsAlpha.col(eV);
    const Eigen::VectorXd& betaCoeffs = sparseEigenVectorsBeta.col(eV);
    double sCeV = sUHF;
    sCeV -= ab_j_iAlpha(spatialOverlap, alphaCoeffs, virtualAlpha, filledBeta, filledAlpha);
    sCeV -= ab_j_iBeta(spatialOverlap, betaCoeffs, virtualBeta, filledAlpha, filledBeta);
    sCeV -= ij_k_aAlpha(spatialOverlap, alphaCoeffs, filledAlpha, virtualAlpha, filledBeta);
    sCeV -= ij_k_aBeta(spatialOverlap, betaCoeffs, filledBeta, virtualBeta, filledAlpha);
    sCeV -= ijab(spatialOverlap, alphaCoeffs, betaCoeffs, filledAlpha, virtualAlpha, filledBeta, virtualBeta);
    spinContamination[eV] = sCeV;
  }

  return spinContamination;
}

double CISSpinContaminator::ab_j_iAlpha(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& alphaCoeffs,
                                        const std::vector<int>& virtualAlpha, const std::vector<int>& occupiedBeta,
                                        const std::vector<int>& occupiedAlpha) {
  double contribution = 0.;
  int alphaVir = virtualAlpha.size();
  int alphaOcc = occupiedAlpha.size();
  for (int a = 0; a < alphaVir; a++) {
    int virA = virtualAlpha[a];
    for (int b = 0; b < alphaVir; b++) {
      int virB = virtualAlpha[b];
      double spatial = 0.;
      for (int occJ : occupiedBeta) {
        spatial += spatialOverlap(virA, occJ) * spatialOverlap(virB, occJ);
      }
      double coeffs = 0.;
      for (int i = 0; i < alphaOcc; i++) {
        coeffs += alphaCoeffs(i * alphaVir + a) * alphaCoeffs(i * alphaVir + b);
      }
      contribution += coeffs * spatial;
    }
  }
  return contribution;
}

double CISSpinContaminator::ab_j_iBeta(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& betaCoeffs,
                                       const std::vector<int>& virtualBeta, const std::vector<int>& occupiedAlpha,
                                       const std::vector<int>& occupiedBeta) {
  double contribution = 0.;
  int betaVir = virtualBeta.size();
  int betOcc = occupiedBeta.size();
  for (int a = 0; a < betaVir; a++) {
    int virA = virtualBeta[a];
    for (int b = 0; b < betaVir; b++) {
      int virB = virtualBeta[b];
      double spatial = 0.;
      for (int occJ : occupiedAlpha) {
        spatial += spatialOverlap(occJ, virA) * spatialOverlap(occJ, virB);
      }
      double coeffs = 0.;
      for (int i = 0; i < betOcc; i++) {
        coeffs += betaCoeffs(i * betaVir + a) * betaCoeffs(i * betaVir + a);
      }
      contribution += coeffs * spatial;
    }
  }
  return contribution;
}

double CISSpinContaminator::ij_k_aAlpha(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& alphaCoeffs,
                                        const std::vector<int>& occupiedAlpha, const std::vector<int>& virtualAlpha,
                                        const std::vector<int>& occupiedBeta) {
  double contribution = 0.;
  int alphaVir = virtualAlpha.size();
  int alphaOcc = occupiedAlpha.size();
  for (int i = 0; i < alphaOcc; i++) {
    int occI = occupiedAlpha[i];
    for (int j = 0; j < alphaOcc; j++) {
      int occJ = occupiedAlpha[j];
      double spatial = 0.;
      for (int occK : occupiedBeta) {
        spatial += spatialOverlap(occI, occK) * spatialOverlap(occJ, occK);
      }
      double coeffs = 0.;
      for (int a = 0; a < alphaVir; a++) {
        coeffs -= alphaCoeffs(i * alphaVir + a) * alphaCoeffs(j * alphaVir + a);
      }
      contribution += coeffs * spatial;
    }
  }
  return contribution;
}

double CISSpinContaminator::ij_k_aBeta(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& betaCoeffs,
                                       const std::vector<int>& occupiedBeta, const std::vector<int>& virtualBeta,
                                       const std::vector<int>& occupiedAlpha) {
  double contribution = 0.;
  int betaVir = virtualBeta.size();
  int betaOcc = occupiedBeta.size();
  for (int i = 0; i < betaOcc; i++) {
    int occI = occupiedBeta[i];
    for (int j = 0; j < betaOcc; j++) {
      int occJ = occupiedBeta[j];
      double spatial = 0.;
      for (int occK : occupiedAlpha) {
        spatial += spatialOverlap(occI, occK) * spatialOverlap(occJ, occK);
      }
      double coeffs = 0.;
      for (int a = 0; a < betaVir; a++) {
        coeffs -= betaCoeffs(i * betaVir + a) * betaCoeffs(j * betaVir + a);
      }
      contribution += coeffs * spatial;
    }
  }
  return contribution;
}

double CISSpinContaminator::ijab(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& alphaCoeffs,
                                 const Eigen::VectorXd& betaCoeffs, const std::vector<int>& occupiedAlpha,
                                 const std::vector<int>& virtualAlpha, const std::vector<int>& occupiedBeta,
                                 const std::vector<int>& virtualBeta) {
  double contribution = 0.;
  int alphaVir = virtualAlpha.size();
  int alphaOcc = occupiedAlpha.size();
  int betaOcc = occupiedBeta.size();
  int betaVir = virtualBeta.size();

  for (int i = 0; i < alphaOcc; i++) {
    int occI = occupiedAlpha[i];
    for (int j = 0; j < betaOcc; j++) {
      int occJ = occupiedBeta[j];
      for (int a = 0; a < alphaVir; a++) {
        int virA = occupiedAlpha[a];
        for (int b = 0; b < betaVir; b++) {
          int virB = virtualBeta[b];
          contribution += spatialOverlap(occI, occJ) * spatialOverlap(virA, virB) * alphaCoeffs(i * alphaVir + a) *
                          betaCoeffs(j * betaVir + b);
        }
      }
    }
  }
  return 2 * contribution;
}
} // namespace Sparrow
} // namespace Scine
