/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "../../../TimeDependent/TimeDependentUtils.h"
#include <Eigen/Core>
namespace Scine {
namespace Sparrow {
#ifndef SPARROW_CISSPINCONTAMINATOR_H
#  define SPARROW_CISSPINCONTAMINATOR_H

class CISSpinContaminator {
 public:
  static Eigen::VectorXd calculateSpinContaminationOpenShell(
      const Utils::MolecularOrbitals& mos, const Eigen::MatrixXd& eigenVectors, const std::vector<int>& filledAlpha,
      const std::vector<int>& filledBeta,
      const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, std::vector<int>>& excitationIndices);
  static double ab_j_iAlpha(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& alphaCoeffs,
                            const std::vector<int>& virtualAlpha, const std::vector<int>& occupiedBeta,
                            const std::vector<int>& occupiedAlpha);
  static double ab_j_iBeta(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& betaCoeffs,
                           const std::vector<int>& virtualBeta, const std::vector<int>& occupiedAlpha,
                           const std::vector<int>& occupiedBeta);
  static double ij_k_aAlpha(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& alphaCoeffs,
                            const std::vector<int>& occupiedAlpha, const std::vector<int>& virtualAlpha,
                            const std::vector<int>& occupiedBeta);
  static double ij_k_aBeta(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& betaCoeffs,
                           const std::vector<int>& occupiedBeta, const std::vector<int>& virtualBeta,
                           const std::vector<int>& occupiedAlpha);
  static double ijab(const Eigen::MatrixXd& spatialOverlap, const Eigen::VectorXd& alphaCoeffs, const Eigen::VectorXd& betaCoeffs,
                     const std::vector<int>& occupiedAlpha, const std::vector<int>& virtualAlpha,
                     const std::vector<int>& occupiedBeta, const std::vector<int>& virtualBeta);
};
}
}

#endif // SPARROW_CISSPINCONTAMINATOR_H
