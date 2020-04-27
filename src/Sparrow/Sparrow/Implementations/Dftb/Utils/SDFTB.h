/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SSPARROW_DFTB_H
#define SSPARROW_DFTB_H

#include "Utils/DataStructures/AtomsOrbitalsIndexes.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Utils/Math/DerivOrderEnum.h>
#include <Eigen/Core>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
class MatrixWithDerivatives;
class SpinAdaptedMatrix;
} // namespace Utils
namespace Sparrow {

namespace dftb {
class SKAtom;
class DFTBCommon;

class SDFTB {
 public:
  explicit SDFTB(const Utils::ElementTypeCollection& elements, const std::vector<std::unique_ptr<SKAtom>>& atomicParameters);

  ~SDFTB() = default;
  void spinPopulationAnalysis(const Eigen::MatrixXd& densityMatrixUp, const Eigen::MatrixXd& densityMatrixDn,
                              const Eigen::MatrixXd& overlapMatrix);

  void initialize(int nAtoms, int nAOs, Utils::AtomsOrbitalsIndexes indexes);
  void calculateSpinContribution();
  void constructSpinHamiltonians(Utils::SpinAdaptedMatrix& H, const Eigen::MatrixXd& overlap) const;
  double spinEnergyContribution() const;
  template<Utils::derivativeType O>
  void addDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivativesContainer,
                      const Utils::MatrixWithDerivatives& overlapDeriv, const Eigen::MatrixXd& pUp,
                      const Eigen::MatrixXd& pDn) const;

 private:
  void spinPopulationAnalysis(const Eigen::MatrixXd& densityMatrix, const Eigen::MatrixXd& overlapMatrix,
                              std::vector<double>& population);

  int nAtoms_;
  int nAOs_;
  std::vector<double> pup, // spinup P vector (in case of spin-polarized calculation)
      pdn,                 // spindown P vector (in case of spin-polarized calculation)
      pdif;                // difference P vector (in case of spin-polarized calculation)
  const std::vector<std::unique_ptr<SKAtom>>& atomParameters; // parameters for atoms
  Utils::AtomsOrbitalsIndexes aoIndexes_;
  Eigen::MatrixXd spinContribution_;
  const Utils::ElementTypeCollection& elementTypes_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SSPARROW_DFTB_H
