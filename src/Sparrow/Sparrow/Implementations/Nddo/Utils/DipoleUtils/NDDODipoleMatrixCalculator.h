/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DIPOLEMATRTIXCALCULATOR_H
#define SPARROW_DIPOLEMATRTIXCALCULATOR_H

#include "Sparrow/Implementations/Nddo/Utils/DipoleUtils/GTODipoleMatrixBlock.h"
#include <Sparrow/Implementations/DipoleMatrixCalculator.h>
#include <Utils/DataStructures/DipoleMatrix.h>
#include <Eigen/Core>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
class AtomsOrbitalsIndexes;
class MolecularOrbitals;
} // namespace Utils

namespace Sparrow {
namespace nddo {
class ElementParameters;
} // namespace nddo

/**
 * @brief Class responsible for the calculation of the dipole matrix.
 * @tparam NDDOMethod An NDDO method, must derive from Utils::ScfMethod
 *         and must have the method getInitializer().
 * getInitializer must return a valid NDDOInitializer instance.
 * This class is responsible for the calculation of the dipole matrix in both
 * atomic and molecular orbital basis. The dipole matrix in atomic orbital basis
 * is calculated by explicitly integrating in closed form or using the highly
 * efficient Obara-Saika scheme the dipole integral <\mu|r|\nu>.
 * The dipole matrix in atomic orbital basis is then transformed in molecular
 * orbital basis by D_{MO} = C^T * D_{AO} * C in restricted formalism,
 * D_{MO} = C_\alpha^T * D_{AO} * C_\alpha + C_\beta^T * D_{AO} * C_\beta in
 * unrestricted formalism.
 */
template<class NDDOMethod>
class NDDODipoleMatrixCalculator : public DipoleMatrixCalculator {
 public:
  static std::unique_ptr<NDDODipoleMatrixCalculator<NDDOMethod>> create(NDDOMethod& method);

  ~NDDODipoleMatrixCalculator() final;
  /**
   * @brief Getter for the dipole matrix in AO.
   */
  const Utils::DipoleMatrix& getAODipoleMatrix() const final;
  /**
   * @brief Getter for the dipole matrix in MO.
   * The dipole matrix in MO basis is transformed from the one in AO basis with every call to this function.
   * @return The dipole matrix in MO basis.
   */
  Utils::DipoleMatrix getMODipoleMatrix() const final;
  /**
   * @brief Setter for the dipole matrix in AO basis.
   */
  void setAODipoleMatrix(Utils::DipoleMatrix dipoleMatrix) final;
  /**
   * @brief Calculates the dipole matrix in AO basis.
   * @param dipoleEvaluationCoordinate The coordinates to consider as the origin for the dipole calculation.
   */
  void fillDipoleMatrix(const Eigen::RowVector3d& dipoleEvaluationCoordinate) final;
  /**
   * @brief Initializes the matrices and invalidates the current dipole matrix.
   */
  void initialize() final;
  /**
   * @brief Sets the dipole integral calculation method.
   * @param IntegralMethod Either calculates is from the evaluation of a closed form formula, or with the
   *        Obara-Saika recursive algorithm.
   */
  void setIntegralMethod(const IntegralMethod& IntegralMethod) final;
  /**
   * @brief Return the validity status of the dipole matrix.
   */
  bool isValid() const final;
  /**
   * @brief Invalidates the underlying dipole matrices and forces a new calculation.
   */
  void invalidate() final;

 private:
  explicit NDDODipoleMatrixCalculator(NDDOMethod& method);
  Utils::DipoleMatrix calculateMODipoleMatrixRestricted() const;
  Utils::DipoleMatrix calculateMODipoleMatrixUnrestricted() const;
  int nAOs_;
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;
  const Utils::ElementTypeCollection& elementTypes_;
  const Utils::PositionCollection& positions_;
  const nddo::ElementParameters& elementParameters_;
  const Eigen::MatrixXd& overlapMatrix_;
  const Utils::MolecularOrbitals& molecularOrbitals_;
  Utils::DipoleMatrix dipoleMatrix_;
  int nAtoms_;
  IntegralMethod integralMethod_{IntegralMethod::ObaraSaika};
  bool valid_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DIPOLEMATRTIXCALCULATOR_H
