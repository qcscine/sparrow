/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_NDDODIPOLEMOMENTCALCULATOR_H
#define SPARROW_NDDODIPOLEMOMENTCALCULATOR_H

#include <Sparrow/Implementations/DipoleMomentCalculator.h>
#include <Utils/DataStructures/DipoleMatrix.h>
#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {
namespace Sparrow {
class DipoleMatrixCalculator;
/**
 * @brief Class resposible for the calculation of the dipole in the NDDO methods.
 * It must be able to calculate the dipole both with the NDDO approximation and with the use of the dipole matrix.
 * @tparam NDDOMethod An NDDO method, must derive from Utils::ScfMethod and must have the method getInitializer.
 *         getInitializer must return a valid NDDOInitializer instance.
 */
template<class NDDOMethod>
class NDDODipoleMomentCalculator : public DipoleMomentCalculator {
 public:
  /**
   * @brief Factory method for the NDDODipoleMomentCalculator class.
   * @param method An NDDO method, i.e. PM6Method, MNDOMethod, AM1Method,...
   * @param dipoleMatrixCalculator An instance to the class calculating the dipole matrix.
   * @return An unique pointer to an instance of NDDODipoleMomentCalculator<NDDOMethod>
   */
  static std::unique_ptr<NDDODipoleMomentCalculator<NDDOMethod>> create(NDDOMethod& method,
                                                                        DipoleMatrixCalculator& dipoleMatrixCalculator);
  ~NDDODipoleMomentCalculator() final;

  /**
   * @brief Calculates the molecular electrical dipole moment.
   * @return An Eigen::Vector3d containing the dipole moment components in atomic units.
   */
  Eigen::RowVector3d calculate() const final;

  /**
   * @brief Sets wether to use the NDDO dipole approximation or calculate the dipole from the dipole matrix.
   * @param useNDDOApprox The flag present in the method wrapper settings.
   */
  void useNDDOApproximation(bool useNDDOApprox);

 private:
  // Implementation of the calculate method.
  Eigen::RowVector3d calculateWithNDDOApproximation(std::vector<double> atomicCharges, Utils::PositionCollection positions,
                                                    Eigen::MatrixXd densityMatrix, Utils::ElementTypeCollection elements,
                                                    std::vector<int> nrAOs, std::vector<double> chargeSeparationSP,
                                                    std::vector<double> chargeSeparationPD) const;
  Eigen::RowVector3d calculateWithDipoleMatrix(std::vector<double> coreCharges, Utils::PositionCollection positions,
                                               Eigen::MatrixXd densityMatrix, Utils::DipoleMatrix dipoleMatrix,
                                               Eigen::MatrixXd overlapMatrix,
                                               Eigen::RowVector3d dipoleEvaluationCoordinate) const;
  // Private constructor.
  NDDODipoleMomentCalculator(NDDOMethod& method, DipoleMatrixCalculator& dipoleMatrixCalculator);
  ;
  NDDOMethod& method_;
  DipoleMatrixCalculator& dipoleMatrixCalculator_;

  bool useNDDOApproximation_{true};
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_NDDODIPOLEMOMENTCALCULATOR_H
