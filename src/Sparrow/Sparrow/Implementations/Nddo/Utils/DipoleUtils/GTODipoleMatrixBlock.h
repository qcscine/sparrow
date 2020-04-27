/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_GTODIPOLEMATRIXBLOCK_H
#define SPARROW_GTODIPOLEMATRIXBLOCK_H

#include <Sparrow/Implementations/DipoleMatrixCalculator.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GTOOverlapMatrixBlock.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <array>

namespace Scine {

namespace Utils {
class GtoExpansion;
} // namespace Utils

namespace Sparrow {

/**
 * @brief Class responsible for calculating a block of the dipole matrix in GTO ao basis.
 * The basis must be one of the STO-nG type.
 */
class GTODipoleMatrixBlock {
 public:
  using oneDimensionalArray = std::array<double, 3>;
  using twoDimensionalArray = std::array<oneDimensionalArray, 3>;
  using threeDimensionalArray = std::array<twoDimensionalArray, 3>;

  /**
   * @brief Constructor. It initializes the angular momenta and the orbital AO indices.
   */
  GTODipoleMatrixBlock();

  /**
   * Sets the integral method to either Obara-Saika or analytical.
   */
  void setIntegralMethod(IntegralMethod method);
  /**
   * @brief Initializes the data structures used throughout the calculation.
   */
  void initialize(const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                  const Eigen::RowVector3d& dipoleEvaluationCoordinate, double expA, double expB,
                  const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB);

  /**
   * @brief Calculates a block of the dipole matrix.
   */
  std::array<Eigen::MatrixXd, 3> createSTOBlock(const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                                const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                                                const Eigen::Vector3d& Rab, const Eigen::Vector3d& dipoleEvaluationCoordinate);

 private:
  void dOrbitalsFromSixCartesianToFiveRealSolidHarmonics();
  static double getNormalizationFactorAndCoefficient(int GaussTypeFunctionIndexA, int GaussTypeFunctionIndexB,
                                                     const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                                     double expCoefficientS_00, const Eigen::Vector3d& Rab);
  void calculateAnalyticalDipoleIntegrals(int GaussTypeFunctionIndexA, int GaussTypeFunctionIndexB,
                                          const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                          const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                                          const Eigen::Vector3d& dipoleEvaluationCoordinate);
  void calculateSingleGTFPair(int GaussTypeFunctionIndexA, int GaussTypeFunctionIndexB, const Utils::GtoExpansion& gtoA,
                              const Utils::GtoExpansion& gtoB, const Eigen::RowVector3d& Ra, const Eigen::RowVector3d& Rb,
                              const Eigen::Vector3d& Rab, const Eigen::Vector3d& dipoleEvaluationCoordinate);
  void calculateSingleGTFPairIntegralOverShell(const Utils::GtoExpansion& gtoA, const Utils::GtoExpansion& gtoB,
                                               double exponentSum);
  void createBlockForOneGTFOverShell(double normalizedCoefficientPrefactor);

  int startAOonA_;
  int startAOonB_;
  int orbitalShellSizeA_;
  int orbitalShellSizeB_;
  std::array<Eigen::Matrix<double, 6, 6>, 3> dipoleComponentsBlocks_;
  std::array<nddo::AngularMomentum, 10> AOMomenta_;
  std::array<int, 10> AOIndexes_;
  oneDimensionalArray PminusA_;
  oneDimensionalArray PminusB_;
  oneDimensionalArray PminusC_;
  threeDimensionalArray overlapIntegral_;
  threeDimensionalArray dipoleIntegral_;
  IntegralMethod method_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_GTODIPOLEMATRIXBLOCK_H
