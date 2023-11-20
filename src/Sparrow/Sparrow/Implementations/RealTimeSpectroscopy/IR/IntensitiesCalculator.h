/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef RTSPECTROSCOPY_INTENSITIESCALCULATOR_H
#define RTSPECTROSCOPY_INTENSITIESCALCULATOR_H

#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <vector>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

class IntensitiesCalculator {
 public:
  const Eigen::VectorXd getAdsorptionCoefficients() const;
  void setSquaredNormalDipoleGradient(const Eigen::VectorXd& newSquaredNormalDipoleGradient);
  static Eigen::VectorXd transformCartesianToSquaredNormalDipoleGradient(const Eigen::MatrixXd& massWeightedEigenvectors,
                                                                         const Eigen::MatrixX3d& cartesianDipoleMoment,
                                                                         const std::vector<double>& masses);

  Eigen::VectorXd calculateRelativeVibrationalIntensities();

 private:
  Eigen::VectorXd squaredNormalDipoleGradient_;
};

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // RTSPECTROSCOPY_INTENSITIESCALCULATOR_H
