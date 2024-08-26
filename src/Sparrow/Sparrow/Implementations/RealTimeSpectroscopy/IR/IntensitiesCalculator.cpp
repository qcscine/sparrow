/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#define _USE_MATH_DEFINES
#include "IntensitiesCalculator.h"
#include "../Utils/LineWidthGenerator.h"
#include <cmath>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {
constexpr double cSquared = 1.6021765 * 1.6021765e-38;
constexpr double coulombForceConstant = 8.9875517873681e9;
constexpr double avogadrosNumber = 6.02214199e23;
constexpr double c = 299792458.;
constexpr double atomicMassUnit = 1.6605387280149467e-27;
constexpr double conversionFactor =
    cSquared * M_PI * coulombForceConstant * avogadrosNumber / (3. * c * c) / atomicMassUnit / 1000.;

const Eigen::VectorXd IntensitiesCalculator::getAdsorptionCoefficients() const {
  return squaredNormalDipoleGradient_ * conversionFactor;
}

Eigen::VectorXd IntensitiesCalculator::calculateRelativeVibrationalIntensities() {
  Eigen::VectorXd relativeIntensities = squaredNormalDipoleGradient_;
  relativeIntensities /= relativeIntensities.maxCoeff();
  return relativeIntensities;
}

void IntensitiesCalculator::setSquaredNormalDipoleGradient(const Eigen::VectorXd& newSquaredNormalDipoleGradient) {
  squaredNormalDipoleGradient_ = newSquaredNormalDipoleGradient;
}

Eigen::VectorXd
IntensitiesCalculator::transformCartesianToSquaredNormalDipoleGradient(const Eigen::MatrixXd& massWeightedEigenvectors,
                                                                       const Eigen::MatrixX3d& cartesianDipoleMoment,
                                                                       const std::vector<double>& masses) {
  int nAtoms = masses.size();
  Eigen::Matrix3Xd massWeightedDipoleGradient = cartesianDipoleMoment.transpose();
  for (int atom = 0; atom < nAtoms; ++atom) {
    massWeightedDipoleGradient.middleCols(3 * atom, 3) /= std::sqrt(masses[atom]);
  }
  // get mass weighted eigenvector matrix A and calculate the three components of du/dQ (ux, uy, uz) as
  // (du/dQ)' = (du/dR)' * A = (A' * du/dR)'. In eigen Matrix3Xd is a 3x3N matrix representing (du/dR)' (' =
  // transposed). calculates the squared norm
  Eigen::VectorXd normalDipoleGradient = (massWeightedDipoleGradient * massWeightedEigenvectors).colwise().squaredNorm();

  return normalDipoleGradient;
}

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
