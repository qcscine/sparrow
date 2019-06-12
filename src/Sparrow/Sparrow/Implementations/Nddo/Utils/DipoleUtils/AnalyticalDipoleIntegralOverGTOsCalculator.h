/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ANALYTICALDIPOLEINTEGRALOVERGTOSCALCULATOR_H
#define SPARROW_ANALYTICALDIPOLEINTEGRALOVERGTOSCALCULATOR_H

#include <Eigen/Core>
#include <array>

namespace Scine {
namespace Sparrow {

class AnalyticalDipoleIntegralOverGTOsCalculator {
 public:
  AnalyticalDipoleIntegralOverGTOsCalculator(int angularMomentumA, int angulatMomentumB, double expA, double expB,
                                             const Eigen::Vector3d& Ra, const Eigen::Vector3d& Rb,
                                             const Eigen::Vector3d& evaluationCoordinate);
  std::array<double, 3> calculateAnalyticalDipoleElement();

 private:
  const double exponentA_, exponentB_;
  double expSum_;
  const int angularMomentumA_, angularMomentumB_;
  const Eigen::Vector3d Ra_;
  const Eigen::Vector3d Rb_;
  const Eigen::Vector3d evaluationCoordinate_;
  Eigen::Vector3d Rab_, weightedSum_;
  double dipoleSS(int dimension);
  double dipoleSP(double expA, double Ra, double Rb, int dimension);
  double dipoleSD(double expA, double expB, double Ra, double Rb, double evaluationCoordinate);
  double dipolePP(double expA, double expB, double Ra, double Rb, double evaluationCoordinate);
  double dipolePD(double expA, double expB, double Ra, double Rb, double evaluationCoordinate);
  double dipoleDD(double expA, double expB, double Ra, double Rb, double evaluationCoordinate);
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_ANALYTICALDIPOLEINTEGRALOVERGTOSCALCULATOR_H
