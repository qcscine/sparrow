/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "AnalyticalDipoleIntegralOverGTOsCalculator.h"

namespace Scine {
namespace Sparrow {

AnalyticalDipoleIntegralOverGTOsCalculator::AnalyticalDipoleIntegralOverGTOsCalculator(
    int angularMomentumA, int angularMomentumB, double expA, double expB, const Eigen::Vector3d& Ra,
    const Eigen::Vector3d& Rb, const Eigen::Vector3d& evaluationCoordinate)
  : exponentA_(expA),
    exponentB_(expB),
    expSum_(expA + expB),
    angularMomentumA_(angularMomentumA),
    angularMomentumB_(angularMomentumB),
    Ra_(Ra),
    Rb_(Rb),
    evaluationCoordinate_(evaluationCoordinate) {
  weightedSum_.resize(3);
  for (int dimension = 0; dimension < 3; ++dimension)
    weightedSum_(dimension) = expA * Ra_(dimension) + expB * Rb_(dimension);
  Rab_ = Ra_ - Rb_;
}

std::array<double, 3> AnalyticalDipoleIntegralOverGTOsCalculator::calculateAnalyticalDipoleElement() {
  assert(angularMomentumA_ <= 2 && angularMomentumA_ >= 0 && "Angular momentum is bigger than 2 or smaller than 0!");
  assert(angularMomentumB_ <= 2 && angularMomentumB_ >= 0 && "Angular momentum is bigger than 2 or smaller than 0!");

  std::array<double, 3> analyticalDipole{};
  for (int dimension = 0; dimension < 3; ++dimension) {
    if (angularMomentumA_ == 0 && angularMomentumB_ == 0)
      analyticalDipole[dimension] = dipoleSS(dimension);
    else if (angularMomentumA_ == 0 && angularMomentumB_ == 1)
      analyticalDipole[dimension] = dipoleSP(exponentA_, Ra_[dimension], Rb_[dimension], dimension);
    else if (angularMomentumA_ == 1 && angularMomentumB_ == 0)
      analyticalDipole[dimension] = dipoleSP(exponentB_, Rb_[dimension], Ra_[dimension], dimension);
    else if (angularMomentumA_ == 1 && angularMomentumB_ == 1)
      analyticalDipole[dimension] =
          dipolePP(exponentA_, exponentB_, Ra_[dimension], Rb_[dimension], evaluationCoordinate_[dimension]);
    else if (angularMomentumA_ == 0 && angularMomentumB_ == 2)
      analyticalDipole[dimension] =
          dipoleSD(exponentA_, exponentB_, Ra_[dimension], Rb_[dimension], evaluationCoordinate_[dimension]);
    else if (angularMomentumA_ == 2 && angularMomentumB_ == 0)
      analyticalDipole[dimension] =
          dipoleSD(exponentB_, exponentA_, Rb_[dimension], Ra_[dimension], evaluationCoordinate_[dimension]);
    else if (angularMomentumA_ == 1 && angularMomentumB_ == 2)
      analyticalDipole[dimension] =
          dipolePD(exponentA_, exponentB_, Ra_[dimension], Rb_[dimension], evaluationCoordinate_[dimension]);
    else if (angularMomentumA_ == 2 && angularMomentumB_ == 1)
      analyticalDipole[dimension] =
          dipolePD(exponentB_, exponentA_, Rb_[dimension], Ra_[dimension], evaluationCoordinate_[dimension]);
    else
      analyticalDipole[dimension] =
          dipoleDD(exponentA_, exponentB_, Ra_[dimension], Rb_[dimension], evaluationCoordinate_[dimension]);
  }
  return analyticalDipole;
}

double AnalyticalDipoleIntegralOverGTOsCalculator::dipoleSS(int dimension) {
  return weightedSum_(dimension) - expSum_ * evaluationCoordinate_[dimension];
}

double AnalyticalDipoleIntegralOverGTOsCalculator::dipoleSP(double expA, double Ra, double Rb, int dimension) {
  double Rab = Ra - Rb;
  return 0.5 + expA * Rab * (weightedSum_(dimension) / expSum_ - evaluationCoordinate_(dimension));
}

double AnalyticalDipoleIntegralOverGTOsCalculator::dipolePP(double expA, double expB, double rA, double rB, double rC) {
  auto const rAB = rA - rB;
  auto const firstTerm = expB * expB * (rA - 2 * rB + rC);
  auto const secondTerm = expA * expA * (2 * rA * (expB * (rAB * rAB) - 1) + rB + rC - 2 * expB * rAB * rAB * rC);
  auto const thirdTerm = expA * expB *
                         (rB * (2 * expB * rB * (rB - rC) - 1) + 2 * expB * rA * rA * (rB - rC) + 2 * rC +
                          rA * (4 * expB * rB * (rC - rB) - 1));
  return -1. * (firstTerm + secondTerm + thirdTerm) / (2. * expSum_ * expSum_);
}

double AnalyticalDipoleIntegralOverGTOsCalculator::dipoleSD(double expA, double expB, double rA, double rB, double rC) {
  auto const rAB = rA - rB;
  auto const firstTerm = expA * expB * (3 * rA - rB - 2 * rC);
  auto const secondTerm = expA * 2 * expA *
                          (expA * rAB * rAB * (rA - rC) + rB * (expB * rB * (rB - rC) - 1) +
                           expB * rA * rA * (rB - rC) - 0.5 * rC + rA * (1.5 + 2 * expB * rB * (rC - rB)));
  auto const thirdTerm = expB * expB * (rB - rC);
  return (firstTerm + secondTerm + thirdTerm) / (2 * expSum_ * expSum_);
}

double AnalyticalDipoleIntegralOverGTOsCalculator::dipolePD(double expA, double expB, double rA, double rB, double rC) {
  auto const rAB = rA - rB;
  auto const rAB2 = rAB * rAB;
  auto const expB2rAB2 = expB * 2. * rAB2;
  auto const firstTerm = -3. * expA * expB * (expB2rAB2 - 2.);
  auto const secondTerm = -2. * expA * expA * expA * rAB * (rA * (expB2rAB2 - 3.) + rB + rC * (2 - expB2rAB2));
  auto const thirdTerm = expA * expA * (3. + 2. * expB * (expB2rAB2 - 3.) * rAB * (rC - rB));
  auto const fourthTerm = expB * expB * (3. - 2. * expB * (rB - rC) * rAB);

  return (firstTerm + secondTerm + thirdTerm + fourthTerm) / (4 * expSum_ * expSum_ * expSum_);
}

double AnalyticalDipoleIntegralOverGTOsCalculator::dipoleDD(double expA, double expB, double rA, double rB, double rC) {
  auto const expSum2 = expSum_ * expSum_;
  auto const rA2 = rA * rA;
  auto const rB2 = rB * rB;
  auto const rB3 = rB2 * rB;
  auto const rA3 = rA2 * rA;
  auto const expA2 = expA * expA;
  auto const expB2 = expB * expB;

  auto const term1 = expB2 * expB *
                     (9 * rB + 2 * expB * rB3 - rA * (6.0 + 4 * expB * rB * (rB - rC)) + 2 * expB * rA2 * (rB - rC) -
                      3 * rC - 2 * expB * rB2 * rC);
  auto const term2 = expA * expB2 *
                     (6 * expB * rA2 * rA + 12 * rB - 10 * expB * rB3 + rA * (2 * expB * rB * (13 * rB - 4 * rC) - 3) -
                      9 * rC + 4 * expB * rB2 * rC + expB * rA2 * (4 * rC - 22 * rB));
  auto const term3 = expA2 * expB *
                     (-3 * rB - 6 * expB * rB3 + 4 * expB2 * rB3 * rB2 - 2 * expB * rA3 * (3 + 8 * expB * rB * (rB - rC)) -
                      2 * rA * (-6 - 3 * expB * rB * (rB - 4 * rC) + 8 * expB2 * rB3 * (rB - rC)) +
                      4 * expB2 * rA2 * rA2 * (rB - rC) - 9 * rC + 12 * expB * rB2 * rC - 4 * expB2 * rB2 * rB2 * rC +
                      6 * expB * rA2 * (rB + 4 * expB * rB3 + 2 * rC - 4 * expB * rB2 * rC));

  auto const term4 =
      expA * expA2 *
      (4 * expB2 * rA2 * rA3 - 6 * rB + 6 * expB * rB3 - 3 * rC + 4 * expB * rB2 * rC - 4 * expB2 * rB2 * rB2 * rC -
       4 * expB2 * rA2 * rA2 * (4 * rB + rC) + 2 * expB * rA2 * (13 * rB - 8 * expB * rB3 + 2 * rC - 12 * expB * rB2 * rC) +
       2 * expB * rA3 * (4 * expB2 * rB * (3 * rB + 2 * rC) - 5) +
       rA * (9 + 4 * expB2 * rB3 * (rB + 4 * rC) - 2 * expB * rB * (11 * rB + 4 * rC)));

  return (term1 + term2 + term3 + term4) / (4 * expSum2 * expSum2);
}

} // namespace Sparrow
} // namespace Scine
