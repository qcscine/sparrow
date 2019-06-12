/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_GLOBAL2C2EMATRIX_H
#define SPARROW_GLOBAL2C2EMATRIX_H

#include "Global2c2eTerms.h"
#include "Local2c2eMatrix.h"
#include "multipoleTypes.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsTypesHelper.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {
class ChargeSeparationParameter;
class KlopmanParameter;

/*!
 * This class calculates the two-center two-electron integrals
 * in the global coordinate system.
 */

class Global2c2eMatrix {
 public:
  using orb_index_t = int;
  using orbPair_index_t = int;

  /**
   *  @brief constructor of the Global2c2eMatrix, constaining the 2 centers 2 electrons integrals between two atoms.
   * @param l1 orbital angular momentum quantum number of the first orbital
   * @param l2 orbital angular momentum quantum number of the second orbital
   * @param D1
   * @param D2
   * @param r1
   * @param r2
   */
  Global2c2eMatrix(int l1, int l2, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& D2,
                   const KlopmanParameter& r1, const KlopmanParameter& r2);

  /**
   * @brief set whether the class refers to integrals within two elements of the same type.
   * @param sym true if the elements are the same
   */
  void setSymmetric(bool sym) {
    localZero_.setSymmetric(sym);
    localOne_.setSymmetric(sym);
    sameElement_ = sym;
  }

  template<Utils::derivOrder O>
  void calculate(const Eigen::Vector3d& Rab);

  template<Utils::derivativeType O>
  Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative(orb_index_t o1, orb_index_t o2, orb_index_t o3,
                                                                   orb_index_t o4) const;
  template<Utils::derivativeType O>
  Utils::AutomaticDifferentiation::DerivativeType<O> getDerivative(orbPair_index_t op1, orbPair_index_t op2) const;
  double get(orb_index_t o1, orb_index_t o2, orb_index_t o3, orb_index_t o4) const;
  double get(orbPair_index_t op1, orbPair_index_t op2) const;
  orbPair_index_t getPairIndex(orb_index_t o1, orb_index_t o2) const;
  void output() const;

 private:
  Eigen::RowVector3d getFirstDerivative(orbPair_index_t op1, orbPair_index_t op2) const;
  Utils::AutomaticDifferentiation::Second3D getSecondDerivative(orbPair_index_t op1, orbPair_index_t op2) const;

  template<Utils::derivOrder O>
  void evaluate();

  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value3DType<O> evaluateMatrixElement(orbPair_index_t op1, orbPair_index_t op2);

  const int d1_, d2_;
  Global2c2eTerms terms_;
  const ChargeSeparationParameter &dist1, &dist2;
  const KlopmanParameter &rho1, &rho2;
  Local2c2eMatrix<Utils::derivOrder::zero> localZero_;
  Local2c2eMatrix<Utils::derivOrder::one> localOne_;
  Local2c2eMatrix<Utils::derivOrder::two> localTwo_;
  OrbitalRotation<Utils::derivOrder::zero> rotationsZero_;
  OrbitalRotation<Utils::derivOrder::one> rotationsOne_;
  OrbitalRotation<Utils::derivOrder::two> rotationsTwo_;
  double R_, RNormX_, RNormY_, RNormZ_;
  Utils::AutomaticDifferentiation::First3D nullDeriv, oneDeriv;
  Eigen::MatrixXd globalMatrix_;
  Eigen::Matrix<Utils::AutomaticDifferentiation::First3D, Eigen::Dynamic, Eigen::Dynamic> globalMatrixOne_;
  Eigen::Matrix<Utils::AutomaticDifferentiation::Second3D, Eigen::Dynamic, Eigen::Dynamic> globalMatrixTwo_;
  bool sameElement_;
  TwoElectronIntegralIndexes pairIndexes_;
};

inline Eigen::RowVector3d Global2c2eMatrix::getFirstDerivative(orbPair_index_t op1, orbPair_index_t op2) const {
  if (op1 == 100 || op2 == 100)
    return Utils::Gradient::Zero();

  return globalMatrixOne_(op1, op2).derivatives();
}

inline Utils::AutomaticDifferentiation::Second3D Global2c2eMatrix::getSecondDerivative(orbPair_index_t op1,
                                                                                       orbPair_index_t op2) const {
  if (op1 == 100 || op2 == 100)
    return {};

  return globalMatrixTwo_(op1, op2);
}

template<Utils::derivativeType O>
Utils::AutomaticDifferentiation::DerivativeType<O> Global2c2eMatrix::getDerivative(orb_index_t o1, orb_index_t o2,
                                                                                   orb_index_t o3, orb_index_t o4) const {
  return getDerivative<O>(getPairIndex(o1, o2), getPairIndex(o3, o4));
}

template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::first>
Global2c2eMatrix::getDerivative<Utils::derivativeType::first>(orbPair_index_t op1, orbPair_index_t op2) const {
  return getFirstDerivative(op1, op2);
}

template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_atomic>
Global2c2eMatrix::getDerivative<Utils::derivativeType::second_atomic>(orbPair_index_t op1, orbPair_index_t op2) const {
  return getSecondDerivative(op1, op2);
}

template<>
inline Utils::AutomaticDifferentiation::DerivativeType<Utils::derivativeType::second_full>
Global2c2eMatrix::getDerivative<Utils::derivativeType::second_full>(orbPair_index_t op1, orbPair_index_t op2) const {
  return getSecondDerivative(op1, op2);
}

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_GLOBAL2C2EMATRIX_H
