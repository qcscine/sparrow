/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Global2c2eMatrix.h"
#include <Utils/Constants.h>
#include <iostream>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {
using namespace GeneralTypes;

namespace multipole {

Global2c2eMatrix::Global2c2eMatrix(int l1, int l2, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& D2,
                                   const KlopmanParameter& r1, const KlopmanParameter& r2)
  : d1_(l1 == 0 ? 1 : l1 == 1 ? 10 : 40),
    d2_(l2 == 0 ? 1 : l2 == 1 ? 10 : 40),
    dist1(D1),
    dist2(D2),
    rho1(r1),
    rho2(r2),
    localZero_(l1, l2, D1, D2, r1, r2),
    localOne_(l1, l2, D1, D2, r1, r2),
    localTwo_(l1, l2, D1, D2, r1, r2),
    rotationsZero_(l1, l2),
    rotationsOne_(l1, l2),
    rotationsTwo_(l1, l2) {
  sameElement_ = false;
  nullDeriv = First3D(0, 0, 0, 0);
  oneDeriv = First3D(1, 0, 0, 0);
  globalMatrix_ = Eigen::MatrixXd::Zero(d1_, d2_);
  globalMatrixOne_ = Eigen::Matrix<First3D, Eigen::Dynamic, Eigen::Dynamic>::Constant(
      d1_, d2_, constant3D<Utils::DerivativeOrder::One>(0.0));
  globalMatrixTwo_ = Eigen::Matrix<Second3D, Eigen::Dynamic, Eigen::Dynamic>::Constant(
      d1_, d2_, constant3D<Utils::DerivativeOrder::Two>(0.0));
}

template<>
Value3DType<Utils::DerivativeOrder::Zero>
Global2c2eMatrix::evaluateMatrixElement<Utils::DerivativeOrder::Zero>(orbPair_index_t op1, orbPair_index_t op2) {
  double element = 0;
  const std::list<RotationTerm>& expr = terms_.getTermList(op1, op2);
  for (const auto& term : expr) {
    element += rotationsZero_.getRotationCoefficient(term.f1_) * rotationsZero_.getRotationCoefficient(term.f2_) *
               rotationsZero_.getRotationCoefficient(term.f3_) * rotationsZero_.getRotationCoefficient(term.f4_) *
               localZero_(term.pair1_, term.pair2_);
  }
  return element;
}
template<>
Value3DType<Utils::DerivativeOrder::One>
Global2c2eMatrix::evaluateMatrixElement<Utils::DerivativeOrder::One>(orbPair_index_t op1, orbPair_index_t op2) {
  First3D element = nullDeriv;
  const std::list<RotationTerm>& expr = terms_.getTermList(op1, op2);
  for (const auto& term : expr) {
    element += rotationsOne_.getRotationCoefficient(term.f1_) * rotationsOne_.getRotationCoefficient(term.f2_) *
               rotationsOne_.getRotationCoefficient(term.f3_) * rotationsOne_.getRotationCoefficient(term.f4_) *
               First3D(localOne_(term.pair1_, term.pair2_).value(), localOne_(term.pair1_, term.pair2_).derivative() * RNormX_,
                       localOne_(term.pair1_, term.pair2_).derivative() * RNormY_,
                       localOne_(term.pair1_, term.pair2_).derivative() * RNormZ_);
  }
  return element;
}
template<>
Value3DType<Utils::DerivativeOrder::Two>
Global2c2eMatrix::evaluateMatrixElement<Utils::DerivativeOrder::Two>(orbPair_index_t op1, orbPair_index_t op2) {
  Second3D element;
  element.setZero();
  const std::list<RotationTerm>& expr = terms_.getTermList(op1, op2);
  for (const auto& term : expr) {
    const Second1D& localElement = localTwo_(term.pair1_, term.pair2_);
    double localElFirstByR = localElement.first() / R_;
    Second3D transformedElement(localElement.value(), localElement.first() * RNormX_, localElement.first() * RNormY_,
                                localElement.first() * RNormZ_,
                                localElFirstByR * (1 - RNormX_ * RNormX_) + RNormX_ * RNormX_ * localElement.second(),
                                localElFirstByR * (1 - RNormY_ * RNormY_) + RNormY_ * RNormY_ * localElement.second(),
                                localElFirstByR * (1 - RNormZ_ * RNormZ_) + RNormZ_ * RNormZ_ * localElement.second(),
                                RNormX_ * RNormY_ * (localElement.second() - localElFirstByR),
                                RNormX_ * RNormZ_ * (localElement.second() - localElFirstByR),
                                RNormY_ * RNormZ_ * (localElement.second() - localElFirstByR));
    element += rotationsTwo_.getRotationCoefficient(term.f1_) * rotationsTwo_.getRotationCoefficient(term.f2_) *
               rotationsTwo_.getRotationCoefficient(term.f3_) * rotationsTwo_.getRotationCoefficient(term.f4_) *
               std::move(transformedElement);
  }
  return element;
}

template<>
void Global2c2eMatrix::evaluate<Utils::DerivativeOrder::Zero>() {
  for (int i = 0; i < d1_; i++)
    for (int j = 0; j < d2_; j++)
      globalMatrix_(i, j) = evaluateMatrixElement<Utils::DerivativeOrder::Zero>(i, j);
}
template<>
void Global2c2eMatrix::evaluate<Utils::DerivativeOrder::One>() {
  for (int i = 0; i < d1_; i++)
    for (int j = 0; j < d2_; j++) {
      globalMatrixOne_(i, j) = evaluateMatrixElement<Utils::DerivativeOrder::One>(i, j);
      globalMatrix_(i, j) = globalMatrixOne_(i, j).value();
    }
}
template<>
void Global2c2eMatrix::evaluate<Utils::DerivativeOrder::Two>() {
  for (int i = 0; i < d1_; i++)
    for (int j = 0; j < d2_; j++) {
      globalMatrixTwo_(i, j) = evaluateMatrixElement<Utils::DerivativeOrder::Two>(i, j);
      globalMatrix_(i, j) = globalMatrixTwo_(i, j).value();
    }
}

template<>
void Global2c2eMatrix::calculate<Utils::DerivativeOrder::Zero>(const Eigen::Vector3d& Rab) {
  R_ = rotationsZero_.setVector(Rab);
  RNormX_ = Rab[0] / R_;
  RNormY_ = Rab[1] / R_;
  RNormZ_ = Rab[2] / R_;

  localZero_.calculate(R_);
  rotationsZero_.evaluate();
  evaluate<Utils::DerivativeOrder::Zero>();
}
template<>
void Global2c2eMatrix::calculate<Utils::DerivativeOrder::One>(const Eigen::Vector3d& Rab) {
  R_ = rotationsOne_.setVector(Rab);
  RNormX_ = Rab[0] / R_;
  RNormY_ = Rab[1] / R_;
  RNormZ_ = Rab[2] / R_;

  localOne_.calculate(R_);
  rotationsOne_.evaluate();
  evaluate<Utils::DerivativeOrder::One>();
}
template<>
void Global2c2eMatrix::calculate<Utils::DerivativeOrder::Two>(const Eigen::Vector3d& Rab) {
  R_ = rotationsTwo_.setVector(Rab);
  RNormX_ = Rab[0] / R_;
  RNormY_ = Rab[1] / R_;
  RNormZ_ = Rab[2] / R_;

  localTwo_.calculate(R_);
  rotationsTwo_.evaluate();
  evaluate<Utils::DerivativeOrder::Two>();
}

double Global2c2eMatrix::get(orb_index_t o1, orb_index_t o2, orb_index_t o3, orb_index_t o4) const {
  return get(getPairIndex(o1, o2), getPairIndex(o3, o4));
}

double Global2c2eMatrix::get(orbPair_index_t op1, orbPair_index_t op2) const {
  if (op1 == 100 || op2 == 100)
    return 0;

  return globalMatrix_(op1, op2);
}

Global2c2eMatrix::orbPair_index_t Global2c2eMatrix::getPairIndex(orb_index_t o1, orb_index_t o2) const {
  return TwoElectronIntegralIndexes::getPairIndex(o1, o2);
}

void Global2c2eMatrix::output() const {
  std::cout << "output: 2c2e\n";
  int max1 = (d1_ == 1) ? 1 : (d1_ == 10) ? 4 : 9;
  int max2 = (d2_ == 1) ? 1 : (d2_ == 10) ? 4 : 9;
  int d = 0;
  std::cout << "n of orbs: " << max1 << " " << max2 << std::endl;
  for (int i = 0; i < max1; i++)
    for (int j = 0; j <= i; j++)
      for (int k = 0; k < max2; k++)
        for (int l = 0; l <= k; l++) {
          d++;
          std::cout << get(i, j, k, l) * Utils::Constants::ev_per_hartree << " ";
          if (d == 10) {
            d = 0;
            std::cout << std::endl;
          }
        }
  std::cout << std::endl;
  std::cout << "end of output." << std::endl;
}

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
