/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Local2c2eMatrix.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>

namespace Scine {
namespace Sparrow {

namespace nddo {
using namespace GeneralTypes;

namespace multipole {

template<Utils::derivOrder O>
Local2c2eMatrix<O>::Local2c2eMatrix(int l1, int l2, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& D2,
                                    const KlopmanParameter& r1, const KlopmanParameter& r2)
  : l1_(l1), l2_(l2), dist1(D1), dist2(D2), rho1(r1), rho2(r2) {
  sameElement_ = false;
  const int d1 = l1_ == 0 ? 1 : l1_ == 1 ? 10 : 40;
  const int d2 = l2_ == 0 ? 1 : l2_ == 1 ? 10 : 40;
  mat = Eigen::Matrix<Utils::AutomaticDifferentiation::Value1DType<O>, Eigen::Dynamic, Eigen::Dynamic>::Constant(
      d1, d2, Utils::AutomaticDifferentiation::constant1D<O>(0.0));
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::calculate(double R) {
  if (sameElement_)
    calculateSym(R);
  else
    calculateAsym(R);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::calculateSym(double R) {
  // ss
  buildSSMatrix(R);

  // sp, ps, pp
  if (l1_ > 0) {
    buildSPMatrix(R);
    buildPPMatrixSym(R);
    buildPSMatrixSym(R);
  }

  // sd,ds, pd,dp
  if (l2_ > 1) {
    buildSDMatrix(R);
    buildPDMatrix(R);
    buildDDMatrixSym(R);
    buildDSMatrixSym(R);
    buildDPMatrixSym(R);
  }
}
template<Utils::derivOrder O>
void Local2c2eMatrix<O>::calculateAsym(double R) {
  // ss
  buildSSMatrix(R);

  // sp, ps
  if (l2_ > 0)
    buildSPMatrix(R);
  if (l1_ > 0)
    buildPSMatrix(R);

  // sd,ds, pd,dp
  if (l2_ > 1) {
    buildSDMatrix(R);
    if (l1_ > 0)
      buildPDMatrix(R);
  }
  if (l1_ > 1) {
    buildDSMatrix(R);
    if (l2_ > 0)
      buildDPMatrix(R);
  }

  // pp, dd
  if (l1_ > 0 && l2_ > 0) {
    buildPPMatrix(R);
    if (l1_ > 1 && l2_ > 1)
      buildDDMatrix(R);
  }
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildSSMatrix(double R) {
  mat(0, 0) = calculator_.getIntegral<O>(0, 0, R, dist1, dist2, rho1, rho2);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildSPMatrix(double R) {
  mat(0, 2) = calculator_.getIntegral<O>(0, 2, R, dist1, dist2, rho1, rho2);
  mat(0, 5) = mat(0, 2);
  mat(0, 6) = calculator_.getIntegral<O>(0, 6, R, dist1, dist2, rho1, rho2);
  mat(0, 9) = calculator_.getIntegral<O>(0, 9, R, dist1, dist2, rho1, rho2);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildPSMatrix(double R) {
  mat(2, 0) = calculator_.getIntegral<O>(2, 0, R, dist1, dist2, rho1, rho2);
  mat(5, 0) = mat(2, 0);
  mat(6, 0) = calculator_.getIntegral<O>(6, 0, R, dist1, dist2, rho1, rho2);
  mat(9, 0) = calculator_.getIntegral<O>(9, 0, R, dist1, dist2, rho1, rho2);
}
template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildPSMatrixSym(double /*R*/) {
  mat(2, 0) = mat(0, 2);
  mat(5, 0) = mat(0, 5);
  mat(6, 0) = -mat(0, 6);
  mat(9, 0) = mat(0, 9);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildPPMatrix(double R) {
  mat(1, 1) = calculator_.getIntegral<O>(1, 1, R, dist1, dist2, rho1, rho2);
  mat(1, 7) = calculator_.getIntegral<O>(1, 7, R, dist1, dist2, rho1, rho2);
  mat(2, 2) = calculator_.getIntegral<O>(2, 2, R, dist1, dist2, rho1, rho2);
  mat(2, 5) = calculator_.getIntegral<O>(2, 5, R, dist1, dist2, rho1, rho2);
  mat(2, 6) = calculator_.getIntegral<O>(2, 6, R, dist1, dist2, rho1, rho2);
  mat(2, 9) = calculator_.getIntegral<O>(2, 9, R, dist1, dist2, rho1, rho2);
  mat(3, 3) = mat(1, 1);
  mat(3, 8) = mat(1, 7);
  // mat(4,4) = methodWrapper_.getIntegral<O>(4,4,R,dist1,dist2,rho1,rho2);
  mat(5, 2) = mat(2, 5);
  mat(5, 5) = mat(2, 2);
  mat(5, 6) = mat(2, 6);
  mat(5, 9) = mat(2, 9);
  mat(6, 2) = calculator_.getIntegral<O>(6, 2, R, dist1, dist2, rho1, rho2);
  mat(6, 5) = mat(6, 2);
  mat(6, 6) = calculator_.getIntegral<O>(6, 6, R, dist1, dist2, rho1, rho2);
  mat(6, 9) = calculator_.getIntegral<O>(6, 9, R, dist1, dist2, rho1, rho2);
  mat(7, 1) = calculator_.getIntegral<O>(7, 1, R, dist1, dist2, rho1, rho2);
  mat(7, 7) = calculator_.getIntegral<O>(7, 7, R, dist1, dist2, rho1, rho2);
  mat(8, 3) = mat(7, 1);
  mat(8, 8) = mat(7, 7);
  mat(9, 2) = calculator_.getIntegral<O>(9, 2, R, dist1, dist2, rho1, rho2);
  mat(9, 5) = mat(9, 2);
  mat(9, 6) = calculator_.getIntegral<O>(9, 6, R, dist1, dist2, rho1, rho2);
  mat(9, 9) = calculator_.getIntegral<O>(9, 9, R, dist1, dist2, rho1, rho2);

  mat(static_cast<int>(twoElIntegral_t::x_y), static_cast<int>(twoElIntegral_t::x_y)) =
      0.5 * (mat(static_cast<int>(twoElIntegral_t::x_x), static_cast<int>(twoElIntegral_t::x_x)) -
             mat(static_cast<int>(twoElIntegral_t::x_x), static_cast<int>(twoElIntegral_t::y_y)));
}
template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildPPMatrixSym(double R) {
  mat(1, 1) = calculator_.getIntegral<O>(1, 1, R, dist1, dist2, rho1, rho2);
  mat(1, 7) = calculator_.getIntegral<O>(1, 7, R, dist1, dist2, rho1, rho2);
  mat(2, 2) = calculator_.getIntegral<O>(2, 2, R, dist1, dist2, rho1, rho2);
  mat(2, 5) = calculator_.getIntegral<O>(2, 5, R, dist1, dist2, rho1, rho2);
  mat(2, 6) = calculator_.getIntegral<O>(2, 6, R, dist1, dist2, rho1, rho2);
  mat(2, 9) = calculator_.getIntegral<O>(2, 9, R, dist1, dist2, rho1, rho2);
  mat(3, 3) = mat(1, 1);
  mat(3, 8) = mat(1, 7);
  // mat(4,4) = methodWrapper_.getIntegral<O>(4,4,R,dist1,dist2,rho1,rho2);
  mat(5, 2) = mat(2, 5);
  mat(5, 5) = mat(2, 2);
  mat(5, 6) = mat(2, 6);
  mat(5, 9) = mat(2, 9);
  mat(6, 2) = -mat(2, 6);
  mat(6, 5) = -mat(2, 6);
  mat(6, 6) = calculator_.getIntegral<O>(6, 6, R, dist1, dist2, rho1, rho2);
  mat(6, 9) = calculator_.getIntegral<O>(6, 9, R, dist1, dist2, rho1, rho2);
  mat(7, 1) = -mat(1, 7);
  mat(7, 7) = calculator_.getIntegral<O>(7, 7, R, dist1, dist2, rho1, rho2);
  mat(8, 3) = -mat(1, 7);
  mat(8, 8) = mat(7, 7);
  mat(9, 2) = mat(2, 9);
  mat(9, 5) = mat(2, 9);
  mat(9, 6) = -mat(6, 9);
  mat(9, 9) = calculator_.getIntegral<O>(9, 9, R, dist1, dist2, rho1, rho2);
  mat(static_cast<int>(twoElIntegral_t::x_y), static_cast<int>(twoElIntegral_t::x_y)) =
      0.5 * (mat(static_cast<int>(twoElIntegral_t::x_x), static_cast<int>(twoElIntegral_t::x_x)) -
             mat(static_cast<int>(twoElIntegral_t::x_x), static_cast<int>(twoElIntegral_t::y_y)));
}
template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildDDMatrixSym(double R) {
  mat(10, 10) = calculator_.getIntegral<O>(10, 10, R, dist1, dist2, rho1, rho2);
  mat(10, 16) = calculator_.getIntegral<O>(10, 16, R, dist1, dist2, rho1, rho2);
  mat(10, 20) = mat(10, 16);
  mat(10, 23) = calculator_.getIntegral<O>(10, 23, R, dist1, dist2, rho1, rho2);
  mat(10, 26) = calculator_.getIntegral<O>(10, 26, R, dist1, dist2, rho1, rho2);
  mat(10, 31) = calculator_.getIntegral<O>(10, 31, R, dist1, dist2, rho1, rho2);
  mat(10, 35) = mat(10, 31);
  mat(10, 38) = calculator_.getIntegral<O>(10, 38, R, dist1, dist2, rho1, rho2);
  mat(10, 39) = mat(10, 38);
  mat(11, 11) = calculator_.getIntegral<O>(11, 11, R, dist1, dist2, rho1, rho2);
  mat(11, 15) = calculator_.getIntegral<O>(11, 15, R, dist1, dist2, rho1, rho2);
  mat(11, 17) = calculator_.getIntegral<O>(11, 17, R, dist1, dist2, rho1, rho2);
  mat(11, 22) = mat(11, 17);
  mat(11, 24) = mat(11, 17);
  mat(11, 27) = calculator_.getIntegral<O>(11, 27, R, dist1, dist2, rho1, rho2);
  mat(11, 33) = calculator_.getIntegral<O>(11, 33, R, dist1, dist2, rho1, rho2);
  mat(11, 37) = mat(11, 33);
  mat(12, 12) = mat(11, 11);
  mat(12, 18) = mat(11, 17);
  mat(12, 19) = mat(11, 15);
  mat(12, 21) = -mat(11, 17);
  mat(12, 25) = mat(11, 17);
  mat(12, 28) = mat(11, 27);
  mat(12, 34) = mat(11, 33);
  mat(12, 36) = -mat(11, 33);
  mat(13, 13) = calculator_.getIntegral<O>(13, 13, R, dist1, dist2, rho1, rho2);
  mat(13, 29) = calculator_.getIntegral<O>(13, 29, R, dist1, dist2, rho1, rho2);
  mat(13, 31) = calculator_.getIntegral<O>(13, 31, R, dist1, dist2, rho1, rho2);
  mat(13, 35) = -mat(13, 31);
  mat(14, 14) = mat(13, 13);
  mat(14, 30) = mat(13, 29);
  mat(14, 32) = mat(13, 31);
  mat(15, 11) = -mat(11, 15);
  mat(15, 15) = calculator_.getIntegral<O>(15, 15, R, dist1, dist2, rho1, rho2);
  mat(15, 17) = calculator_.getIntegral<O>(15, 17, R, dist1, dist2, rho1, rho2);
  mat(15, 22) = mat(15, 17);
  mat(15, 24) = mat(15, 17);
  mat(15, 27) = calculator_.getIntegral<O>(15, 27, R, dist1, dist2, rho1, rho2);
  mat(15, 33) = calculator_.getIntegral<O>(15, 33, R, dist1, dist2, rho1, rho2);
  mat(15, 37) = mat(15, 33);
  mat(16, 10) = -mat(10, 16);
  mat(16, 16) = calculator_.getIntegral<O>(16, 16, R, dist1, dist2, rho1, rho2);
  mat(16, 20) = mat(16, 16);
  mat(16, 23) = calculator_.getIntegral<O>(16, 23, R, dist1, dist2, rho1, rho2);
  mat(16, 26) = calculator_.getIntegral<O>(16, 26, R, dist1, dist2, rho1, rho2);
  mat(16, 31) = calculator_.getIntegral<O>(16, 31, R, dist1, dist2, rho1, rho2);
  mat(16, 35) = mat(16, 31);
  mat(16, 38) = calculator_.getIntegral<O>(16, 38, R, dist1, dist2, rho1, rho2);
  mat(16, 39) = mat(16, 38);
  mat(17, 11) = -mat(11, 17);
  mat(17, 15) = mat(15, 17);
  mat(17, 17) = calculator_.getIntegral<O>(17, 17, R, dist1, dist2, rho1, rho2);
  mat(17, 22) = mat(17, 17);
  mat(17, 24) = mat(17, 17);
  mat(17, 27) = -mat(15, 33);
  mat(17, 33) = calculator_.getIntegral<O>(17, 33, R, dist1, dist2, rho1, rho2);
  mat(17, 37) = mat(17, 33);
  mat(18, 12) = -mat(11, 17);
  mat(18, 18) = mat(17, 17);
  mat(18, 19) = mat(15, 17);
  mat(18, 21) = -mat(17, 17);
  mat(18, 25) = mat(17, 17);
  mat(18, 28) = -mat(15, 33);
  mat(18, 34) = mat(17, 33);
  mat(18, 36) = -mat(17, 33);
  mat(19, 12) = -mat(11, 15);
  mat(19, 18) = mat(15, 17);
  mat(19, 19) = mat(15, 15);
  mat(19, 21) = -mat(15, 17);
  mat(19, 25) = mat(15, 17);
  mat(19, 28) = mat(15, 27);
  mat(19, 34) = mat(15, 33);
  mat(19, 36) = -mat(15, 33);
  mat(20, 10) = -mat(10, 16);
  mat(20, 16) = mat(16, 16);
  mat(20, 20) = mat(16, 16);
  mat(20, 23) = mat(16, 23);
  mat(20, 26) = mat(16, 26);
  mat(20, 31) = mat(16, 31);
  mat(20, 35) = mat(16, 31);
  mat(20, 38) = mat(16, 38);
  mat(20, 39) = mat(16, 38);
  mat(21, 12) = mat(11, 17);
  mat(21, 18) = -mat(17, 17);
  mat(21, 19) = -mat(15, 17);
  mat(21, 21) = mat(17, 17);
  mat(21, 25) = -mat(17, 17);
  mat(21, 28) = mat(15, 33);
  mat(21, 34) = -mat(17, 33);
  mat(21, 36) = mat(17, 33);
  mat(22, 11) = -mat(11, 17);
  mat(22, 15) = mat(15, 17);
  mat(22, 17) = mat(17, 17);
  mat(22, 22) = mat(17, 17);
  mat(22, 24) = mat(17, 17);
  mat(22, 27) = -mat(15, 33);
  mat(22, 33) = mat(17, 33);
  mat(22, 37) = mat(17, 33);
  mat(23, 10) = -mat(10, 23);
  mat(23, 16) = mat(16, 23);
  mat(23, 20) = mat(16, 23);
  mat(23, 23) = calculator_.getIntegral<O>(23, 23, R, dist1, dist2, rho1, rho2);
  mat(23, 26) = calculator_.getIntegral<O>(23, 26, R, dist1, dist2, rho1, rho2);
  mat(23, 31) = calculator_.getIntegral<O>(23, 31, R, dist1, dist2, rho1, rho2);
  mat(23, 35) = mat(23, 31);
  mat(23, 38) = calculator_.getIntegral<O>(23, 38, R, dist1, dist2, rho1, rho2);
  mat(23, 39) = mat(23, 38);
  mat(24, 11) = -mat(11, 17);
  mat(24, 15) = mat(15, 17);
  mat(24, 17) = mat(17, 17);
  mat(24, 22) = mat(17, 17);
  mat(24, 24) = mat(17, 17);
  mat(24, 27) = -mat(15, 33);
  mat(24, 33) = mat(17, 33);
  mat(24, 37) = mat(17, 33);
  mat(25, 12) = -mat(11, 17);
  mat(25, 18) = mat(17, 17);
  mat(25, 19) = mat(15, 17);
  mat(25, 21) = -mat(17, 17);
  mat(25, 25) = mat(17, 17);
  mat(25, 28) = -mat(15, 33);
  mat(25, 34) = mat(17, 33);
  mat(25, 36) = -mat(17, 33);
  mat(26, 10) = mat(10, 26);
  mat(26, 16) = -mat(16, 26);
  mat(26, 20) = -mat(16, 26);
  mat(26, 23) = -mat(23, 26);
  mat(26, 26) = calculator_.getIntegral<O>(26, 26, R, dist1, dist2, rho1, rho2);
  mat(26, 31) = calculator_.getIntegral<O>(26, 31, R, dist1, dist2, rho1, rho2);
  mat(26, 35) = mat(26, 31);
  mat(26, 38) = calculator_.getIntegral<O>(26, 38, R, dist1, dist2, rho1, rho2);
  mat(26, 39) = mat(26, 38);
  mat(27, 11) = mat(11, 27);
  mat(27, 15) = -mat(15, 27);
  mat(27, 17) = mat(15, 33);
  mat(27, 22) = mat(15, 33);
  mat(27, 24) = mat(15, 33);
  mat(27, 27) = calculator_.getIntegral<O>(27, 27, R, dist1, dist2, rho1, rho2);
  mat(27, 33) = calculator_.getIntegral<O>(27, 33, R, dist1, dist2, rho1, rho2);
  mat(27, 37) = mat(27, 33);
  mat(28, 12) = mat(11, 27);
  mat(28, 18) = mat(15, 33);
  mat(28, 19) = -mat(15, 27);
  mat(28, 21) = -mat(15, 33);
  mat(28, 25) = mat(15, 33);
  mat(28, 28) = mat(27, 27);
  mat(28, 34) = mat(27, 33);
  mat(28, 36) = -mat(27, 33);
  mat(29, 13) = mat(13, 29);
  mat(29, 29) = calculator_.getIntegral<O>(29, 29, R, dist1, dist2, rho1, rho2);
  mat(29, 31) = calculator_.getIntegral<O>(29, 31, R, dist1, dist2, rho1, rho2);
  mat(29, 35) = -mat(29, 31);
  mat(30, 14) = mat(13, 29);
  mat(30, 30) = mat(29, 29);
  mat(30, 32) = mat(29, 31);
  mat(31, 10) = mat(10, 31);
  mat(31, 13) = mat(13, 31);
  mat(31, 16) = -mat(16, 31);
  mat(31, 20) = -mat(16, 31);
  mat(31, 23) = -mat(23, 31);
  mat(31, 26) = mat(26, 31);
  mat(31, 29) = mat(29, 31);
  mat(31, 31) = calculator_.getIntegral<O>(31, 31, R, dist1, dist2, rho1, rho2);
  mat(31, 35) = calculator_.getIntegral<O>(31, 35, R, dist1, dist2, rho1, rho2);
  mat(31, 38) = calculator_.getIntegral<O>(31, 38, R, dist1, dist2, rho1, rho2);
  mat(31, 39) = mat(31, 38);
  mat(32, 14) = mat(13, 31);
  mat(32, 30) = mat(29, 31);
  mat(32, 32) = calculator_.getIntegral<O>(32, 32, R, dist1, dist2, rho1, rho2);
  mat(33, 11) = mat(11, 33);
  mat(33, 15) = -mat(15, 33);
  mat(33, 17) = -mat(17, 33);
  mat(33, 22) = -mat(17, 33);
  mat(33, 24) = -mat(17, 33);
  mat(33, 27) = mat(27, 33);
  mat(33, 33) = calculator_.getIntegral<O>(33, 33, R, dist1, dist2, rho1, rho2);
  mat(33, 37) = mat(33, 33);
  mat(34, 12) = mat(11, 33);
  mat(34, 18) = -mat(17, 33);
  mat(34, 19) = -mat(15, 33);
  mat(34, 21) = mat(17, 33);
  mat(34, 25) = -mat(17, 33);
  mat(34, 28) = mat(27, 33);
  mat(34, 34) = mat(33, 33);
  mat(34, 36) = -mat(33, 33);
  mat(35, 10) = mat(10, 31);
  mat(35, 13) = -mat(13, 31);
  mat(35, 16) = -mat(16, 31);
  mat(35, 20) = -mat(16, 31);
  mat(35, 23) = -mat(23, 31);
  mat(35, 26) = mat(26, 31);
  mat(35, 29) = -mat(29, 31);
  mat(35, 31) = mat(31, 35);
  mat(35, 35) = mat(31, 31);
  mat(35, 38) = mat(31, 38);
  mat(35, 39) = mat(31, 38);
  mat(36, 12) = -mat(11, 33);
  mat(36, 18) = mat(17, 33);
  mat(36, 19) = mat(15, 33);
  mat(36, 21) = -mat(17, 33);
  mat(36, 25) = mat(17, 33);
  mat(36, 28) = -mat(27, 33);
  mat(36, 34) = -mat(33, 33);
  mat(36, 36) = mat(33, 33);
  mat(37, 11) = mat(11, 33);
  mat(37, 15) = -mat(15, 33);
  mat(37, 17) = -mat(17, 33);
  mat(37, 22) = -mat(17, 33);
  mat(37, 24) = -mat(17, 33);
  mat(37, 27) = mat(27, 33);
  mat(37, 33) = mat(33, 33);
  mat(37, 37) = mat(33, 33);
  mat(38, 10) = mat(10, 38);
  mat(38, 16) = -mat(16, 38);
  mat(38, 20) = -mat(16, 38);
  mat(38, 23) = -mat(23, 38);
  mat(38, 26) = mat(26, 38);
  mat(38, 31) = mat(31, 38);
  mat(38, 35) = mat(31, 38);
  mat(38, 38) = calculator_.getIntegral<O>(38, 38, R, dist1, dist2, rho1, rho2);
  mat(38, 39) = mat(38, 38);
  mat(39, 10) = mat(10, 38);
  mat(39, 16) = -mat(16, 38);
  mat(39, 20) = -mat(16, 38);
  mat(39, 23) = -mat(23, 38);
  mat(39, 26) = mat(26, 38);
  mat(39, 31) = mat(31, 38);
  mat(39, 35) = mat(31, 38);
  mat(39, 38) = mat(38, 38);
  mat(39, 39) = mat(38, 38);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildSDMatrix(double R) {
  mat(0, 10) = calculator_.getIntegral<O>(0, 10, R, dist1, dist2, rho1, rho2);
  mat(0, 16) = calculator_.getIntegral<O>(0, 16, R, dist1, dist2, rho1, rho2);
  mat(0, 20) = mat(0, 16);
  mat(0, 23) = calculator_.getIntegral<O>(0, 23, R, dist1, dist2, rho1, rho2);
  mat(0, 26) = calculator_.getIntegral<O>(0, 26, R, dist1, dist2, rho1, rho2);
  mat(0, 31) = calculator_.getIntegral<O>(0, 31, R, dist1, dist2, rho1, rho2);
  mat(0, 35) = mat(0, 31);
  mat(0, 38) = calculator_.getIntegral<O>(0, 38, R, dist1, dist2, rho1, rho2);
  mat(0, 39) = mat(0, 38);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildDSMatrix(double R) {
  mat(10, 0) = calculator_.getIntegral<O>(10, 0, R, dist1, dist2, rho1, rho2);
  mat(16, 0) = calculator_.getIntegral<O>(16, 0, R, dist1, dist2, rho1, rho2);
  mat(20, 0) = mat(16, 0);
  mat(23, 0) = calculator_.getIntegral<O>(23, 0, R, dist1, dist2, rho1, rho2);
  mat(26, 0) = calculator_.getIntegral<O>(26, 0, R, dist1, dist2, rho1, rho2);
  mat(31, 0) = calculator_.getIntegral<O>(31, 0, R, dist1, dist2, rho1, rho2);
  mat(35, 0) = mat(31, 0);
  mat(38, 0) = calculator_.getIntegral<O>(38, 0, R, dist1, dist2, rho1, rho2);
  mat(39, 0) = mat(38, 0);
}
template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildDSMatrixSym(double /*R*/) {
  mat(10, 0) = mat(0, 10);
  mat(16, 0) = -mat(0, 16);
  mat(20, 0) = -mat(0, 20);
  mat(23, 0) = -mat(0, 23);
  mat(26, 0) = mat(0, 26);
  mat(31, 0) = mat(0, 31);
  mat(35, 0) = mat(0, 35);
  mat(38, 0) = mat(0, 38);
  mat(39, 0) = mat(0, 39);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildPDMatrix(double R) {
  mat(1, 11) = calculator_.getIntegral<O>(1, 11, R, dist1, dist2, rho1, rho2);
  mat(1, 15) = calculator_.getIntegral<O>(1, 15, R, dist1, dist2, rho1, rho2);
  mat(1, 17) = calculator_.getIntegral<O>(1, 17, R, dist1, dist2, rho1, rho2);
  mat(1, 22) = mat(1, 17);
  mat(1, 24) = mat(1, 17);
  mat(1, 27) = calculator_.getIntegral<O>(1, 27, R, dist1, dist2, rho1, rho2);
  mat(1, 33) = calculator_.getIntegral<O>(1, 33, R, dist1, dist2, rho1, rho2);
  mat(1, 37) = mat(1, 33);
  mat(2, 10) = calculator_.getIntegral<O>(2, 10, R, dist1, dist2, rho1, rho2);
  mat(2, 13) = calculator_.getIntegral<O>(2, 13, R, dist1, dist2, rho1, rho2);
  mat(2, 16) = calculator_.getIntegral<O>(2, 16, R, dist1, dist2, rho1, rho2);
  mat(2, 20) = mat(2, 16);
  mat(2, 23) = calculator_.getIntegral<O>(2, 23, R, dist1, dist2, rho1, rho2);
  mat(2, 26) = calculator_.getIntegral<O>(2, 26, R, dist1, dist2, rho1, rho2);
  mat(2, 29) = calculator_.getIntegral<O>(2, 29, R, dist1, dist2, rho1, rho2);
  mat(2, 31) = calculator_.getIntegral<O>(2, 31, R, dist1, dist2, rho1, rho2);
  mat(2, 35) = calculator_.getIntegral<O>(2, 35, R, dist1, dist2, rho1, rho2);
  mat(2, 38) = calculator_.getIntegral<O>(2, 38, R, dist1, dist2, rho1, rho2);
  mat(2, 39) = mat(2, 38);
  mat(3, 12) = mat(1, 11);
  mat(3, 18) = mat(1, 17);
  mat(3, 19) = mat(1, 15);
  mat(3, 21) = -mat(1, 17);
  mat(3, 25) = mat(1, 17);
  mat(3, 28) = mat(1, 27);
  mat(3, 34) = mat(1, 33);
  mat(3, 36) = -mat(1, 33);
  mat(4, 14) = mat(2, 13);
  mat(4, 30) = mat(2, 29);
  mat(4, 32) = calculator_.getIntegral<O>(4, 32, R, dist1, dist2, rho1, rho2);
  mat(5, 10) = mat(2, 10);
  mat(5, 13) = -mat(2, 13);
  mat(5, 16) = mat(2, 16);
  mat(5, 20) = mat(2, 16);
  mat(5, 23) = mat(2, 23);
  mat(5, 26) = mat(2, 26);
  mat(5, 29) = -mat(2, 29);
  mat(5, 31) = mat(2, 35);
  mat(5, 35) = mat(2, 31);
  mat(5, 38) = mat(2, 38);
  mat(5, 39) = mat(2, 38);
  mat(6, 10) = calculator_.getIntegral<O>(6, 10, R, dist1, dist2, rho1, rho2);
  mat(6, 16) = calculator_.getIntegral<O>(6, 16, R, dist1, dist2, rho1, rho2);
  mat(6, 20) = mat(6, 16);
  mat(6, 23) = calculator_.getIntegral<O>(6, 23, R, dist1, dist2, rho1, rho2);
  mat(6, 26) = calculator_.getIntegral<O>(6, 26, R, dist1, dist2, rho1, rho2);
  mat(6, 31) = calculator_.getIntegral<O>(6, 31, R, dist1, dist2, rho1, rho2);
  mat(6, 35) = mat(6, 31);
  mat(6, 38) = calculator_.getIntegral<O>(6, 38, R, dist1, dist2, rho1, rho2);
  mat(6, 39) = mat(6, 38);
  mat(7, 11) = calculator_.getIntegral<O>(7, 11, R, dist1, dist2, rho1, rho2);
  mat(7, 15) = calculator_.getIntegral<O>(7, 15, R, dist1, dist2, rho1, rho2);
  mat(7, 17) = calculator_.getIntegral<O>(7, 17, R, dist1, dist2, rho1, rho2);
  mat(7, 22) = mat(7, 17);
  mat(7, 24) = mat(7, 17);
  mat(7, 27) = calculator_.getIntegral<O>(7, 27, R, dist1, dist2, rho1, rho2);
  mat(7, 33) = calculator_.getIntegral<O>(7, 33, R, dist1, dist2, rho1, rho2);
  mat(7, 37) = mat(7, 33);
  mat(8, 12) = mat(7, 11);
  mat(8, 18) = mat(7, 17);
  mat(8, 19) = mat(7, 15);
  mat(8, 21) = -mat(7, 17);
  mat(8, 25) = mat(7, 17);
  mat(8, 28) = mat(7, 27);
  mat(8, 34) = mat(7, 33);
  mat(8, 36) = -mat(7, 33);
  mat(9, 10) = calculator_.getIntegral<O>(9, 10, R, dist1, dist2, rho1, rho2);
  mat(9, 16) = calculator_.getIntegral<O>(9, 16, R, dist1, dist2, rho1, rho2);
  mat(9, 20) = mat(9, 16);
  mat(9, 23) = calculator_.getIntegral<O>(9, 23, R, dist1, dist2, rho1, rho2);
  mat(9, 26) = calculator_.getIntegral<O>(9, 26, R, dist1, dist2, rho1, rho2);
  mat(9, 31) = calculator_.getIntegral<O>(9, 31, R, dist1, dist2, rho1, rho2);
  mat(9, 35) = mat(9, 31);
  mat(9, 38) = calculator_.getIntegral<O>(9, 38, R, dist1, dist2, rho1, rho2);
  mat(9, 39) = mat(9, 38);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildDPMatrix(double R) {
  mat(10, 2) = calculator_.getIntegral<O>(10, 2, R, dist1, dist2, rho1, rho2);
  mat(10, 5) = mat(10, 2);
  mat(10, 6) = calculator_.getIntegral<O>(10, 6, R, dist1, dist2, rho1, rho2);
  mat(10, 9) = calculator_.getIntegral<O>(10, 9, R, dist1, dist2, rho1, rho2);
  mat(11, 1) = calculator_.getIntegral<O>(11, 1, R, dist1, dist2, rho1, rho2);
  mat(11, 7) = calculator_.getIntegral<O>(11, 7, R, dist1, dist2, rho1, rho2);
  mat(12, 3) = mat(11, 1);
  mat(12, 8) = mat(11, 7);
  mat(13, 2) = calculator_.getIntegral<O>(13, 2, R, dist1, dist2, rho1, rho2);
  mat(13, 5) = -mat(13, 2);
  mat(14, 4) = mat(13, 2);
  mat(15, 1) = calculator_.getIntegral<O>(15, 1, R, dist1, dist2, rho1, rho2);
  mat(15, 7) = calculator_.getIntegral<O>(15, 7, R, dist1, dist2, rho1, rho2);
  mat(16, 2) = calculator_.getIntegral<O>(16, 2, R, dist1, dist2, rho1, rho2);
  mat(16, 5) = mat(16, 2);
  mat(16, 6) = calculator_.getIntegral<O>(16, 6, R, dist1, dist2, rho1, rho2);
  mat(16, 9) = calculator_.getIntegral<O>(16, 9, R, dist1, dist2, rho1, rho2);
  mat(17, 1) = calculator_.getIntegral<O>(17, 1, R, dist1, dist2, rho1, rho2);
  mat(17, 7) = calculator_.getIntegral<O>(17, 7, R, dist1, dist2, rho1, rho2);
  mat(18, 3) = mat(17, 1);
  mat(18, 8) = mat(17, 7);
  mat(19, 3) = mat(15, 1);
  mat(19, 8) = mat(15, 7);
  mat(20, 2) = mat(16, 2);
  mat(20, 5) = mat(16, 2);
  mat(20, 6) = mat(16, 6);
  mat(20, 9) = mat(16, 9);
  mat(21, 3) = -mat(17, 1);
  mat(21, 8) = -mat(17, 7);
  mat(22, 1) = mat(17, 1);
  mat(22, 7) = mat(17, 7);
  mat(23, 2) = calculator_.getIntegral<O>(23, 2, R, dist1, dist2, rho1, rho2);
  mat(23, 5) = mat(23, 2);
  mat(23, 6) = calculator_.getIntegral<O>(23, 6, R, dist1, dist2, rho1, rho2);
  mat(23, 9) = calculator_.getIntegral<O>(23, 9, R, dist1, dist2, rho1, rho2);
  mat(24, 1) = mat(17, 1);
  mat(24, 7) = mat(17, 7);
  mat(25, 3) = mat(17, 1);
  mat(25, 8) = mat(17, 7);
  mat(26, 2) = calculator_.getIntegral<O>(26, 2, R, dist1, dist2, rho1, rho2);
  mat(26, 5) = mat(26, 2);
  mat(26, 6) = calculator_.getIntegral<O>(26, 6, R, dist1, dist2, rho1, rho2);
  mat(26, 9) = calculator_.getIntegral<O>(26, 9, R, dist1, dist2, rho1, rho2);
  mat(27, 1) = calculator_.getIntegral<O>(27, 1, R, dist1, dist2, rho1, rho2);
  mat(27, 7) = calculator_.getIntegral<O>(27, 7, R, dist1, dist2, rho1, rho2);
  mat(28, 3) = mat(27, 1);
  mat(28, 8) = mat(27, 7);
  mat(29, 2) = calculator_.getIntegral<O>(29, 2, R, dist1, dist2, rho1, rho2);
  mat(29, 5) = -mat(29, 2);
  mat(30, 4) = mat(29, 2);
  mat(31, 2) = calculator_.getIntegral<O>(31, 2, R, dist1, dist2, rho1, rho2);
  mat(31, 5) = calculator_.getIntegral<O>(31, 5, R, dist1, dist2, rho1, rho2);
  mat(31, 6) = calculator_.getIntegral<O>(31, 6, R, dist1, dist2, rho1, rho2);
  mat(31, 9) = calculator_.getIntegral<O>(31, 9, R, dist1, dist2, rho1, rho2);
  mat(32, 4) = calculator_.getIntegral<O>(32, 4, R, dist1, dist2, rho1, rho2);
  mat(33, 1) = calculator_.getIntegral<O>(33, 1, R, dist1, dist2, rho1, rho2);
  mat(33, 7) = calculator_.getIntegral<O>(33, 7, R, dist1, dist2, rho1, rho2);
  mat(34, 3) = mat(33, 1);
  mat(34, 8) = mat(33, 7);
  mat(35, 2) = mat(31, 5);
  mat(35, 5) = mat(31, 2);
  mat(35, 6) = mat(31, 6);
  mat(35, 9) = mat(31, 9);
  mat(36, 3) = -mat(33, 1);
  mat(36, 8) = -mat(33, 7);
  mat(37, 1) = mat(33, 1);
  mat(37, 7) = mat(33, 7);
  mat(38, 2) = calculator_.getIntegral<O>(38, 2, R, dist1, dist2, rho1, rho2);
  mat(38, 5) = mat(38, 2);
  mat(38, 6) = calculator_.getIntegral<O>(38, 6, R, dist1, dist2, rho1, rho2);
  mat(38, 9) = calculator_.getIntegral<O>(38, 9, R, dist1, dist2, rho1, rho2);
  mat(39, 2) = mat(38, 2);
  mat(39, 5) = mat(38, 2);
  mat(39, 6) = mat(38, 6);
  mat(39, 9) = mat(38, 9);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildDPMatrixSym(double /*R*/) {
  mat(10, 2) = mat(2, 10);
  mat(10, 5) = mat(5, 10);
  mat(10, 6) = -mat(6, 10);
  mat(10, 9) = mat(9, 10);
  mat(11, 1) = -mat(1, 11);
  mat(11, 7) = mat(7, 11);
  mat(12, 3) = -mat(3, 12);
  mat(12, 8) = mat(8, 12);
  mat(13, 2) = mat(2, 13);
  mat(13, 5) = mat(5, 13);
  mat(14, 4) = mat(4, 14);
  mat(15, 1) = mat(1, 15);
  mat(15, 7) = -mat(7, 15);
  mat(16, 2) = -mat(2, 16);
  mat(16, 5) = -mat(5, 16);
  mat(16, 6) = mat(6, 16);
  mat(16, 9) = -mat(9, 16);
  mat(17, 1) = mat(1, 17);
  mat(17, 7) = -mat(7, 17);
  mat(18, 3) = mat(3, 18);
  mat(18, 8) = -mat(8, 18);
  mat(19, 3) = mat(3, 19);
  mat(19, 8) = -mat(8, 19);
  mat(20, 2) = -mat(2, 20);
  mat(20, 5) = -mat(5, 20);
  mat(20, 6) = mat(6, 20);
  mat(20, 9) = -mat(9, 20);
  mat(21, 3) = mat(3, 21);
  mat(21, 8) = -mat(8, 21);
  mat(22, 1) = mat(1, 22);
  mat(22, 7) = -mat(7, 22);
  mat(23, 2) = -mat(2, 23);
  mat(23, 5) = -mat(5, 23);
  mat(23, 6) = mat(6, 23);
  mat(23, 9) = -mat(9, 23);
  mat(24, 1) = mat(1, 24);
  mat(24, 7) = -mat(7, 24);
  mat(25, 3) = mat(3, 25);
  mat(25, 8) = -mat(8, 25);
  mat(26, 2) = mat(2, 26);
  mat(26, 5) = mat(5, 26);
  mat(26, 6) = -mat(6, 26);
  mat(26, 9) = mat(9, 26);
  mat(27, 1) = -mat(1, 27);
  mat(27, 7) = mat(7, 27);
  mat(28, 3) = -mat(3, 28);
  mat(28, 8) = mat(8, 28);
  mat(29, 2) = mat(2, 29);
  mat(29, 5) = mat(5, 29);
  mat(30, 4) = mat(4, 30);
  mat(31, 2) = mat(2, 31);
  mat(31, 5) = mat(5, 31);
  mat(31, 6) = -mat(6, 31);
  mat(31, 9) = mat(9, 31);
  mat(32, 4) = mat(4, 32);
  mat(33, 1) = -mat(1, 33);
  mat(33, 7) = mat(7, 33);
  mat(34, 3) = -mat(3, 34);
  mat(34, 8) = mat(8, 34);
  mat(35, 2) = mat(2, 35);
  mat(35, 5) = mat(5, 35);
  mat(35, 6) = -mat(6, 35);
  mat(35, 9) = mat(9, 35);
  mat(36, 3) = -mat(3, 36);
  mat(36, 8) = mat(8, 36);
  mat(37, 1) = -mat(1, 37);
  mat(37, 7) = mat(7, 37);
  mat(38, 2) = mat(2, 38);
  mat(38, 5) = mat(5, 38);
  mat(38, 6) = -mat(6, 38);
  mat(38, 9) = mat(9, 38);
  mat(39, 2) = mat(2, 39);
  mat(39, 5) = mat(5, 39);
  mat(39, 6) = -mat(6, 39);
  mat(39, 9) = mat(9, 39);
}

template<Utils::derivOrder O>
void Local2c2eMatrix<O>::buildDDMatrix(double R) {
  mat(10, 10) = calculator_.getIntegral<O>(10, 10, R, dist1, dist2, rho1, rho2);
  mat(10, 16) = calculator_.getIntegral<O>(10, 16, R, dist1, dist2, rho1, rho2);
  mat(10, 20) = mat(10, 16);
  mat(10, 23) = calculator_.getIntegral<O>(10, 23, R, dist1, dist2, rho1, rho2);
  mat(10, 26) = calculator_.getIntegral<O>(10, 26, R, dist1, dist2, rho1, rho2);
  mat(10, 31) = calculator_.getIntegral<O>(10, 31, R, dist1, dist2, rho1, rho2);
  mat(10, 35) = mat(10, 31);
  mat(10, 38) = calculator_.getIntegral<O>(10, 38, R, dist1, dist2, rho1, rho2);
  mat(10, 39) = mat(10, 38);
  mat(11, 11) = calculator_.getIntegral<O>(11, 11, R, dist1, dist2, rho1, rho2);
  mat(11, 15) = calculator_.getIntegral<O>(11, 15, R, dist1, dist2, rho1, rho2);
  mat(11, 17) = calculator_.getIntegral<O>(11, 17, R, dist1, dist2, rho1, rho2);
  mat(11, 22) = mat(11, 17);
  mat(11, 24) = mat(11, 17);
  mat(11, 27) = calculator_.getIntegral<O>(11, 27, R, dist1, dist2, rho1, rho2);
  mat(11, 33) = calculator_.getIntegral<O>(11, 33, R, dist1, dist2, rho1, rho2);
  mat(11, 37) = mat(11, 33);
  mat(12, 12) = mat(11, 11);
  mat(12, 18) = mat(11, 17);
  mat(12, 19) = mat(11, 15);
  mat(12, 21) = -mat(11, 17);
  mat(12, 25) = mat(11, 17);
  mat(12, 28) = mat(11, 27);
  mat(12, 34) = mat(11, 33);
  mat(12, 36) = -mat(11, 33);
  mat(13, 13) = calculator_.getIntegral<O>(13, 13, R, dist1, dist2, rho1, rho2);
  mat(13, 29) = calculator_.getIntegral<O>(13, 29, R, dist1, dist2, rho1, rho2);
  mat(13, 31) = calculator_.getIntegral<O>(13, 31, R, dist1, dist2, rho1, rho2);
  mat(13, 35) = -mat(13, 31);
  mat(14, 14) = mat(13, 13);
  mat(14, 30) = mat(13, 29);
  mat(14, 32) = mat(13, 31);
  mat(15, 11) = calculator_.getIntegral<O>(15, 11, R, dist1, dist2, rho1, rho2);
  mat(15, 15) = calculator_.getIntegral<O>(15, 15, R, dist1, dist2, rho1, rho2);
  mat(15, 17) = calculator_.getIntegral<O>(15, 17, R, dist1, dist2, rho1, rho2);
  mat(15, 22) = mat(15, 17);
  mat(15, 24) = mat(15, 17);
  mat(15, 27) = calculator_.getIntegral<O>(15, 27, R, dist1, dist2, rho1, rho2);
  mat(15, 33) = calculator_.getIntegral<O>(15, 33, R, dist1, dist2, rho1, rho2);
  mat(15, 37) = mat(15, 33);
  mat(16, 10) = calculator_.getIntegral<O>(16, 10, R, dist1, dist2, rho1, rho2);
  mat(16, 16) = calculator_.getIntegral<O>(16, 16, R, dist1, dist2, rho1, rho2);
  mat(16, 20) = mat(16, 16);
  mat(16, 23) = calculator_.getIntegral<O>(16, 23, R, dist1, dist2, rho1, rho2);
  mat(16, 26) = calculator_.getIntegral<O>(16, 26, R, dist1, dist2, rho1, rho2);
  mat(16, 31) = calculator_.getIntegral<O>(16, 31, R, dist1, dist2, rho1, rho2);
  mat(16, 35) = mat(16, 31);
  mat(16, 38) = calculator_.getIntegral<O>(16, 38, R, dist1, dist2, rho1, rho2);
  mat(16, 39) = mat(16, 38);
  mat(17, 11) = calculator_.getIntegral<O>(17, 11, R, dist1, dist2, rho1, rho2);
  mat(17, 15) = mat(15, 17);
  mat(17, 17) = calculator_.getIntegral<O>(17, 17, R, dist1, dist2, rho1, rho2);
  mat(17, 22) = mat(17, 17);
  mat(17, 24) = mat(17, 17);
  mat(17, 27) = -mat(15, 33);
  mat(17, 33) = calculator_.getIntegral<O>(17, 33, R, dist1, dist2, rho1, rho2);
  mat(17, 37) = mat(17, 33);
  mat(18, 12) = mat(17, 11);
  mat(18, 18) = mat(17, 17);
  mat(18, 19) = mat(15, 17);
  mat(18, 21) = -mat(17, 17);
  mat(18, 25) = mat(17, 17);
  mat(18, 28) = -mat(15, 33);
  mat(18, 34) = mat(17, 33);
  mat(18, 36) = -mat(17, 33);
  mat(19, 12) = mat(15, 11);
  mat(19, 18) = mat(15, 17);
  mat(19, 19) = mat(15, 15);
  mat(19, 21) = -mat(15, 17);
  mat(19, 25) = mat(15, 17);
  mat(19, 28) = mat(15, 27);
  mat(19, 34) = mat(15, 33);
  mat(19, 36) = -mat(15, 33);
  mat(20, 10) = mat(16, 10);
  mat(20, 16) = mat(16, 16);
  mat(20, 20) = mat(16, 16);
  mat(20, 23) = mat(16, 23);
  mat(20, 26) = mat(16, 26);
  mat(20, 31) = mat(16, 31);
  mat(20, 35) = mat(16, 31);
  mat(20, 38) = mat(16, 38);
  mat(20, 39) = mat(16, 38);
  mat(21, 12) = -mat(17, 11);
  mat(21, 18) = -mat(17, 17);
  mat(21, 19) = -mat(15, 17);
  mat(21, 21) = mat(17, 17);
  mat(21, 25) = -mat(17, 17);
  mat(21, 28) = mat(15, 33);
  mat(21, 34) = -mat(17, 33);
  mat(21, 36) = mat(17, 33);
  mat(22, 11) = mat(17, 11);
  mat(22, 15) = mat(15, 17);
  mat(22, 17) = mat(17, 17);
  mat(22, 22) = mat(17, 17);
  mat(22, 24) = mat(17, 17);
  mat(22, 27) = -mat(15, 33);
  mat(22, 33) = mat(17, 33);
  mat(22, 37) = mat(17, 33);
  mat(23, 10) = calculator_.getIntegral<O>(23, 10, R, dist1, dist2, rho1, rho2);
  mat(23, 16) = mat(16, 23);
  mat(23, 20) = mat(16, 23);
  mat(23, 23) = calculator_.getIntegral<O>(23, 23, R, dist1, dist2, rho1, rho2);
  mat(23, 26) = calculator_.getIntegral<O>(23, 26, R, dist1, dist2, rho1, rho2);
  mat(23, 31) = calculator_.getIntegral<O>(23, 31, R, dist1, dist2, rho1, rho2);
  mat(23, 35) = mat(23, 31);
  mat(23, 38) = calculator_.getIntegral<O>(23, 38, R, dist1, dist2, rho1, rho2);
  mat(23, 39) = mat(23, 38);
  mat(24, 11) = mat(17, 11);
  mat(24, 15) = mat(15, 17);
  mat(24, 17) = mat(17, 17);
  mat(24, 22) = mat(17, 17);
  mat(24, 24) = mat(17, 17);
  mat(24, 27) = -mat(15, 33);
  mat(24, 33) = mat(17, 33);
  mat(24, 37) = mat(17, 33);
  mat(25, 12) = mat(17, 11);
  mat(25, 18) = mat(17, 17);
  mat(25, 19) = mat(15, 17);
  mat(25, 21) = -mat(17, 17);
  mat(25, 25) = mat(17, 17);
  mat(25, 28) = -mat(15, 33);
  mat(25, 34) = mat(17, 33);
  mat(25, 36) = -mat(17, 33);
  mat(26, 10) = calculator_.getIntegral<O>(26, 10, R, dist1, dist2, rho1, rho2);
  mat(26, 16) = calculator_.getIntegral<O>(26, 16, R, dist1, dist2, rho1, rho2);
  mat(26, 20) = mat(26, 16);
  mat(26, 23) = calculator_.getIntegral<O>(26, 23, R, dist1, dist2, rho1, rho2);
  mat(26, 26) = calculator_.getIntegral<O>(26, 26, R, dist1, dist2, rho1, rho2);
  mat(26, 31) = calculator_.getIntegral<O>(26, 31, R, dist1, dist2, rho1, rho2);
  mat(26, 35) = mat(26, 31);
  mat(26, 38) = calculator_.getIntegral<O>(26, 38, R, dist1, dist2, rho1, rho2);
  mat(26, 39) = mat(26, 38);
  mat(27, 11) = calculator_.getIntegral<O>(27, 11, R, dist1, dist2, rho1, rho2);
  mat(27, 15) = calculator_.getIntegral<O>(27, 15, R, dist1, dist2, rho1, rho2);
  mat(27, 17) = calculator_.getIntegral<O>(27, 17, R, dist1, dist2, rho1, rho2);
  mat(27, 22) = mat(27, 17);
  mat(27, 24) = mat(27, 17);
  mat(27, 27) = calculator_.getIntegral<O>(27, 27, R, dist1, dist2, rho1, rho2);
  mat(27, 33) = calculator_.getIntegral<O>(27, 33, R, dist1, dist2, rho1, rho2);
  mat(27, 37) = mat(27, 33);
  mat(28, 12) = mat(27, 11);
  mat(28, 18) = mat(27, 17);
  mat(28, 19) = mat(27, 15);
  mat(28, 21) = -mat(27, 17);
  mat(28, 25) = mat(27, 17);
  mat(28, 28) = mat(27, 27);
  mat(28, 34) = mat(27, 33);
  mat(28, 36) = -mat(27, 33);
  mat(29, 13) = calculator_.getIntegral<O>(29, 13, R, dist1, dist2, rho1, rho2);
  mat(29, 29) = calculator_.getIntegral<O>(29, 29, R, dist1, dist2, rho1, rho2);
  mat(29, 31) = calculator_.getIntegral<O>(29, 31, R, dist1, dist2, rho1, rho2);
  mat(29, 35) = -mat(29, 31);
  mat(30, 14) = mat(29, 13);
  mat(30, 30) = mat(29, 29);
  mat(30, 32) = mat(29, 31);
  mat(31, 10) = calculator_.getIntegral<O>(31, 10, R, dist1, dist2, rho1, rho2);
  mat(31, 13) = calculator_.getIntegral<O>(31, 13, R, dist1, dist2, rho1, rho2);
  mat(31, 16) = calculator_.getIntegral<O>(31, 16, R, dist1, dist2, rho1, rho2);
  mat(31, 20) = mat(31, 16);
  mat(31, 23) = calculator_.getIntegral<O>(31, 23, R, dist1, dist2, rho1, rho2);
  mat(31, 26) = calculator_.getIntegral<O>(31, 26, R, dist1, dist2, rho1, rho2);
  mat(31, 29) = mat(29, 31);
  mat(31, 31) = calculator_.getIntegral<O>(31, 31, R, dist1, dist2, rho1, rho2);
  mat(31, 35) = calculator_.getIntegral<O>(31, 35, R, dist1, dist2, rho1, rho2);
  mat(31, 38) = calculator_.getIntegral<O>(31, 38, R, dist1, dist2, rho1, rho2);
  mat(31, 39) = mat(31, 38);
  mat(32, 14) = mat(31, 13);
  mat(32, 30) = mat(29, 31);
  mat(32, 32) = calculator_.getIntegral<O>(32, 32, R, dist1, dist2, rho1, rho2);
  mat(33, 11) = calculator_.getIntegral<O>(33, 11, R, dist1, dist2, rho1, rho2);
  mat(33, 15) = -mat(27, 17);
  mat(33, 17) = calculator_.getIntegral<O>(33, 17, R, dist1, dist2, rho1, rho2);
  mat(33, 22) = mat(33, 17);
  mat(33, 24) = mat(33, 17);
  mat(33, 27) = mat(27, 33);
  mat(33, 33) = calculator_.getIntegral<O>(33, 33, R, dist1, dist2, rho1, rho2);
  mat(33, 37) = mat(33, 33);
  mat(34, 12) = mat(33, 11);
  mat(34, 18) = mat(33, 17);
  mat(34, 19) = -mat(27, 17);
  mat(34, 21) = -mat(33, 17);
  mat(34, 25) = mat(33, 17);
  mat(34, 28) = mat(27, 33);
  mat(34, 34) = mat(33, 33);
  mat(34, 36) = -mat(33, 33);
  mat(35, 10) = mat(31, 10);
  mat(35, 13) = -mat(31, 13);
  mat(35, 16) = mat(31, 16);
  mat(35, 20) = mat(31, 16);
  mat(35, 23) = mat(31, 23);
  mat(35, 26) = mat(31, 26);
  mat(35, 29) = -mat(29, 31);
  mat(35, 31) = mat(31, 35);
  mat(35, 35) = mat(31, 31);
  mat(35, 38) = mat(31, 38);
  mat(35, 39) = mat(31, 38);
  mat(36, 12) = -mat(33, 11);
  mat(36, 18) = -mat(33, 17);
  mat(36, 19) = mat(27, 17);
  mat(36, 21) = mat(33, 17);
  mat(36, 25) = -mat(33, 17);
  mat(36, 28) = -mat(27, 33);
  mat(36, 34) = -mat(33, 33);
  mat(36, 36) = mat(33, 33);
  mat(37, 11) = mat(33, 11);
  mat(37, 15) = -mat(27, 17);
  mat(37, 17) = mat(33, 17);
  mat(37, 22) = mat(33, 17);
  mat(37, 24) = mat(33, 17);
  mat(37, 27) = mat(27, 33);
  mat(37, 33) = mat(33, 33);
  mat(37, 37) = mat(33, 33);
  mat(38, 10) = calculator_.getIntegral<O>(38, 10, R, dist1, dist2, rho1, rho2);
  mat(38, 16) = calculator_.getIntegral<O>(38, 16, R, dist1, dist2, rho1, rho2);
  mat(38, 20) = mat(38, 16);
  mat(38, 23) = calculator_.getIntegral<O>(38, 23, R, dist1, dist2, rho1, rho2);
  mat(38, 26) = calculator_.getIntegral<O>(38, 26, R, dist1, dist2, rho1, rho2);
  mat(38, 31) = calculator_.getIntegral<O>(38, 31, R, dist1, dist2, rho1, rho2);
  mat(38, 35) = mat(38, 31);
  mat(38, 38) = calculator_.getIntegral<O>(38, 38, R, dist1, dist2, rho1, rho2);
  mat(38, 39) = mat(38, 38);
  mat(39, 10) = mat(38, 10);
  mat(39, 16) = mat(38, 16);
  mat(39, 20) = mat(38, 16);
  mat(39, 23) = mat(38, 23);
  mat(39, 26) = mat(38, 26);
  mat(39, 31) = mat(38, 31);
  mat(39, 35) = mat(38, 31);
  mat(39, 38) = mat(38, 38);
  mat(39, 39) = mat(38, 38);
}

template class Local2c2eMatrix<Utils::derivOrder::zero>;
template class Local2c2eMatrix<Utils::derivOrder::one>;
template class Local2c2eMatrix<Utils::derivOrder::two>;

} // namespace multipole

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
