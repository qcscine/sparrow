/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef LOCAL2C2EMATRIX_H
#define LOCAL2C2EMATRIX_H

#include "Local2c2eIntegralCalculator.h"
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/*!
 * This class creates the local two-center two-electron matrix for
 * an atom pair.
 */

template<Utils::derivOrder O>
class Local2c2eMatrix {
 public:
  Local2c2eMatrix(int l1, int l2, const ChargeSeparationParameter& D1, const ChargeSeparationParameter& D2,
                  const KlopmanParameter& r1, const KlopmanParameter& r2);
  void setSymmetric(bool sym) {
    sameElement_ = sym;
  }
  void calculate(double R);
  /*! @brief Calculates the two-center two-electron matrix for two identical elements */
  void calculateSym(double R);
  /*! @brief Calculates the two-center two-electron matrix for two different elements */
  void calculateAsym(double R);
  const Utils::AutomaticDifferentiation::Value1DType<O>& operator()(unsigned int i, unsigned int j) const {
    return mat(i, j);
  }

 private:
  void buildSSMatrix(double R);
  void buildSPMatrix(double R);
  void buildPSMatrix(double R);
  void buildPPMatrix(double R);
  void buildSDMatrix(double R);
  void buildDSMatrix(double R);
  void buildPDMatrix(double R);
  void buildDPMatrix(double R);
  void buildDDMatrix(double R);
  void buildPPMatrixSym(double R);
  void buildDDMatrixSym(double R);
  void buildPSMatrixSym(double R);
  void buildDSMatrixSym(double R);
  void buildDPMatrixSym(double R);

  const int l1_, l2_;
  bool sameElement_;
  const ChargeSeparationParameter &dist1, &dist2;
  const KlopmanParameter &rho1, &rho2;
  Local2c2eIntegralCalculator calculator_;
  Eigen::Matrix<Utils::AutomaticDifferentiation::Value1DType<O>, Eigen::Dynamic, Eigen::Dynamic> mat;
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // LOCAL2C2EMATRIX_H
