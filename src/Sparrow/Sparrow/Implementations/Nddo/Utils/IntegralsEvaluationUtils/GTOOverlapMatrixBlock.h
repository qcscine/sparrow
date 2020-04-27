/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_GTOOVERLAPMATRIXBLOCK_H
#define SPARROW_GTOOVERLAPMATRIXBLOCK_H

#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Eigen/Core>

namespace Scine {

namespace Utils {
class MatrixWithDerivatives;
class GtoExpansion;
} // namespace Utils

namespace Sparrow {

namespace nddo {

struct AngularMomentum {
  explicit AngularMomentum(int lx = 0, int ly = 0, int lz = 0) : x(lx), y(ly), z(lz) {
  }
  int x;
  int y;
  int z;
};

/**
 * @brief This class calculates the overlap matrix block and its derivatives for two groups
 *        of orbitals on different atoms, each group of which shares the
 *        same angular momentum. F.i. s-s, s-p, p-d, d-d, ... according to the Obara-Saika method.
 */

template<Utils::derivOrder O>
class GTOOverlapMatrixBlock {
 public:
  using Value3D = Utils::AutomaticDifferentiation::Value3DType<O>;
  using Value1D = Utils::AutomaticDifferentiation::Value1DType<O>;
  //! @brief Constructor initializing the angular momenta and indices of the orbitals used for the calculation.
  GTOOverlapMatrixBlock();
  //! @brief Getter for the matrix block calculated.
  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic>
  getMatrixBlock(const Utils::GtoExpansion& gA, const Utils::GtoExpansion& gB, const Eigen::Vector3d& Rab);

 private:
  void addGTFContribution(int GTFA, int GTFB, const Utils::GtoExpansion& gA, const Utils::GtoExpansion& gB,
                          const Eigen::Vector3d& Rab);
  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic> computeBlock(const Utils::GtoExpansion& gA,
                                                                      const Utils::GtoExpansion& gB);
  Value3D directionsProduct(const Value1D& X, const Value1D& Y, const Value1D& Z);
  Value3D getKBase(double x, double y, double z, double fac);
  void setRArray(double x, double y, double z, double expSum);

  const double sqrt3, pi;

  const Value1D nullValue;
  int startGTFA_, startGTFB_;
  int nGTFsA_, nGTFsB_;
  Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic> result;
  AngularMomentum momenta[10];
  int AOIndexes[10];
  Value1D r[3], p[3];
  Value1D SDo[3][3][3];
  Value3D KD, contributionD;
  Value1D Dist1, Dist2;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_GTOOVERLAPMATRIXBLOCK_H
