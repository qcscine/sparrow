/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "GTOOverlapMatrixBlock.h"
#include "GeneralTypes.h"
#include <Utils/DataStructures/GtoExpansion.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {
using namespace GeneralTypes;

using std::exp;
using std::sqrt;

template<Utils::derivOrder O>
GTOOverlapMatrixBlock<O>::GTOOverlapMatrixBlock() : sqrt3(sqrt(3)), pi(4.0 * atan(1)), nullValue(constant1D<O>(0)) {
  momenta[0] = AngularMomentum(0, 0, 0);
  momenta[1] = AngularMomentum(1, 0, 0);
  momenta[2] = AngularMomentum(0, 1, 0);
  momenta[3] = AngularMomentum(0, 0, 1);
  momenta[4] = AngularMomentum(1, 0, 1);
  momenta[5] = AngularMomentum(0, 1, 1);
  momenta[6] = AngularMomentum(1, 1, 0);
  momenta[7] = AngularMomentum(0, 0, 2); // z2
  momenta[8] = AngularMomentum(0, 2, 0); // y2
  momenta[9] = AngularMomentum(2, 0, 0); // x2

  AOIndexes[0] = static_cast<int>(orb_t::s);
  AOIndexes[1] = static_cast<int>(orb_t::x) - 1;
  AOIndexes[2] = static_cast<int>(orb_t::y) - 1;
  AOIndexes[3] = static_cast<int>(orb_t::z) - 1;
  AOIndexes[4] = static_cast<int>(orb_t::xz) - 4;
  AOIndexes[5] = static_cast<int>(orb_t::yz) - 4;
  AOIndexes[6] = static_cast<int>(orb_t::xy) - 4;
  AOIndexes[7] = static_cast<int>(orb_t::z2) - 4;
  AOIndexes[8] = static_cast<int>(orb_t::x2y2) - 4;
  AOIndexes[9] = 0;
}

template<Utils::derivOrder O>
Eigen::Matrix<typename GTOOverlapMatrixBlock<O>::Value3D, Eigen::Dynamic, Eigen::Dynamic>
GTOOverlapMatrixBlock<O>::getMatrixBlock(const Utils::GtoExpansion& gA, const Utils::GtoExpansion& gB,
                                         const Eigen::Vector3d& Rab) {
  startGTFA_ = (gA.nAOs() == 1) ? 0 : (gA.nAOs() == 3) ? 1 : 4;
  startGTFB_ = (gB.nAOs() == 1) ? 0 : (gB.nAOs() == 3) ? 1 : 4;
  nGTFsA_ = (gA.nAOs() == 5) ? 6 : gA.nAOs();
  nGTFsB_ = (gB.nAOs() == 5) ? 6 : gB.nAOs();

  result.resize(nGTFsA_, nGTFsB_);
  result.setConstant(constant3D<O>(0.0));

  for (int i = 0; i < gA.size(); i++) {
    for (int j = 0; j < gB.size(); j++) {
      addGTFContribution(i, j, gA, gB, Rab);
    }
  }

  return computeBlock(gA, gB);
}
template<Utils::derivOrder O>
void GTOOverlapMatrixBlock<O>::addGTFContribution(int GTFA, int GTFB, const Utils::GtoExpansion& gA,
                                                  const Utils::GtoExpansion& gB, const Eigen::Vector3d& Rab) {
  // Get GTO exponents and compute their sum
  double expA = gA.getExponent(GTFA);
  double expB = gB.getExponent(GTFB);
  double expSum = expA + expB;

  // Calculate the center of gravity of the two gaussian orbitals
  double x = Rab[0], y = Rab[1], z = Rab[2];
  double factorForP = expB / expSum;
  double fac = -factorForP * expA;
  setRArray(x, y, z, expSum);
  p[0] = r[0] * factorForP;
  p[1] = r[1] * factorForP;
  p[2] = r[2] * factorForP;

  // Compute preexponential factor and multiply it with the contraction coefficients and normalization factors
  KD = getKBase(x, y, z, fac) * exp((x * x + y * y + z * z) * fac) *
       (gA.getNormalizedCoefficient(GTFA) * gB.getNormalizedCoefficient(GTFB));

  // S_00 is known: done in setRArray
  SDo[1][0][0] = SDo[0][0][0];
  SDo[2][0][0] = SDo[0][0][0];

  for (int d = 0; d < 3; d++) { // Loop over directions (x,y,z)
    // Distances between center of gravity and atoms
    // Dist1 = p[d];
    Dist2 = p[d] - r[d];

    // There are two relations for the Obara-Saika construction of So values:
    // S_i,j = X_PA * S_i-1,j + 1/(2expSum) * ((i-1)*S_i-2,j + j*S_i-1,j-1)
    // S_i,j = X_PB * S_i,j-1 + 1/(2expSum) * (i*S_i-1,j-1 + (j-1)*S_i,j-2)
    // for either k or l equal to zero, only one of them can work
    for (int k = 0; k <= gA.angularMomentum(); k++) {
      for (int l = 0; l <= gB.angularMomentum(); l++) {
        if (k != 0 || l != 0) {
          SDo[d][k][l] = nullValue;
          if (l == 0) { // use first formula
            if (k > 1)
              SDo[d][k][l] = SDo[d][k - 2][l] * ((k - 1.0) / (2 * expSum));
            SDo[d][k][l] += SDo[d][k - 1][l] * p[d]; // Dist1 = p[d]
          }
          else { // i.e. l>0, use second formula
            if (l > 1)
              SDo[d][k][l] = SDo[d][k][l - 2] * ((l - 1.0) / (2 * expSum));
            if (k > 0)
              SDo[d][k][l] += SDo[d][k - 1][l - 1] * (k / (2 * expSum));
            SDo[d][k][l] += SDo[d][k][l - 1] * Dist2;
          }
        }
      } // End loop l
    }   // End loop k
  }     // End loop d

  // Combine directions and add to S matrix...
  for (int k = 0; k < nGTFsA_; k++) {   // Loop over AOs of orbitals on A
    for (int l = 0; l < nGTFsB_; l++) { // Loop over AOs of orbitals on B
      const Value1D X = SDo[0][momenta[k + startGTFA_].x][momenta[l + startGTFB_].x];
      const Value1D Y = SDo[1][momenta[k + startGTFA_].y][momenta[l + startGTFB_].y];
      const Value1D Z = SDo[2][momenta[k + startGTFA_].z][momenta[l + startGTFB_].z];
      contributionD = KD * directionsProduct(X, Y, Z);
      result(k, l) += contributionD;
    } // End loop l
  }   // End loop k
}

template<Utils::derivOrder O>
Eigen::Matrix<typename GTOOverlapMatrixBlock<O>::Value3D, Eigen::Dynamic, Eigen::Dynamic>
GTOOverlapMatrixBlock<O>::computeBlock(const Utils::GtoExpansion& gA, const Utils::GtoExpansion& gB) {
  // go from 6 d-type GTFs back to 5 d-type atomic orbitals:
  // dz2 = 1/2*(gzz-gxx-gyy)
  // dx2y2 = sqrt(3/4)*(gxx-gyy)

  if (nGTFsA_ == 6 || nGTFsB_ == 6) {
    Eigen::Matrix<Value3D, Eigen::Dynamic, Eigen::Dynamic> solidHarmonicsResult(gA.nAOs(), gB.nAOs());
    if (nGTFsA_ == 6) {
      for (int l = 0; l < nGTFsB_; l++) {
        result(3, l) = (result(3, l) - 0.5 * result(4, l) - 0.5 * result(5, l)) / sqrt3;
        result(4, l) = 0.5 * (result(5, l) - result(4, l));
      }
    }
    if (nGTFsB_ == 6) {
      for (int k = 0; k < nGTFsA_; k++) {
        result(k, 3) = (result(k, 3) - 0.5 * result(k, 4) - 0.5 * result(k, 5)) / sqrt3;
        result(k, 4) = 0.5 * (result(k, 5) - result(k, 4));
      }
    }
    for (int l = 0; l < gB.nAOs(); ++l)
      for (int k = 0; k < gA.nAOs(); ++k)
        solidHarmonicsResult(AOIndexes[startGTFA_ + k], AOIndexes[startGTFB_ + l]) = result(k, l);
    return solidHarmonicsResult;
  }
  else {
    return result;
  }
}

template<>
inline GTOOverlapMatrixBlock<Utils::derivOrder::zero>::Value3D
GTOOverlapMatrixBlock<Utils::derivOrder::zero>::directionsProduct(const Value1D& X, const Value1D& Y, const Value1D& Z) {
  return X * Y * Z;
}
template<>
inline GTOOverlapMatrixBlock<Utils::derivOrder::one>::Value3D
GTOOverlapMatrixBlock<Utils::derivOrder::one>::directionsProduct(const Value1D& X, const Value1D& Y, const Value1D& Z) {
  return {X.value() * Y.value() * Z.value(), X.derivative() * Y.value() * Z.value(),
          Y.derivative() * X.value() * Z.value(), Z.derivative() * X.value() * Y.value()};
}
template<>
inline GTOOverlapMatrixBlock<Utils::derivOrder::two>::Value3D
GTOOverlapMatrixBlock<Utils::derivOrder::two>::directionsProduct(const Value1D& X, const Value1D& Y, const Value1D& Z) {
  return {X.value() * Y.value() * Z.value(),  X.first() * Y.value() * Z.value(),  X.value() * Y.first() * Z.value(),
          X.value() * Y.value() * Z.first(),  X.second() * Y.value() * Z.value(), X.value() * Y.second() * Z.value(),
          X.value() * Y.value() * Z.second(), X.first() * Y.first() * Z.value(),  X.first() * Y.value() * Z.first(),
          X.value() * Y.first() * Z.first()};
}

template<>
inline GTOOverlapMatrixBlock<Utils::derivOrder::zero>::Value3D
GTOOverlapMatrixBlock<Utils::derivOrder::zero>::getKBase(double /*x*/, double /*y*/, double /*z*/, double /*fac*/) {
  return 1;
}
template<>
inline GTOOverlapMatrixBlock<Utils::derivOrder::one>::Value3D
GTOOverlapMatrixBlock<Utils::derivOrder::one>::getKBase(double x, double y, double z, double fac) {
  return {1, 2 * x * fac, 2 * y * fac, 2 * z * fac};
}
template<>
inline GTOOverlapMatrixBlock<Utils::derivOrder::two>::Value3D
GTOOverlapMatrixBlock<Utils::derivOrder::two>::getKBase(double x, double y, double z, double fac) {
  double xfac2 = 2 * x * fac;
  double yfac2 = 2 * y * fac;
  double zfac2 = 2 * z * fac;
  return {1,
          xfac2,
          yfac2,
          zfac2,
          xfac2 * xfac2 + 2 * fac,
          yfac2 * yfac2 + 2 * fac,
          zfac2 * zfac2 + 2 * fac,
          xfac2 * yfac2,
          xfac2 * zfac2,
          yfac2 * zfac2};
}

template<Utils::derivOrder O>
inline void GTOOverlapMatrixBlock<O>::setRArray(double x, double y, double z, double expSum) {
  SDo[0][0][0] = constant1D<O>(sqrt(pi / expSum));
  r[0] = getFromFull<O>(x, 1, 0);
  r[1] = getFromFull<O>(y, 1, 0);
  r[2] = getFromFull<O>(z, 1, 0);
}

template class GTOOverlapMatrixBlock<Utils::derivOrder::zero>;
template class GTOOverlapMatrixBlock<Utils::derivOrder::one>;
template class GTOOverlapMatrixBlock<Utils::derivOrder::two>;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
