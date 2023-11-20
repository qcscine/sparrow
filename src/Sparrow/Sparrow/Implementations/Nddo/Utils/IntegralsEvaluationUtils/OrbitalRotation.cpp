/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "OrbitalRotation.h"
#include "GeneralTypes.h"
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace nddo {
using namespace GeneralTypes;

using std::sqrt;

template<Utils::DerivativeOrder O>
OrbitalRotation<O>::OrbitalRotation(int l1, int l2)
  : needsP(l1 > 0 || l2 > 0),
    needsD(l1 > 1 || l2 > 1),
    sqrt3(std::sqrt(3)),
    nullDeriv(constant3D<O>(0.0)),
    oneDeriv(constant3D<O>(1.0)) {
}

template<>
double OrbitalRotation<Utils::DerivativeOrder::Zero>::setUpOrderDependentValues(const Eigen::Vector3d& Rab) {
  X = Rab(0);
  Y = Rab(1);
  Z = Rab(2);
  R = Rab.norm();
  return R;
}

template<>
double OrbitalRotation<Utils::DerivativeOrder::One>::setUpOrderDependentValues(const Eigen::Vector3d& Rab) {
  X = First3D(Rab(0), 1, 0, 0);
  Y = First3D(Rab(1), 0, 1, 0);
  Z = First3D(Rab(2), 0, 0, 1);
  R = sqrt(X * X + Y * Y + Z * Z);
  return R.value();
}
template<>
double OrbitalRotation<Utils::DerivativeOrder::Two>::setUpOrderDependentValues(const Eigen::Vector3d& Rab) {
  X = Second3D(Rab(0), 1, 0, 0, 0, 0, 0, 0, 0, 0);
  Y = Second3D(Rab(1), 0, 1, 0, 0, 0, 0, 0, 0, 0);
  Z = Second3D(Rab(2), 0, 0, 1, 0, 0, 0, 0, 0, 0);
  R = sqrt(X * X + Y * Y + Z * Z);
  return R.value();
}

template<Utils::DerivativeOrder O>
inline bool isOnlyZComponent(const Utils::AutomaticDifferentiation::Value3DType<O>& planeProjection) {
  return planeProjection == 0;
}
template<>
inline bool isOnlyZComponent<Utils::DerivativeOrder::Zero>(
    const Utils::AutomaticDifferentiation::Value3DType<Utils::DerivativeOrder::Zero>& planeProjection) {
  return planeProjection == 0;
}
template<>
inline bool isOnlyZComponent<Utils::DerivativeOrder::One>(
    const Utils::AutomaticDifferentiation::Value3DType<Utils::DerivativeOrder::One>& planeProjection) {
  return planeProjection.value() == 0;
}
template<>
inline bool isOnlyZComponent<Utils::DerivativeOrder::Two>(
    const Utils::AutomaticDifferentiation::Value3DType<Utils::DerivativeOrder::Two>& planeProjection) {
  return planeProjection.value() == 0;
}

template<Utils::DerivativeOrder O>
double OrbitalRotation<O>::setVector(const Eigen::Vector3d& Rab) {
  double Rabs = setUpOrderDependentValues(Rab);
  Utils::AutomaticDifferentiation::Value3DType<O> planeProjection;

  if (needsP) {
    planeProjection = nullDeriv;
    if (Rab.head(2).norm() != 0) {
      planeProjection = sqrt(X * X + Y * Y);
    }
    if (!isOnlyZComponent<O>(planeProjection)) {
      cosa = (X) / planeProjection;
      sina = (Y) / planeProjection;
    }
    else {
      cosa = oneDeriv;
      sina = nullDeriv;
    }
    cosb = (Z / R);
    sinb = planeProjection / R;

    if (needsD) {
      sin2b = sinb * sinb;
      sin2a = sina * sina;
      cos2b = cosb * cosb;
      cos2a = cosa * cosa;
      saca = sina * cosa;
      sbcb = sinb * cosb;
      dc2am1 = 2 * cos2a - oneDeriv;
    }
  }

  return Rabs;
}
template<Utils::DerivativeOrder O>
void OrbitalRotation<O>::fillRotVector() {
  rotVector[static_cast<int>(rotationOrbitalPair::s_s)] = oneDeriv;

  if (needsP) {
    rotVector[static_cast<int>(rotationOrbitalPair::x_x)] = cosa * cosb;
    rotVector[static_cast<int>(rotationOrbitalPair::y_x)] = sina * cosb;
    rotVector[static_cast<int>(rotationOrbitalPair::z_x)] = -sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::x_y)] = -sina;
    rotVector[static_cast<int>(rotationOrbitalPair::y_y)] = cosa;
    // rotVector[static_cast<int>(rotationOrbitalPair::z_y)] = nullDeriv; //not needed
    rotVector[static_cast<int>(rotationOrbitalPair::x_z)] = cosa * sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::y_z)] = sina * sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::z_z)] = cosb;
  }

  if (needsD) {
    rotVector[static_cast<int>(rotationOrbitalPair::x2y2_x2y2)] = dc2am1 * cos2b + 0.5 * dc2am1 * sin2b;
    rotVector[static_cast<int>(rotationOrbitalPair::x2y2_xz)] = dc2am1 * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::x2y2_z2)] = 0.5 * sqrt3 * dc2am1 * sin2b;
    rotVector[static_cast<int>(rotationOrbitalPair::x2y2_yz)] = -2 * saca * sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::x2y2_xy)] = -2 * saca * cosb;
    rotVector[static_cast<int>(rotationOrbitalPair::xz_x2y2)] = -cosa * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::xz_xz)] = cosa * (2 * cos2b - oneDeriv);
    rotVector[static_cast<int>(rotationOrbitalPair::xz_z2)] = sqrt3 * cosa * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::xz_yz)] = -sina * cosb;
    rotVector[static_cast<int>(rotationOrbitalPair::xz_xy)] = sina * sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::z2_x2y2)] = 0.5 * sqrt3 * sin2b;
    rotVector[static_cast<int>(rotationOrbitalPair::z2_xz)] = -sqrt3 * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::z2_z2)] = cos2b - 0.5 * sin2b;
    // rotVector[static_cast<int>(rotationOrbitalPair::z2_yz)] = nullDeriv; // not needed
    // rotVector[static_cast<int>(rotationOrbitalPair::z2_xy)] = nullDeriv; // not needed
    rotVector[static_cast<int>(rotationOrbitalPair::yz_x2y2)] = -sina * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::yz_xz)] = sina * (2 * cos2b - oneDeriv);
    rotVector[static_cast<int>(rotationOrbitalPair::yz_z2)] = sqrt3 * sina * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::yz_yz)] = cosa * cosb;
    rotVector[static_cast<int>(rotationOrbitalPair::yz_xy)] = -cosa * sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::xy_x2y2)] = 2 * saca * cos2b + saca * sin2b;
    rotVector[static_cast<int>(rotationOrbitalPair::xy_xz)] = 2 * saca * sbcb;
    rotVector[static_cast<int>(rotationOrbitalPair::xy_z2)] = sqrt3 * saca * sin2b;
    rotVector[static_cast<int>(rotationOrbitalPair::xy_yz)] = dc2am1 * sinb;
    rotVector[static_cast<int>(rotationOrbitalPair::xy_xy)] = dc2am1 * cosb;
  }
}

template class OrbitalRotation<Utils::DerivativeOrder::Zero>;
template class OrbitalRotation<Utils::DerivativeOrder::One>;
template class OrbitalRotation<Utils::DerivativeOrder::Two>;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
