/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

template<Utils::derivOrder O>
OrbitalRotation<O>::OrbitalRotation(int l1, int l2)
  : needsP(l1 > 0 || l2 > 0),
    needsD(l1 > 1 || l2 > 1),
    sqrt3(std::sqrt(3)),
    nullDeriv(constant3D<O>(0.0)),
    oneDeriv(constant3D<O>(1.0)) {
}

template<>
double OrbitalRotation<Utils::derivOrder::zero>::setUpOrderDependentValues(const Eigen::Vector3d& Rab) {
  X = Rab(0);
  Y = Rab(1);
  Z = Rab(2);
  R = Rab.norm();
  return R;
}

template<>
double OrbitalRotation<Utils::derivOrder::one>::setUpOrderDependentValues(const Eigen::Vector3d& Rab) {
  X = First3D(Rab(0), 1, 0, 0);
  Y = First3D(Rab(1), 0, 1, 0);
  Z = First3D(Rab(2), 0, 0, 1);
  R = sqrt(X * X + Y * Y + Z * Z);
  return R.value();
}
template<>
double OrbitalRotation<Utils::derivOrder::two>::setUpOrderDependentValues(const Eigen::Vector3d& Rab) {
  X = Second3D(Rab(0), 1, 0, 0, 0, 0, 0, 0, 0, 0);
  Y = Second3D(Rab(1), 0, 1, 0, 0, 0, 0, 0, 0, 0);
  Z = Second3D(Rab(2), 0, 0, 1, 0, 0, 0, 0, 0, 0);
  R = sqrt(X * X + Y * Y + Z * Z);
  return R.value();
}

template<>
bool OrbitalRotation<Utils::derivOrder::zero>::onlyZComponent() const {
  return planeProj == 0;
}
template<>
bool OrbitalRotation<Utils::derivOrder::one>::onlyZComponent() const {
  return planeProj.value() == 0;
}
template<>
bool OrbitalRotation<Utils::derivOrder::two>::onlyZComponent() const {
  return planeProj.value() == 0;
}

template<Utils::derivOrder O>
double OrbitalRotation<O>::setVector(const Eigen::Vector3d& Rab) {
  double Rabs = setUpOrderDependentValues(Rab);

  if (needsP) {
    planeProj = nullDeriv;
    if (Rab(0) * Rab(0) + Rab(1) * Rab(1) != 0)
      planeProj = sqrt(X * X + Y * Y);

    if (!onlyZComponent()) {
      cosa = (X) / planeProj;
      sina = (Y) / planeProj;
    }
    else {
      cosa = oneDeriv;
      sina = nullDeriv;
    }
    cosb = (Z / R);
    sinb = planeProj / R;

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
template<Utils::derivOrder O>
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

template class OrbitalRotation<Utils::derivOrder::zero>;
template class OrbitalRotation<Utils::derivOrder::one>;
template class OrbitalRotation<Utils::derivOrder::two>;

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
