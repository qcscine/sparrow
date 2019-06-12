/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ORBITALROTATION_H
#define SPARROW_ORBITALROTATION_H

#include "GeneralTypes.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

namespace nddo {

template<Utils::derivOrder O>
class OrbitalRotation {
 public:
  OrbitalRotation(int l1, int l2);
  double setVector(const Eigen::Vector3d& Rab);
  void evaluate() {
    fillRotVector();
  }
  const Utils::AutomaticDifferentiation::Value3DType<O>& getRotationCoefficient(GeneralTypes::rotationOrbitalPair p) const {
    return rotVector[static_cast<int>(p)];
  }

 private:
  double setUpOrderDependentValues(const Eigen::Vector3d& Rab);
  bool onlyZComponent() const;
  void fillRotVector();

  Utils::AutomaticDifferentiation::Value3DType<O> rotVector[35];

  const bool needsP, needsD;
  const double sqrt3;
  const Utils::AutomaticDifferentiation::Value3DType<O> nullDeriv, oneDeriv;
  Utils::AutomaticDifferentiation::Value3DType<O> X, Y, Z, R;
  Utils::AutomaticDifferentiation::Value3DType<O> planeProj, cosa, cosb, sina, sinb, sin2b, sin2a, cos2b, cos2a, saca,
      sbcb, dc2am1;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // ORBITALROTATION_H