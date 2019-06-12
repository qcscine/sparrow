/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_PM6DIATOMICPARAMETERS_H
#define SPARROW_PM6DIATOMICPARAMETERS_H

#include "DiatomicParameters.h"

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * Class for the storage of pairwise parameters in the
 * PM6 method. (only those needed at runtime)
 */

class PM6DiatomicParameters : public DiatomicParameters {
 public:
  PM6DiatomicParameters(Utils::ElementType e1 = Utils::ElementType::none, Utils::ElementType e2 = Utils::ElementType::none) {
    setFirstElement(e1);
    setSecondElement(e2);
  };

  void setAlpha(double a) {
    alpha_ = a;
  } // set a in bohr^-1; bohr^-2 for O-H and N-H!
  double alpha() const {
    return alpha_;
  }

  void setX(double x) {
    x_ = x;
  }
  double x() const {
    return x_;
  }

 private:
  double alpha_{0.0};
  double x_{0.0};
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_PM6DIATOMICPARAMETERS_H
