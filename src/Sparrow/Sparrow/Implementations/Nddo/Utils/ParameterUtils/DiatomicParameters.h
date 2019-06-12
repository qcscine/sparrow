/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DIATOMICPARAMETERS_H
#define SPARROW_DIATOMICPARAMETERS_H

#include <Utils/Geometry/ElementTypes.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

/*!
 * Class for the storage of pairwise parameters in the
 * PM6 method. (only those needed at runtime)
 */

class DiatomicParameters {
 public:
  DiatomicParameters() = default;
  virtual ~DiatomicParameters() = default;

  bool isValid() const {
    return e1_ != Utils::ElementType::none && e2_ != Utils::ElementType::none;
  }

  void setFirstElement(Utils::ElementType e) {
    e1_ = e;
  }
  void setSecondElement(Utils::ElementType e) {
    e2_ = e;
  }
  Utils::ElementType firstElement() const {
    return e1_;
  }
  Utils::ElementType secondElement() const {
    return e2_;
  }

 private:
  Utils::ElementType e1_{Utils::ElementType::none}, e2_{Utils::ElementType::none};
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DIATOMICPARAMETERS_H
