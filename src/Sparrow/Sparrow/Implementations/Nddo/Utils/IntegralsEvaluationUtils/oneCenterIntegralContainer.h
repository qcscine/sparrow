/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ONECENTERINTEGRALCONTAINER_H
#define SPARROW_ONECENTERINTEGRALCONTAINER_H

#include <Utils/Geometry/ElementInfo.h>
#include <memory>
#include <vector>

namespace Scine {

namespace Utils {
enum class ElementType : unsigned;
}

namespace Sparrow {

namespace nddo {
class OneCenterTwoElectronIntegrals;

/*!
 * This class contains smart pointers to one-center two-electron
 * matrices for the atoms in the AtomCollection.
 * Since they are equal for an element type, they are calculated only
 * once per element type.
 */
class OneCenterIntegralContainer {
 public:
  using integralMatrix_t = std::unique_ptr<OneCenterTwoElectronIntegrals>;
  using Container = std::vector<integralMatrix_t>;

  OneCenterIntegralContainer();
  ~OneCenterIntegralContainer();

  /*! Clear all the contained objects (replaces the pointers by nullptr and deletes the objects). */
  void clear();

  void set(Utils::ElementType e, integralMatrix_t mat);
  const OneCenterTwoElectronIntegrals& get(Utils::ElementType e) const;

 private:
  Container matrices_;
};

inline const OneCenterTwoElectronIntegrals& OneCenterIntegralContainer::get(Utils::ElementType e) const {
  return *matrices_[Utils::ElementInfo::Z(e)];
}

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ONECENTERINTEGRALCONTAINER_H
