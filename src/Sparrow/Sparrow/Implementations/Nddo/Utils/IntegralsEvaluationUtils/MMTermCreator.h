/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MULTIPOLEMULTIPOLETERMCREATOR_H
#define SPARROW_MULTIPOLEMULTIPOLETERMCREATOR_H

#include "MultipoleChargePair.h"
#include "MultipoleMultipoleInteraction.h"
#include <list>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

class MMTermCreator {
 public:
  void add(const MultipoleChargePair& cp);
  MultipoleMultipoleInteraction computeList() const;

  unsigned int size() const {
    return static_cast<unsigned int>(terms_.size());
  }

 private:
  bool equivalent(const MultipoleChargePair& p1, const MultipoleChargePair& p2);
  double fromEnum(ChargeDistance d) const;
  std::list<MultipoleChargePair> terms_;
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MULTIPOLEMULTIPOLETERMCREATOR_H
