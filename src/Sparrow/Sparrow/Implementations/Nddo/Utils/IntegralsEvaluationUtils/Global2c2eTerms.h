/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_GLOBAL2C2ETERMS_H
#define SPARROW_GLOBAL2C2ETERMS_H

#include "OrbitalRotation.h"
#include "TwoElectronIntegralIndexes.h"
#include <array>
#include <list>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

struct RotationTerm {
  RotationTerm(unsigned int p1, unsigned int p2, GeneralTypes::rotationOrbitalPair f1, GeneralTypes::rotationOrbitalPair f2,
               GeneralTypes::rotationOrbitalPair f3, GeneralTypes::rotationOrbitalPair f4)
    : f1_(f1), f2_(f2), f3_(f3), f4_(f4), pair1_(p1), pair2_(p2) {
  }

  GeneralTypes::rotationOrbitalPair f1_, f2_, f3_, f4_; // NB: or store directly the uint ?
  unsigned int pair1_, pair2_;
};

class Global2c2eTerms {
 public:
  using RotationTerms = std::list<RotationTerm>;
  using RotationTermsArray = std::vector<std::vector<RotationTerms>>;

  using orb_index_t = int;
  using orbPair_index_t = int;

  const RotationTerms& getTermList(orbPair_index_t op1, orbPair_index_t op2) const {
    static RotationTermsArray expressions = createRotationTerms();
    return expressions[op1][op2];
  }

 private:
  static RotationTermsArray createRotationTerms();
  static void createTerm(RotationTermsArray& expressions, std::array<int, 8> i);
  static bool compatibleOrbitals(int a, int b);
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_GLOBAL2C2ETERMS_H