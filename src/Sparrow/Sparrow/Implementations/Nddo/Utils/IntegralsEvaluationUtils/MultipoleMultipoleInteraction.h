/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MULTIPOLEMULTIPOLEINTERACTION_H
#define SPARROW_MULTIPOLEMULTIPOLEINTERACTION_H

#include "MultipoleMultipoleTerm.h"
#include <Utils/Math/AutomaticDifferentiation/AutomaticDifferentiationHelpers.h>
#include <list>

namespace Scine {
namespace Sparrow {

namespace nddo {

namespace multipole {

/**
 * @brief This header-only class performs the actual calculation of the multipole-multipole interaction.
 *
 * First all the charge-charge configurations between two multipoles (MultipoleMultipoleTerm)
 * are inferred and stored in a list, then they are calculated by calling MultipoleMultipoleTerm::calculate(...)
 */

class MultipoleMultipoleInteraction {
 public:
  /**
   * @brief calculates the complete interaction between two multipoles whose single interactions have been stored in
   *        the private member std::list<MultipoleMultipoleTerms> terms_
   * @tparam O order of the desired derivative, for the value this is 0.
   * @param R the distance between the centers of the multipoles
   * @param D1 the charge separation on the first multipole
   * @param d2 the charge separation on the second multipole
   * @param squaredRhos the Klopman--Ohno parameters to compute the repulsion in the Klopman approximation
   * @return an Utils::AutomaticDifferentiation::Value1D<O>, i.e. a collection of all the derivative orders up to the
   * O-th for a 1-dimensional object
   */
  template<Utils::derivOrder O>
  Utils::AutomaticDifferentiation::Value1DType<O> calculate(double R, double D1, double d2, double squaredRhos) const {
    auto sum = Utils::AutomaticDifferentiation::constant1D<O>(0);
    for (const auto& t : terms_) {
      sum += t.calculate<O>(R, D1, d2, squaredRhos);
    }
    return sum;
  }

  /**
   * @brief add an interaction between two charges in the multipoles to the total interaction
   * @param m a MultipoleMultipoleTerm, i.e. a configuration of two charges in two multipoles to be calulated
   */
  void add(const MultipoleMultipoleTerm& m) {
    terms_.push_back(m);
  }
  //! @brief gets the current size of the std::list<MultipoleMultipoleTerm> containing the charge configurations
  unsigned int size() const {
    return static_cast<unsigned int>(terms_.size());
  }

 private:
  std::list<MultipoleMultipoleTerm> terms_;
};

} // namespace multipole

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_MULTIPOLEMULTIPOLEINTERACTION_H
