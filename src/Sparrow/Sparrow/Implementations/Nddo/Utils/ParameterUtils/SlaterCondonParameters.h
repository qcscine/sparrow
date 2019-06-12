/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_SLATERCONDONPARAMETERS_H
#define SPARROW_SLATERCONDONPARAMETERS_H

#include "Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterSlaterIntegral.h"
#include <Utils/Geometry/ElementTypes.h>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace nddo {

enum sc_t {
  F0ss,
  F0pp,
  F0dd,
  F0sp,
  F0sd,
  F0pd,
  F2pp,
  F2dd,
  F2pd,
  F4dd,
  G1sp,
  G1pd,
  G2sd,
  G3pd,
  R1sppd,
  R2sdpp,
  R2sddd
};

/*!
 * This class is the container for the Slater-Condon parameters.
 */
class SlaterCondonParameters {
 public:
  SlaterCondonParameters();
  void setElement(Utils::ElementType element);
  void setExponents(double es, double ep = 1.0, double ed = 1.0); // Default value in case atom doesn't have those
                                                                  // orbitals
  void set(sc_t type, double value);
  double get(sc_t type) const;
  void calculate();

 private:
  unsigned int getIndex(sc_t type) const {
    return static_cast<unsigned int>(type);
  }
  double getUlValue(unsigned int l, unsigned int a, unsigned int b, unsigned int c, unsigned int d);

  Utils::ElementType element_;
  unsigned int n_[3];
  double exp_[3];
  std::vector<double> parameters_;
  std::vector<bool> alreadyGiven_;
  OneCenterSlaterIntegral Ul_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_SLATERCONDONPARAMETERS_H
