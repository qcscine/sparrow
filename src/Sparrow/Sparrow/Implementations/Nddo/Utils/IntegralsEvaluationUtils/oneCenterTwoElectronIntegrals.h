/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_ONECENTERTWOELECTRONINTEGRALS_H
#define SPARROW_ONECENTERTWOELECTRONINTEGRALS_H

#include "OneCenterTwoElectronCalculator.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PrincipalQuantumNumbers.h>
#include <Utils/Geometry/ElementTypes.h>
#include <exception>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace nddo {

class NoElementSetException : public std::exception {};
class NoExponentsSetException : public std::exception {};
class SlaterCondonParameters;
enum class BasisFunctions;

/*!
 * This class calculates all the unique one-center two-electron
 * integrals for a given element.
 */
class OneCenterTwoElectronIntegrals {
 public:
  using orb_index_t = int;

  void setElement(Utils::ElementType e, BasisFunctions basisFunctions = BasisFunctions::spd);
  void calculateIntegrals();
  int getNumberIntegrals() const;

  void set(orb_index_t a, orb_index_t b, orb_index_t c, orb_index_t d, double value);
  double get(orb_index_t a, orb_index_t b, orb_index_t c, orb_index_t d) const;
  double get(int index) const;

  void setSlaterCondonParameters(const SlaterCondonParameters* slaterCondonParameters);

 private:
  void initialize(BasisFunctions basisFunctions);
  void setExchangePIntegral();
  void calculateIntegralsFromExponents();
  bool needsIntegralsFromSlaterCondonParameters();

  Utils::ElementType element_;
  int nIntegrals_;
  std::vector<double> integrals_;
  std::vector<bool> alreadyGiven_;
  const SlaterCondonParameters* scPar_{nullptr};
  OneCenterTwoElectronCalculator calculator_;
  bool elementHasBeenSet_{false};
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_ONECENTERTWOELECTRONINTEGRALS_H
