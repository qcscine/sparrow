/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "oneCenterTwoElectronIntegrals.h"
#include "Sparrow/Implementations/Nddo/Utils/ParameterUtils/SlaterCondonParameters.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PrincipalQuantumNumbers.h>
#include <cassert>

namespace Scine {
namespace Sparrow {

namespace nddo {

void OneCenterTwoElectronIntegrals::setElement(Utils::ElementType e, BasisFunctions basisFunctions) {
  element_ = e;
  initialize(basisFunctions);
  elementHasBeenSet_ = true;
}

int OneCenterTwoElectronIntegrals::getNumberIntegrals() const {
  return nIntegrals_;
}

void OneCenterTwoElectronIntegrals::initialize(BasisFunctions basisFunctions) {
  nIntegrals_ = PM6Elements::getNumberOneCenterTwoElectronIntegrals(element_, basisFunctions);
  integrals_.resize(nIntegrals_);
  alreadyGiven_.resize(nIntegrals_, false);
  nddo::OneCenterTwoElectronCalculator::setIndexes();
}

void OneCenterTwoElectronIntegrals::calculateIntegrals() {
  setExchangePIntegral();

  if (needsIntegralsFromSlaterCondonParameters())
    calculateIntegralsFromExponents();
}

bool OneCenterTwoElectronIntegrals::needsIntegralsFromSlaterCondonParameters() {
  bool exponentsNeeded = false;
  for (int i = 0; i < nIntegrals_; i++) {
    if (!alreadyGiven_[i])
      exponentsNeeded = true;
  }
  return exponentsNeeded;
}

void OneCenterTwoElectronIntegrals::calculateIntegralsFromExponents() {
  for (int i = 0; i < nIntegrals_; i++) {
    if (!alreadyGiven_[i]) {
      assert(scPar_ && "Slater-Condon parameters are needed for the calculation of integrals and have not been set.");
      integrals_[i] = calculator_.getIntegral(i, scPar_);
    }
  }
}

void OneCenterTwoElectronIntegrals::set(orb_index_t a, orb_index_t b, orb_index_t c, orb_index_t d, double value) {
  int index = OneCenterTwoElectronCalculator::getIndex(a, b, c, d);
  integrals_[index] = value;
  alreadyGiven_[index] = true;
}

void OneCenterTwoElectronIntegrals::setExchangePIntegral() {
  //<pp'|pp'> can be expressed by <pp|pp> and <pp|p'p'>:
  // http://openmopac.net/manual/1c2e.html
  if (nIntegrals_ >= 6) {
    if (alreadyGiven_[nddo::OneCenterTwoElectronCalculator::getIndex(1, 1, 1, 1)] &&
        alreadyGiven_[nddo::OneCenterTwoElectronCalculator::getIndex(1, 1, 2, 2)])
      set(1, 2, 1, 2, 0.5 * (integrals_[3] - integrals_[4]));
  }
}

double OneCenterTwoElectronIntegrals::get(orb_index_t a, orb_index_t b, orb_index_t c, orb_index_t d) const {
  return get(nddo::OneCenterTwoElectronCalculator::getIndex(a, b, c, d));
}

double OneCenterTwoElectronIntegrals::get(int index) const {
  if (index >= nIntegrals_)
    return 0;

  return integrals_[index];
}

void OneCenterTwoElectronIntegrals::setSlaterCondonParameters(const SlaterCondonParameters* slaterCondonParameters) {
  scPar_ = slaterCondonParameters;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
