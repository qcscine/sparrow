/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SlaterCondonParameters.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PrincipalQuantumNumbers.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

SlaterCondonParameters::SlaterCondonParameters()
  : parameters_(std::vector<double>(17, 0.0)), alreadyGiven_(std::vector<bool>(17, false)) {
}

void SlaterCondonParameters::setElement(Utils::ElementType element) {
  element_ = element;
  n_[0] = PM6Elements::getQuantumNumberForSOrbital(element_);
  n_[1] = PM6Elements::getQuantumNumberForPOrbital(element_);
  n_[2] = PM6Elements::getQuantumNumberForDOrbital(element_);
}

void SlaterCondonParameters::setExponents(double es, double ep, double ed) {
  exp_[0] = es;
  exp_[1] = ep;
  exp_[2] = ed;
}

void SlaterCondonParameters::set(sc_t type, double value) {
  parameters_[getIndex(type)] = value;
  alreadyGiven_[getIndex(type)] = true;
}

double SlaterCondonParameters::get(sc_t type) const {
  unsigned int index = getIndex(type);
  return parameters_[index];
}

void SlaterCondonParameters::calculate() {
  std::vector<double> calculatedPar(17);
  calculatedPar[getIndex(F0ss)] = getUlValue(0, 0, 0, 0, 0);
  calculatedPar[getIndex(F0pp)] = getUlValue(0, 1, 1, 1, 1);
  calculatedPar[getIndex(F0dd)] = getUlValue(0, 2, 2, 2, 2);
  calculatedPar[getIndex(F0sp)] = getUlValue(0, 0, 0, 1, 1);
  calculatedPar[getIndex(F0sd)] = getUlValue(0, 0, 0, 2, 2);
  calculatedPar[getIndex(F0pd)] = getUlValue(0, 1, 1, 2, 2);
  calculatedPar[getIndex(F2pp)] = getUlValue(2, 1, 1, 1, 1);
  calculatedPar[getIndex(F2dd)] = getUlValue(2, 2, 2, 2, 2);
  calculatedPar[getIndex(F2pd)] = getUlValue(2, 1, 1, 2, 2);
  calculatedPar[getIndex(F4dd)] = getUlValue(4, 2, 2, 2, 2);
  calculatedPar[getIndex(G1sp)] = getUlValue(1, 0, 1, 1, 0);
  calculatedPar[getIndex(G1pd)] = getUlValue(1, 1, 2, 2, 1);
  calculatedPar[getIndex(G2sd)] = getUlValue(2, 0, 2, 2, 0);
  calculatedPar[getIndex(G3pd)] = getUlValue(3, 1, 2, 2, 1);
  calculatedPar[getIndex(R2sddd)] = getUlValue(2, 0, 2, 2, 2);
  calculatedPar[getIndex(R2sdpp)] = getUlValue(2, 0, 2, 1, 1);
  calculatedPar[getIndex(R1sppd)] = getUlValue(1, 0, 1, 1, 2);

  for (int i = 0; i < 17; i++) {
    if (!alreadyGiven_[i])
      parameters_[i] = calculatedPar[i];
  }
}

double SlaterCondonParameters::getUlValue(unsigned int l, unsigned int a, unsigned int b, unsigned int c, unsigned int d) {
  unsigned int ang[4] = {a, b, c, d};
  unsigned int n[4];
  double exp[4];
  for (int i = 0; i < 4; i++) {
    n[i] = n_[ang[i]];
    exp[i] = exp_[ang[i]];
  }

  Ul_.setPrincipal(n[0], n[1], n[2], n[3]);
  Ul_.setAngular(ang[0], ang[1], ang[2], ang[3]);
  Ul_.setExponents(exp[0], exp[1], exp[2], exp[3]);
  return Ul_.calculate(l);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
