/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "SKAtom.h"

namespace Scine {
namespace Sparrow {

namespace dftb {

SKAtom::SKAtom(Utils::ElementType el) : allowsSpin(false), allowsDFTB3(false) {
  if (el <= Utils::ElementType::He) {
    highestOrbital = s;
    nAOs = 1;
  }
  else if (el <= Utils::ElementType::Ne) {
    highestOrbital = p;
    nAOs = 4;
  }
  else {
    highestOrbital = d;
    nAOs = 9;
  }
}

void SKAtom::setEnergies(double es, double ep, double ed) {
  Es = es;
  Ep = ep;
  Ed = ed;
}

void SKAtom::setOccupations(int fs, int fp, int fd) {
  Fs = fs;
  Fp = fp;
  Fd = fd;
  totalOccupation = Fs + Fp + Fd;
}

void SKAtom::setHubbardParameter(double us, double up, double ud) {
  Us = us;
  Up = up;
  Ud = ud;
}

double SKAtom::getOrbitalEnergy(int orbital) {
  if (orbital == 0)
    return Es;
  if (orbital < 4)
    return Ep;
  return Ed;
}

double SKAtom::getEnergy() {
  return Es * Fs + Ep * Fp + Ed * Fd;
}

double SKAtom::getHubbardParameter() {
  return Us;
  // Remark: all three Us, Up and Ud are given in the Slater-Koster files, but they are identical.
}

void SKAtom::setSpinConstants(double spinConstants[][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      sc[i][j] = spinConstants[i][j];
  allowsSpin = true;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
