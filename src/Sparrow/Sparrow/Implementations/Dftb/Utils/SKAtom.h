/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_SKATOM_H
#define SPARROW_SKATOM_H

#include <Utils/Geometry/ElementTypes.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

class SKAtom {
 public:
  enum orbital { s, p, d };

  explicit SKAtom(Utils::ElementType el);
  void setEnergies(double es, double ep, double ed);
  void setOccupations(int fs, int fp, int fd);
  void setHubbardParameter(double us, double up, double ud);
  void setSpinConstants(double spinConstants[][3]);
  bool hasSpinConstants() {
    return allowsSpin;
  }
  void setHubbardDerivative(double hubbard) {
    hubbardDerivative = hubbard;
    allowsDFTB3 = true;
  }
  double getHubbardDerivative() const {
    return hubbardDerivative;
  }
  bool hasHubbardDerivative() const {
    return allowsDFTB3;
  }
  int getnAOs() const {
    return nAOs;
  }
  int getOccupation() const {
    return totalOccupation;
  }
  double getHubbardParameter() const;
  double getOrbitalEnergy(int orbital) const;
  double getEnergy() const;
  double getSpinConstant(int i, int j) {
    return sc[i][j];
  }
  orbital getHighestOrbital() const {
    return highestOrbital;
  }

  bool operator==(const SKAtom& rhs) const {
    return element_ == rhs.element_;
  }

 private:
  Utils::ElementType element_;
  int nAOs;
  orbital highestOrbital;
  double Es, Ep, Ed;               // Energies of orbitals
  int totalOccupation, Fs, Fp, Fd; // Occupations of orbitals
  double Us, Up, Ud;               // Hubbard parameters of orbitals

  // For SDFTB
  double sc[3][3]; // Array containing the spin constants; first index is for angular momentum of first orbital, second
                   // index for second one.
  bool allowsSpin;

  // For DFTB3
  double hubbardDerivative;
  bool allowsDFTB3;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_SKATOM_H
