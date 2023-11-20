/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_SKATOM_H
#define SPARROW_SKATOM_H

#include <Utils/Geometry/ElementInfo.h>
#include <algorithm>
#include <array>
#include <vector>

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
  void setSpinConstants(std::array<std::array<double, 3>, 3> arr);
  bool hasSpinConstants() const {
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
  /**
   * @brief Function returning the spin constants (magnetic Hubbard)
   * Spin constants are stored as follows:
   * sc[0][0] : W_s
   * sc[0][1] : W_{s,p}
   * sc[0][2] : W_{s,d}
   * sc[1][0] : W_{p,s}
   * sc[1][1] : W_p
   * sc[1][2] : W_{p,d}
   * sc[2][0] : W_{d,s}
   * sc[2][1] : W_{d,p}
   * sc[2][2] : W_d
   */
  double getSpinConstant(int i, int j) {
    return sc[i][j];
  }
  double getAtomicResolvedSpinConstant() const;
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
  std::array<std::array<double, 3>, 3> sc; // Array containing the spin constants; first index is for angular momentum
                                           // of first orbital, second index for second one.
  bool allowsSpin;

  // For DFTB3
  double hubbardDerivative;
  bool allowsDFTB3;
};

inline double SKAtom::getAtomicResolvedSpinConstant() const {
  if (!hasSpinConstants())
    throw std::runtime_error("No spin constant available for element " + Utils::ElementInfo::symbol(element_));
  if (highestOrbital == orbital::s) { // s orbital is HOMO
    return sc[0][0];                  // return element W_{ss}
  }
  else if (highestOrbital == orbital::p) { // p orbital is HOMO
    return sc[1][1];                       // return element W_{pp}
  }
  else {             // d orbital is HOMO
    return sc[2][2]; // return element W_{dd}
  }
}
} // namespace dftb
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_SKATOM_H
