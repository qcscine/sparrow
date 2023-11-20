/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DensityGuess.h"
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/DensityMatrix.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

DensityGuess::DensityGuess(const Utils::AtomsOrbitalsIndexes& aoIndexes, const std::vector<double>& coreCharges,
                           const int& nElectrons)
  : aoIndexes_(aoIndexes), coreCharges_(coreCharges), nElectrons_(nElectrons) {
}

Utils::DensityMatrix DensityGuess::calculateGuess() const {
  auto nAOs = aoIndexes_.getNAtomicOrbitals();
  auto nAtoms = aoIndexes_.getNAtoms();

  // Set the density matrix in order for the atomic charges to be zero
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nAOs, nAOs);
  for (int i = 0; i < nAtoms; ++i) {
    auto aoIndex = aoIndexes_.getFirstOrbitalIndex(i);
    auto nAOsForAtom = aoIndexes_.getNOrbitals(i);
    for (int a = 0; a < nAOsForAtom; ++a) {
      P(aoIndex + a, aoIndex + a) = coreCharges_[i] / nAOsForAtom;
    }
  }

  Utils::DensityMatrix d;
  d.setDensity(std::move(P), nElectrons_);
  return d;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
