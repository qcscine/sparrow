/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DFTBDipoleMomentCalculator.h"
#include <Sparrow/Implementations/Dftb/Dftb0/DFTB0.h>
#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
#include <Sparrow/Implementations/Dftb/Dftb3/DFTB3.h>

namespace Scine {
namespace Sparrow {

template<class DFTBMethod>
DFTBDipoleMomentCalculator<DFTBMethod>::~DFTBDipoleMomentCalculator() = default;

template<class DFTBMethod>
DFTBDipoleMomentCalculator<DFTBMethod>::DFTBDipoleMomentCalculator(const DFTBMethod& method) : method_(method) {
}

template<class DFTBMethod>
Eigen::RowVector3d Scine::Sparrow::DFTBDipoleMomentCalculator<DFTBMethod>::calculate() const {
  const auto& atomicCharges = method_.getAtomicCharges();
  const auto& positions = method_.getPositions();
  Eigen::Vector3d dipole = Eigen::Vector3d::Zero();

  assert(atomicCharges.size() == positions.rows() && "Not same amount of atomic charges and positions.");

  for (int atom = 0; atom < atomicCharges.size(); ++atom) {
    dipole += positions.row(atom) * atomicCharges[atom];
  }

  return dipole;
}

template class DFTBDipoleMomentCalculator<dftb::DFTB0>;
template class DFTBDipoleMomentCalculator<dftb::DFTB2>;
template class DFTBDipoleMomentCalculator<dftb::DFTB3>;

} // namespace Sparrow
} // namespace Scine
