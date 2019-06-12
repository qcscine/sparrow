/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_DFTBDIPOLEMOMENTCALCULATOR_H
#define SPARROW_DFTBDIPOLEMOMENTCALCULATOR_H

#include <Sparrow/Implementations/DipoleMomentCalculator.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace Sparrow {

/**
 * @brief This class calculates the electrical dipole moment for the DFTB methods.
 * Right now it calculates the dipole through a Mulliken population analysis. New methods
 * can then be implemented if deemed necessary.
 * In particular, with this approximation transition dipoles from orbitals
 * on the same atom can be vastly underestimated, which is of importance for TD-DFTB.
 * @tparam DFTBMethod One of the DFTB methods type, i.e. DFTB0, DFTB2, DFTB3.
 */
template<class DFTBMethod>
class DFTBDipoleMomentCalculator : public DipoleMomentCalculator {
 public:
  explicit DFTBDipoleMomentCalculator(const DFTBMethod& method);
  ~DFTBDipoleMomentCalculator() final;
  /**
   * @brief Calculates the dipole with a Mulliken population analysis.
   */
  Eigen::RowVector3d calculate() const final;

 private:
  const DFTBMethod& method_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DFTBDIPOLEMOMENTCALCULATOR_H
