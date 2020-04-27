/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_DENSITYGUESS_H
#define SPARROW_DFTB_DENSITYGUESS_H

#include <Utils/Scf/MethodInterfaces/DensityMatrixGuessCalculator.h>
#include <vector>

namespace Scine {
namespace Utils {
class AtomsOrbitalsIndexes;
}
namespace Sparrow {

namespace dftb {

class DensityGuess : public Utils::DensityMatrixGuessCalculator {
 public:
  DensityGuess(const Utils::AtomsOrbitalsIndexes& aoIndexes, const std::vector<double>& coreCharges, const unsigned& nElectrons);

  Utils::DensityMatrix calculateGuess() const override;

  void setNElectrons(int nElectrons) final;

 private:
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;
  const std::vector<double>& coreCharges_;
  int nElectrons_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_DENSITYGUESS_H