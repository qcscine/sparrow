/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
  DensityGuess(const Utils::AtomsOrbitalsIndexes& aoIndexes, const std::vector<double>& coreCharges, const int& nElectrons);

  Utils::DensityMatrix calculateGuess() const override;

 private:
  const Utils::AtomsOrbitalsIndexes& aoIndexes_;
  const std::vector<double>& coreCharges_;
  const int& nElectrons_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_DENSITYGUESS_H
