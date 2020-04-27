/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_NDDODENSITYGUESS_H
#define SPARROW_NDDODENSITYGUESS_H

#include <Utils/Scf/MethodInterfaces/DensityMatrixGuessCalculator.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace Utils {
class OverlapCalculator;
}
namespace Sparrow {

namespace nddo {
class ElementParameters;

/*!
 * Implementation of DensityMatrixGuessCalculator for the NDDO methods.
 */
class NDDODensityGuess : public Utils::DensityMatrixGuessCalculator {
 public:
  NDDODensityGuess(const Utils::ElementTypeCollection& elements, const ElementParameters& elementParameters,
                   Utils::OverlapCalculator& overlapCalculator, const int& nElectrons, const int& nAOs);

  Utils::DensityMatrix calculateGuess() const override;
  void setNElectrons(int nElectrons) final;

 private:
  const Utils::ElementTypeCollection& elements_;
  const ElementParameters& elementParameters_;
  Utils::OverlapCalculator& overlapCalculator_;
  int nElectrons_;
  const int& nAOs_;
};

} // namespace nddo

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_PM6DENSITYGUESS_H
