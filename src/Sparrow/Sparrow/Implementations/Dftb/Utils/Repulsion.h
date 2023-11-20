/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB_REPULSION_H
#define SPARROW_DFTB_REPULSION_H

#include "DFTBCommon.h"
#include <Utils/Scf/MethodInterfaces/RepulsionCalculator.h>
#include <memory>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace dftb {
class PairwiseRepulsion;

class Repulsion : public Utils::RepulsionCalculator {
 public:
  Repulsion(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
            const DFTBCommon::DiatomicParameterContainer& diatomicParameters);
  ~Repulsion() override;

  using PairRepulsion = std::unique_ptr<dftb::PairwiseRepulsion>;
  using Container = std::vector<std::vector<PairRepulsion>>;

  void initialize() override;
  void calculateRepulsion(Utils::DerivativeOrder order) override;
  double getRepulsionEnergy() const override;
  void addRepulsionDerivatives(Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::First>& derivatives) const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::SecondAtomic>& derivatives) const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::Derivative::SecondFull>& derivatives) const override;

 private:
  template<Utils::Derivative O>
  void addRepulsionDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const;
  void calculatePairRepulsion(int i, int j, Utils::DerivativeOrder order);
  void initializePair(int i, int j);

  int nAtoms_;
  Container pairRepulsions_;
  const DFTBCommon::DiatomicParameterContainer& diatomicParameters_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_REPULSION_H
