/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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
  void calculateRepulsion(Utils::derivOrder order) override;
  double getRepulsionEnergy() const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const override;
  void addRepulsionDerivatives(
      Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const override;

 private:
  template<Utils::derivativeType O>
  void addRepulsionDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const;
  void calculatePairRepulsion(int i, int j, Utils::derivOrder order);
  void initializePair(int i, int j);

  int nAtoms_;
  Container pairRepulsions_;
  const DFTBCommon::DiatomicParameterContainer& diatomicParameters_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTB_REPULSION_H