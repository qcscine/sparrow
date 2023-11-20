/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_CISAOFOCKBUILDERFACTORY_H
#define SPARROW_CISAOFOCKBUILDERFACTORY_H

#include "CISMatrixAOFockBuilder.h"
namespace Scine {
namespace Sparrow {

template<Utils::Reference restrictedness>
class CISMatrixAOFockBuilderFactory {
 public:
  static std::shared_ptr<CISMatrixAOFockBuilderBase<restrictedness>>
  createAOFockBuilder(const Utils::SpinTransition spinBlock, CISData cisData, const ExcitedStatesParam& excitedStatesParam) {
    switch (spinBlock) {
      case Utils::SpinTransition::Singlet:
        return std::make_shared<CISMatrixAOFockBuilder<restrictedness, Utils::SpinTransition::Singlet>>(
            std::move(cisData), excitedStatesParam);
      case Utils::SpinTransition::Triplet:
        return std::make_shared<CISMatrixAOFockBuilder<restrictedness, Utils::SpinTransition::Triplet>>(
            std::move(cisData), excitedStatesParam);
      default:
        throw std::runtime_error("Invalid spin-block argument in CISMatrixAOFockBuilderFactory.");
    }
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_CISAOFOCKBUILDERFACTORY_H
