/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTBMethodWrapper.h"
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Sparrow/StatesHandling/DFTBStatesHandler.h>

namespace Scine {
namespace Sparrow {

DFTBMethodWrapper::DFTBMethodWrapper() {
  this->statesHandler_ = std::make_unique<DFTBStatesHandler>(*this);
}

DFTBMethodWrapper::~DFTBMethodWrapper() = default;

Utils::PropertyList DFTBMethodWrapper::possibleProperties() const {
  return Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Dipole |
         Utils::Property::DipoleGradient | Utils::Property::DipoleMatrixMO;
}

Utils::Results DFTBMethodWrapper::assembleResults(const std::string& description) const {
  bool MODipoleMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::DipoleMatrixMO) &&
                                possibleProperties().containsSubSet(Utils::Property::DipoleMatrixMO);
  auto results = GenericMethodWrapper::assembleResults(description);

  auto dipoleEvaluationCoordinate = Utils::Position::Zero();

  if (MODipoleMatrixRequired) {
    if (!dipoleMatrixCalculator_->isValid()) {
      dipoleMatrixCalculator_->fillDipoleMatrix(dipoleEvaluationCoordinate);
    }
    results.setMODipoleMatrix(dipoleMatrixCalculator_->getMODipoleMatrix());
  }
  return results;
}

} // namespace Sparrow
} // namespace Scine
