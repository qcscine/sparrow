/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTBMethodWrapper.h"
#include <Sparrow/Implementations/Dftb/Utils/DipoleUtils/DFTBDipoleMatrixCalculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>

namespace Scine {
namespace Sparrow {

DFTBMethodWrapper::DFTBMethodWrapper() = default;

DFTBMethodWrapper::~DFTBMethodWrapper() = default;

Utils::PropertyList DFTBMethodWrapper::possibleProperties() const {
  auto properties = GenericMethodWrapper::possibleProperties();
  properties.addProperty(Utils::Property::DipoleGradient);
  properties.addProperty(Utils::Property::DipoleMatrixMO);

  return properties;
}

void DFTBMethodWrapper::assembleResults(const std::string& description) {
  bool MODipoleMatrixRequired = requiredProperties_.containsSubSet(Utils::Property::DipoleMatrixMO) &&
                                possibleProperties().containsSubSet(Utils::Property::DipoleMatrixMO);
  GenericMethodWrapper::assembleResults(description);

  auto dipoleEvaluationCoordinate = Utils::Position::Zero();

  if (MODipoleMatrixRequired) {
    if (!dipoleMatrixCalculator_->isValid()) {
      dipoleMatrixCalculator_->fillDipoleMatrix(dipoleEvaluationCoordinate);
    }
    results_.set<Utils::Property::DipoleMatrixMO>(dipoleMatrixCalculator_->getMODipoleMatrix());
  }
}

bool DFTBMethodWrapper::getZPVEInclusion() const {
  return false;
}

} // namespace Sparrow
} // namespace Scine
