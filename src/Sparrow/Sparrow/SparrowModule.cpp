/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

/* Internal Includes */
#include "SparrowModule.h"
#include "Sparrow/Implementations/Nddo/Am1/Wrapper/AM1TypeMethodWrapper.h"
#include "Sparrow/Implementations/Nddo/Mndo/Wrapper/MNDOMethodWrapper.h"
#include "Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h"
/* External Includes */
#include <Core/DerivedModule.h>
#include <Core/Exceptions.h>
#include <Core/Interfaces/Calculator.h>
#include <Sparrow/Implementations/Dftb/Dftb0/Wrapper/DFTB0MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb3/Wrapper/DFTB3MethodWrapper.h>
#include <Utils/CalculatorBasics/StatesHandler.h>
#include <Utils/Settings.h>
#include <algorithm>
#include <memory>
#include <string>

namespace Scine {
namespace Sparrow {

SparrowModule::SparrowModule() noexcept
  : calculator{PM6MethodWrapper::model,   AM1MethodWrapper::model,  RM1MethodWrapper::model,
               PM3MethodWrapper::model,   MNDOMethodWrapper::model, DFTB0MethodWrapper::model,
               DFTB2MethodWrapper::model, DFTB3MethodWrapper::model} {
  Core::DerivedModule::checkInvariants(*this);
}

std::string SparrowModule::name() const noexcept {
  return "Sparrow";
};

boost::any SparrowModule::get(const std::string& interface, const std::string& model) const {
  using CalculatorPtr = std::shared_ptr<Core::Calculator>;
  auto calculatorName = model;
  std::transform(calculatorName.begin(), calculatorName.end(), calculatorName.begin(), ::toupper);
  if (interface == Core::Calculator::interface) {
    if (calculatorName == PM6MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<PM6MethodWrapper>());
    }
    else if (calculatorName == AM1MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<AM1MethodWrapper>());
    }
    else if (calculatorName == RM1MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<RM1MethodWrapper>());
    }
    else if (calculatorName == PM3MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<PM3MethodWrapper>());
    }
    else if (calculatorName == MNDOMethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<MNDOMethodWrapper>());
    }
    else if (calculatorName == DFTB0MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<DFTB0MethodWrapper>());
    }
    else if (calculatorName == DFTB2MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<DFTB2MethodWrapper>());
    }
    else if (calculatorName == DFTB3MethodWrapper::model) {
      return static_cast<CalculatorPtr>(std::make_shared<DFTB3MethodWrapper>());
    }
  }
  // Throw an exception if we cannot match a model
  throw Core::ClassNotImplementedError();
};

bool SparrowModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Core::DerivedModule::has(interface, model, *this);
}

std::vector<std::string> SparrowModule::announceInterfaces() const noexcept {
  return Core::DerivedModule::announceInterfaces(*this);
}

std::vector<std::string> SparrowModule::announceModels(const std::string& concept) const noexcept {
  return Core::DerivedModule::announceModels(concept, *this);
}

std::shared_ptr<Core::Module> SparrowModule::make() {
  return std::make_shared<SparrowModule>();
}

} /* namespace Sparrow */
} /* namespace Scine */
