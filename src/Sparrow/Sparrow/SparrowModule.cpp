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
#include <Sparrow/Implementations/Dftb/Dftb0/Wrapper/DFTB0MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb2/Wrapper/DFTB2MethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb3/Wrapper/DFTB3MethodWrapper.h>
#include <Sparrow/Implementations/MoldenFileGenerator.h>
/* External Includes */
#include <Core/DerivedModule.h>
#include <Core/Exceptions.h>
#include <Core/Interfaces/WavefunctionOutputGenerator.h>
#include <Utils/Settings.h>

namespace Scine {
namespace Sparrow {

using InterfaceModelMap = boost::mpl::map<
    boost::mpl::pair<Core::Calculator, boost::mpl::vector<PM6MethodWrapper, AM1MethodWrapper, RM1MethodWrapper, PM3MethodWrapper, MNDOMethodWrapper,
                                                          DFTB0MethodWrapper, DFTB2MethodWrapper, DFTB3MethodWrapper>>,
    boost::mpl::pair<Core::WavefunctionOutputGenerator, boost::mpl::vector<PM6MethodWrapper, AM1MethodWrapper, RM1MethodWrapper, PM3MethodWrapper, MNDOMethodWrapper,
                                                                           DFTB0MethodWrapper, DFTB2MethodWrapper, DFTB3MethodWrapper>>>;

std::string SparrowModule::name() const noexcept {
  return "Sparrow";
};

boost::any SparrowModule::get(const std::string& interface, const std::string& model) const {
  boost::any resolved = Core::DerivedModule::resolve<InterfaceModelMap>(interface, model);

  // Throw an exception if we could not match an interface or model
  if (resolved.empty()) {
    throw Core::ClassNotImplementedError();
  }

  return resolved;
};

bool SparrowModule::has(const std::string& interface, const std::string& model) const noexcept {
  return Core::DerivedModule::has<InterfaceModelMap>(interface, model);
}

std::vector<std::string> SparrowModule::announceInterfaces() const noexcept {
  return Core::DerivedModule::announceInterfaces<InterfaceModelMap>();
}

std::vector<std::string> SparrowModule::announceModels(const std::string& interface) const noexcept {
  return Core::DerivedModule::announceModels<InterfaceModelMap>(interface);
}

std::shared_ptr<Core::Module> SparrowModule::make() {
  return std::make_shared<SparrowModule>();
}

std::vector<std::shared_ptr<Core::Module>> moduleFactory() {
  return {SparrowModule::make()};
}
} /* namespace Sparrow */
} /* namespace Scine */
