/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include <Core/ModuleManager.h>
#include <gmock/gmock.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/dll/shared_library.hpp>

class SparrowModuleLoad : public ::testing::Environment {
 public:
  ~SparrowModuleLoad() override {
  }

  void SetUp() override {
    auto& moduleManager = Scine::Core::ModuleManager::getInstance();
    if (!moduleManager.moduleLoaded("Sparrow")) {
      const auto programPath = boost::dll::program_location();
      auto libPath = programPath.parent_path() / "sparrow.module";
      libPath += boost::dll::shared_library::suffix();
      moduleManager.load(libPath);
    }
  }
};

/* Add SparrowModuleLoad to global testing environments at global variable init
 * time (before main()), see gtest docs regarding global setup and teardown
 */
::testing::Environment* const moduleLoadEnvironment = ::testing::AddGlobalTestEnvironment(new SparrowModuleLoad);
