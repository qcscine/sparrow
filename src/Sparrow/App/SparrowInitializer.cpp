/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SparrowInitializer.h"
#include <Core/Exceptions.h>
#include <Core/ModuleManager.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/dll/shared_library.hpp>
#include <boost/filesystem.hpp>

namespace Scine {
namespace Sparrow {

void SparrowInitializer::initialize() {
  namespace fs = boost::filesystem;
  const std::string moduleName = "sparrow.module" + boost::dll::shared_library::suffix().string();

  fs::path appPath = fs::canonical(boost::dll::program_location());
  fs::path basePath = appPath.parent_path();

  auto& manager = Core::ModuleManager::getInstance();
  if (!manager.moduleLoaded("Sparrow")) {
    std::vector<fs::path> trialPaths{appPath.parent_path() / moduleName, basePath / moduleName, basePath / "lib" / moduleName};

    for (const auto& trialPath : trialPaths) {
      if (fs::exists(trialPath) && fs::is_regular_file(trialPath)) {
        manager.load(trialPath);
        break;
      }
    }

    if (!manager.moduleLoaded("Sparrow")) {
      throw std::runtime_error(
          "Sparrow module could not be found. Please add the path to the sparrow module to SCINE_MODULE_PATH.");
    }
  }
}

} // namespace Sparrow
} // namespace Scine
