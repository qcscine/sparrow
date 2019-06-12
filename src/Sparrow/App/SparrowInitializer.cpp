/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SparrowInitializer.h"
#include <Core/Exceptions.h>
#include <Core/ModuleManager.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <boost/filesystem.hpp>

namespace Scine {
namespace Sparrow {

std::string SparrowInitializer::baseDirectory_;
std::string SparrowInitializer::libDirectory_;
std::string SparrowInitializer::resourceDirectory_;

SparrowInitializer::SparrowInitializer() : manager_(Core::ModuleManager::getInstance()) {
}

void SparrowInitializer::initialize() {
  setDirectories();

  auto libraryPath = boost::filesystem::path{libDirectory_} / "sparrow";

  if (!manager_.moduleLoaded("Sparrow")) {
    manager_.load(libraryPath);
  }
}

void SparrowInitializer::setDirectories() {
  auto pathToBinary = boost::dll::program_location();
  namespace fs = boost::filesystem;
  auto stripped_path_to_binary = fs::canonical(pathToBinary);

  auto pathToBinDir = stripped_path_to_binary.parent_path();

  auto pathToBuildDir = pathToBinDir.parent_path();

  // If the project is installed in Debug mode adapt the path accordingly
  auto installationRoot = pathToBuildDir;
  if (pathToBuildDir.filename() == "Debug") {
    installationRoot = pathToBuildDir.parent_path();
  }

  auto pathToResourceDir = installationRoot / "resources" / "Parameters";
  auto pathToLibDir = installationRoot / "lib";

  // If the project is not installed point to the right library location
  if (pathToBinDir.filename() != "bin") {
    pathToLibDir = pathToBinDir;
    pathToResourceDir = pathToBuildDir.parent_path().parent_path() / "src" / "Sparrow" / "Sparrow" / "Resources";
  }

  baseDirectory_ = installationRoot.string();
  libDirectory_ = pathToLibDir.string();
  resourceDirectory_ = pathToResourceDir.string();
}

std::string SparrowInitializer::getResourceDirectory() {
  if (resourceDirectory_.empty()) {
    throw(Core::InitializationException("Resource directory empty. The initializer did not run."));
  }
  return resourceDirectory_;
}

std::string SparrowInitializer::getBaseDirectory() {
  if (baseDirectory_.empty()) {
    throw(Core::InitializationException("Base directory empty. The initializer did not run."));
  }
  return baseDirectory_;
}

std::string SparrowInitializer::getLibDirectory() {
  if (baseDirectory_.empty()) {
    throw(Core::InitializationException("lib directory empty. The initializer did not run."));
  }
  return libDirectory_;
}

std::vector<std::string> SparrowInitializer::getAvailableMethods() const {
  return manager_.getLoadedModels("calculator");
}

Core::ModuleManager& SparrowInitializer::getManager() const {
  return manager_;
}

} // namespace Sparrow
} // namespace Scine