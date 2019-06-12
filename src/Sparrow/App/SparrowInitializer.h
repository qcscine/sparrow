/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_SPARROWINITIALIZER_H
#define SPARROW_SPARROWINITIALIZER_H

#include <string>
#include <vector>

namespace Scine {

namespace Core {
class ModuleManager;
}

namespace Sparrow {

/**
 * @brief This class initializes the resource and base directories from the position of the program.
 *
 * This initialization is particularly important for finding the default parameters in the resource folder
 * in the installation path.
 */
class SparrowInitializer {
 public:
  SparrowInitializer();
  /**
   * @brief Initializes the base (root of the installation/build),
   *        the resource and the lib directory (for the installation).
   * Initializes the relevant paths and loads the Sparrow library. Only works for the installed project.
   */
  void initialize();

  /** @brief getter for the resource directory */
  static std::string getResourceDirectory();
  /** @brief getter fot the base directory */
  static std::string getBaseDirectory();
  /** @brief getter fot the base directory */
  static std::string getLibDirectory();

  /** @brief getter for the list of available methods. */
  std::vector<std::string> getAvailableMethods() const;

  Core::ModuleManager& getManager() const;

 private:
  void setDirectories();
  static std::string baseDirectory_;
  static std::string libDirectory_;
  static std::string resourceDirectory_;
  Core::ModuleManager& manager_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_SPARROWINITIALIZER_H
