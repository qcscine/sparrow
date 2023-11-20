/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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
  /**
   * @brief Initializes the base (root of the installation/build),
   *        the resource and the lib directory (for the installation).
   * Initializes the relevant paths and loads the Sparrow library. Only works for the installed project.
   */
  static void initialize();
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_SPARROWINITIALIZER_H
