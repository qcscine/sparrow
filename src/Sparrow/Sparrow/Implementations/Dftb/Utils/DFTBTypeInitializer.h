/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTBTYPEINITIALIZER_H
#define SPARROW_DFTBTYPEINITIALIZER_H

#include <Utils/MethodEssentials/MethodFactories/MethodInitializer.h>

namespace Scine {
namespace Sparrow {

namespace dftb {

/*!
 * Common interface for the initializer of DFTB type methods.
 */

class DFTBTypeInitializer : public Utils::MethodInitializer {
 public:
  explicit DFTBTypeInitializer(std::string set);

  ///*! Set the folder containing the DFTB parameter sets. */
  // void setParameterFolder(std::string path);
  /*! Get the folder name for the DFTB parameter sets. */
  std::string getParameterFolder() const;
  /*! Set the name of the parameter set to use. */
  void setParameterSet(std::string set);
  /*! Set the name of the parameter set to use. */
  std::string getParameterSet() const;

 protected:
  std::string parameterSet_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_DFTBTYPEINITIALIZER_H