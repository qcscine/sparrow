/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_SPARROW_IMPLEMENTATIONS_DFTB_PARAMETER_SET_H
#define INCLUDE_SPARROW_IMPLEMENTATIONS_DFTB_PARAMETER_SET_H

#include "Sparrow/Implementations/Dftb/Utils/SkfParser.h"
#include "boost/functional/hash.hpp"

namespace Scine {
namespace Sparrow {
namespace dftb {

struct ParameterSet {
  using Key = std::pair<int, int>;
  using DiatomicParameters = std::unordered_map<Key, SkfData, boost::hash<Key>>;

  /*! @brief Collects parameters from a directory
   *
   * @param directory Directory to scour for parameters
   * @param elements Unique atomic numbers of elements for which to collect
   *   parameters
   *
   * Non-recursively traverses @p directory. If @p elements is empty, reads all
   * files with extension .skf. If @p elements is non-empty, reads only .skf
   * files of the required elements, trying names with and without a dash
   * separating element names.
   *
   * The "hubbard.dat" and "spin.dat" files are sought unconditionally.
   *
   * @throws std::runtime_error If an skf file cannot be parsed, a required
   * element pair's skf file cannot be found, or during directory traversal if
   * an skf file is encountered whose name canot be parsed into element names.
   */
  static ParameterSet collect(const std::string& directory, const std::vector<int>& elements = {});

  //! Overrides stored parameters with those from another parameter set
  ParameterSet& patch(ParameterSet patchSet);

  DiatomicParameters pairData;
  boost::optional<SkfSpinConstants> spin;
  boost::optional<SkfHubbardDerivatives> hubbard;
};

} // namespace dftb
} // namespace Sparrow
} // namespace Scine

#endif
