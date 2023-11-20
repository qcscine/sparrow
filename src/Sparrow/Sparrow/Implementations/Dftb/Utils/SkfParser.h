/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef INCLUDE_SPARROW_IMPLEMENTATIONS_DFTB_UTILS_SKF_PARSER_H
#define INCLUDE_SPARROW_IMPLEMENTATIONS_DFTB_UTILS_SKF_PARSER_H

#include "Sparrow/Implementations/Dftb/Utils/RepulsionParameters.h"
#include "boost/optional.hpp"
#include <array>
#include <string>
#include <unordered_map>

namespace Scine {
namespace Sparrow {
namespace dftb {

struct SkfData {
  struct SameElementLine {
    // On-site energies for angular momenta d, p, s
    double Ed;
    double Ep;
    double Es;
    // spin polarisation error for calculating formation energies
    double SPE;
    // Hubbard U values for appropriate angular momenta
    double Ud;
    double Up;
    double Us;
    // Occupations for the neutral atom
    unsigned fd;
    unsigned fp;
    unsigned fs;
  };

  using DoublesList = std::vector<double>;
  using IntegralTable = std::array<DoublesList, 28>;

  /*! @brief Parses a SKF file
   *
   * Resolves forwarding.
   */
  static SkfData read(const std::string& filename);

  double gridDistance;
  boost::optional<SameElementLine> atomicParameters;
  IntegralTable integralTable;
  RepulsionParameters repulsion;
};

struct SkfSpinConstants {
  using MatrixType = std::array<std::array<double, 3>, 3>;
  using MapType = std::unordered_map<int, MatrixType>;
  MapType map;

  static SkfSpinConstants read(const std::string& filename);

  void patch(SkfSpinConstants other);
};

struct SkfHubbardDerivatives {
  using MapType = std::unordered_map<int, double>;
  MapType map;

  static SkfHubbardDerivatives read(const std::string& filename);

  void patch(SkfHubbardDerivatives other);
};

} // namespace dftb
} // namespace Sparrow
} // namespace Scine

#endif
