/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Basisfile.h"
#include "Common.h"
#include <Utils/IO/TurbomoleMinimalBasisfile.h>

namespace Scine {
namespace Utils {

std::ostream& operator<<(std::ostream& os, const GtoExpansion& expansion) {
  os << "Utils::GtoExpansion {" << expansion.angularMomentum << sep << "{{";
  for (auto iter = std::begin(expansion.gtfs); iter != std::end(expansion.gtfs); ++iter) {
    os << "{" << join(expansion.angularMomentum, iter->exponent, iter->coefficient) << "}";

    if (iter != std::end(expansion.gtfs) - 1) {
      os << sep;
    }
  }
  os << "}}}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const AtomicGtos& atomics) {
  os << "{" << join(wrapOptional(atomics.s), wrapOptional(atomics.p), wrapOptional(atomics.d)) << "}";
  return os;
}

} // namespace Utils
} // namespace Scine

void basisfile(boost::filesystem::path source) {
  std::cout << "Writing CPP for turbomole basisfile " << source.string() << "...\n";
  auto parameters = Scine::Utils::readTurbomoleBasisfile(source.string());

  CppFileStructure structure{{"Utils/DataStructures/AtomicGtos.h", "<unordered_map>"}, {"Scine", "Sparrow", "Sto6g"}};

  std::ofstream outfile(source.replace_extension("cpp").string());
  structure.writeHeader(outfile);

  outfile << std::hexfloat;
  outfile << "std::unordered_map<int, Utils::AtomicGtos> method() {" << nl;
  outfile << sp(2) << R"delim(// clang-format off)delim" << nl;
  outfile << sp(2) << "return {" << nl;
  for (const auto& value : parameters) {
    outfile << sp(4) << "{" << join(value.first, value.second) << "}," << nl;
  }
  outfile << sp(2) << "};" << nl;
  outfile << sp(2) << R"delim(// clang-format on)delim" << nl;
  outfile << "}" << nl;

  structure.writeFooter(outfile);
}
