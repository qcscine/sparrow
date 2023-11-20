/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "DftbParameterSet.h"
#include "Common.h"
#include "Indent.h"
#include "Sparrow/Implementations/Dftb/ParameterSet.h"
#include <Utils/Geometry/ElementInfo.h>

namespace Scine {
namespace Sparrow {
namespace dftb {

/* Inline ostream writers */

std::ostream& operator<<(std::ostream& os, const SkfData::SameElementLine& line) {
  os << "SkfData::SameElementLine {"
     << join(line.Ed, line.Ep, line.Es, line.SPE, line.Ud, line.Up, line.Us, line.fd, line.fp, line.fs) << "}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const RepulsionParameters::Spline& spline) {
  os << "{" << join(spline.start, spline.end, spline.c0, spline.c1, spline.c2, spline.c3) << "}";
  return os;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine

namespace {

std::ostream& operator<<(std::ostream& os, const std::pair<int, int>& key) {
  // Take the safe route in case key ordering changes
  os << "std::make_pair(" << join(key.first, key.second) << ")";
  return os;
}

using Hubbard = Scine::Sparrow::dftb::SkfHubbardDerivatives;
using Spin = Scine::Sparrow::dftb::SkfSpinConstants;
using ParameterSet = Scine::Sparrow::dftb::ParameterSet;
using ParameterSetKey = Scine::Sparrow::dftb::ParameterSet::Key;
using PairData = Scine::Sparrow::dftb::SkfData;
using Repulsion = Scine::Sparrow::dftb::RepulsionParameters;

std::string setIdentifier(std::string a) {
  std::transform(std::begin(a), std::end(a), std::begin(a), [](const char x) {
    // Replace dashes with underscores
    if (x == '-') {
      return '_';
    }

    // Replace anything non-alphanumeric with x
    if (!std::isalnum(x)) {
      return 'x';
    }

    return x;
  });

  // Prepend param in all cases, avoids numeric character at front
  a = "params_" + a;
  return a;
}

std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr) {
  os << "{{" << arr[0] << sep << arr[1] << sep << arr[2] << "}}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const std::array<std::array<double, 3>, 3>& arr) {
  os << "{{" << arr[0] << sep << arr[1] << sep << arr[2] << "}}";
  return os;
}

/* Multiline ostream writers - Indent ostream specializations */
template<typename T>
std::ostream& operator<<(std::ostream& os, const Indented<T>& indentedType) {
  return os;
}

std::ostream& operator<<(std::ostream& os, const Indented<PairData::IntegralTable>& indent) {
  os << indent(0) << "{{" << nl;
  for (const auto& col : indent.bound) {
    if (std::all_of(std::begin(col), std::end(col), [](const double v) { return v == 0.0; })) {
      os << indent(1) << "std::vector<double>(" << col.size() << ", 0.0)" << sep << nl;
    }
    else {
      os << indent(1) << "{{";
      for (const auto& v : col) {
        os << v << sep;
      }
      os << "}}" << sep << nl;
    }
  }
  os << indent(0) << "}}" << sep << nl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const Indented<Repulsion>& indent) {
  const auto& rep = indent.bound;

  os << indent(0) << "{" << nl;
  os << indent(1) << join(rep.nSplineInts, rep.cutoff, rep.a1, rep.a2, rep.a3) << sep << nl;
  os << indent(1) << "{{" << nl;
  for (const Repulsion::Spline& spline : rep.splines) {
    os << indent(2) << spline << sep << nl;
  }
  os << indent(1) << "}}" << sep << nl;
  os << indent(1) << join(rep.c4, rep.c5) << nl;
  os << indent(0) << "}" << nl;
  return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const Indented<Spin>& indent) {
  const auto& spin = indent.bound;

  os << indent(0) << "return {{" << nl;
  for (const auto& mapPair : spin.map) {
    os << indent(1) << "{" << mapPair.first << sep << mapPair.second << "}," << nl;
  }
  os << indent(0) << "}};" << nl;
  return os;
}

template<>
std::ostream& operator<<(std::ostream& os, const Indented<Hubbard>& indent) {
  const auto& hubbard = indent.bound;

  os << indent(0) << "return {{" << nl;
  for (const auto& mapPair : hubbard.map) {
    os << indent(1) << "{" << join(mapPair.first, mapPair.second) << "}," << nl;
  }
  os << indent(0) << "}};" << nl;

  return os;
}

/* Identifier and generator fns */
std::string keyFnIdentifier(const ParameterSetKey& key, const std::string& setName) {
  const std::string e1 = Scine::Utils::ElementInfo::symbol(Scine::Utils::ElementInfo::element(key.first));
  const std::string e2 = Scine::Utils::ElementInfo::symbol(Scine::Utils::ElementInfo::element(key.second));
  return lower(setName + "_" + e1 + "_" + e2);
}

std::string keyFnSignature(const ParameterSetKey& key, const std::string& setName) {
  return "SkfData " + keyFnIdentifier(key, setName) + "()";
}

void writeParameterSetGenerator(std::ostream& os, const ParameterSet& parameterSet, const std::string& setName) {
  os << R"delim(namespace {

std::vector<bool> remapElements(const std::vector<int>& elements) {
  std::vector<bool> required(112, false);
  for(int Z : elements) {
    required.at(Z) = true;
  }
  return required;
}

} // namespace

)delim";
  os << "ParameterSet " << setName << "(const std::vector<int>& elements) {" << nl;
  os << Indent::level(1) << "ParameterSet data;" << nl << nl;
  os << Indent::level(1) << "if(elements.empty()) {" << nl;
  for (const auto& entry : parameterSet.pairData) {
    const auto& key = entry.first;
    os << Indent::level(2) << "data.pairData.emplace(" << key << ", " << keyFnIdentifier(key, setName) << "());" << nl;
  }
  os << Indent::level(1) << "} else {" << nl;
  os << Indent::level(2) << "const auto required = remapElements(elements);" << nl;
  for (const auto& entry : parameterSet.pairData) {
    const auto& key = entry.first;

    if (key.first == key.second) {
      os << Indent::level(2) << "if(required.at(" << key.first << ")) {" << nl;
    }
    else {
      os << Indent::level(2) << "if(required.at(" << key.first << ") && required.at(" << key.second << ")) {" << nl;
    }

    os << Indent::level(3) << "data.pairData.emplace(" << key << ", " << keyFnIdentifier(key, setName) << "());" << nl;
    os << Indent::level(2) << "}" << nl;
  }
  os << Indent::level(1) << "}" << nl;

  if (parameterSet.spin) {
    os << Indent::level(1) << "data.spin = " << setName + "_spin();" << nl;
  }
  if (parameterSet.hubbard) {
    os << Indent::level(1) << "data.hubbard = " << setName + "_hubbard();" << nl;
  }
  os << nl << Indent::level(1) << "return data;" << nl;
  os << "}" << nl << nl;
}

std::string spinFnIdentifier(const std::string& setName) {
  return setName + "_spin";
}

std::string spinFnSignature(const std::string& setName) {
  return "SkfSpinConstants " + spinFnIdentifier(setName) + "()";
}

std::string hubbardFnIdentifier(const std::string& setName) {
  return setName + "_hubbard";
}

std::string hubbardFnSignature(const std::string& setName) {
  return "SkfHubbardDerivatives " + hubbardFnIdentifier(setName) + "()";
}

std::string masterHeaderIncludePath(const boost::filesystem::path& headerPath) {
  boost::filesystem::path decomposition = headerPath;
  while (decomposition.filename() != "Sparrow" && decomposition != "") {
    decomposition = decomposition.parent_path();
  }
  return boost::filesystem::relative(headerPath, decomposition.parent_path()).string();
}

} // namespace

void dftbParameterSet(const boost::filesystem::path& directory) {
  using namespace Scine;
  using namespace Sparrow::dftb;
  namespace fs = boost::filesystem;

  std::cout << "Collecting all parameters at " << directory.string() << "...\n";
  const auto parameterSet = ParameterSet::collect(directory.string());

  std::string originalName;
  std::string setName;
  if (directory.filename_is_dot()) {
    originalName = directory.parent_path().stem().string();
    setName = setIdentifier(originalName);
  }
  else {
    originalName = directory.filename().string();
    setName = setIdentifier(originalName);
  }

  // Transform the directory name into a passable C++ identifier
  std::cout << "Writing '" << setName << "' generators ...\n";

  /* Master header: Declarations of
   *
   * - ParameterSet generator fn
   * - All individual atom-pair generators
   * - Spin and Hubbard generators (if present in the set)
   */
  CppFileStructure headerStructure{{"Sparrow/Implementations/Dftb/ParameterSet.h"},
                                   {"Scine", "Sparrow", "dftb"},
                                   "INCLUDE_SPARROW_IMPLEMENTATIONS_DFTB_" + upper(setName) + "_PARAMETER_SET_H"};

  const fs::path headerPath = directory / "Parameters.h";
  const fs::path implPath = directory / "Parameters.cpp";
  std::ofstream masterHeader(headerPath.string());
  headerStructure.writeHeader(masterHeader);
  masterHeader << R"delim(/*! @brief Dftb )delim" << originalName << R"delim( parameter set generator
 *
 * @param elements Elements for which to load parameters. If empty, loads all
 *   available parameters.
 */
)delim";
  masterHeader << "ParameterSet " + setName + "(const std::vector<int>& elements = {});" << nl << nl;

  masterHeader << R"delim(// NOTE: Atom-pair data generator fns)delim" << nl;
  for (const auto& entry : parameterSet.pairData) {
    masterHeader << keyFnSignature(entry.first, setName) << ";" << nl;
  }

  masterHeader << nl;
  if (parameterSet.spin) {
    masterHeader << R"delim(// Spin constants)delim" << nl;
    masterHeader << spinFnSignature(setName) << ";" << nl << nl;
  }

  if (parameterSet.hubbard) {
    masterHeader << R"delim(// Hubbard derivatives)delim" << nl;
    masterHeader << hubbardFnSignature(setName) << ";" << nl;
  }

  headerStructure.writeFooter(masterHeader);

  /* Master impl: Definitions of
   *
   * - ParameterSet generator fn
   * - Spin and Hubbard generators (if present in the set)
   *
   * and then separate cpps for each parameter pair fn
   */
  CppFileStructure implStructure{{"Sparrow/Implementations/Dftb/ParameterSet.h", masterHeaderIncludePath(headerPath)},
                                 {"Scine", "Sparrow", "dftb"}};
  std::ofstream masterImpl(implPath.string());
  implStructure.writeHeader(masterImpl);
  masterImpl << std::hexfloat;
  writeParameterSetGenerator(masterImpl, parameterSet, setName);
  if (parameterSet.spin) {
    masterImpl << spinFnSignature(setName) << " {" << nl;
    masterImpl << Indent::level(1) << R"delim(// clang-format off)delim" << nl;
    masterImpl << Indent::bind(1, parameterSet.spin.value());
    masterImpl << Indent::level(1) << R"delim(// clang-format on)delim" << nl;
    masterImpl << "}" << nl;
  }
  if (parameterSet.hubbard) {
    masterImpl << hubbardFnSignature(setName) << " {" << nl;
    masterImpl << Indent::level(1) << R"delim(// clang-format off)delim" << nl;
    masterImpl << Indent::bind(1, parameterSet.hubbard.value());
    masterImpl << Indent::level(1) << R"delim(// clang-format on)delim" << nl;
    masterImpl << "}" << nl;
  }
  /* All atom pair generator impls */
  for (const auto& entryPair : parameterSet.pairData) {
    const auto& key = entryPair.first;
    const auto& entry = entryPair.second;
    masterImpl << std::hexfloat;
    masterImpl << keyFnSignature(key, setName) << " {" << nl;
    masterImpl << Indent::level(1) << R"delim(// clang-format off)delim" << nl;
    masterImpl << Indent::level(1) << "return {" << nl;
    masterImpl << Indent::level(2) << entry.gridDistance << sep << nl;
    masterImpl << Indent::level(2) << wrapOptional(entry.atomicParameters) << sep << nl;
    masterImpl << Indent::bind(2, entry.integralTable);
    masterImpl << Indent::bind(2, entry.repulsion);
    masterImpl << Indent::level(1) << "};" << nl;
    masterImpl << Indent::level(1) << R"delim(// clang-format on)delim" << nl;
    masterImpl << "}" << nl;
  }
  implStructure.writeFooter(masterImpl);
}
