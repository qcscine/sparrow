/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Sparrow/Implementations/Dftb/ParameterSet.h"
#include "Utils/Geometry/ElementInfo.h"
#include "boost/filesystem.hpp"
#include "boost/fusion/adapted/std_pair.hpp"
#include "boost/spirit/include/qi.hpp"

namespace Scine {
namespace Sparrow {
namespace dftb {
namespace {

// Mapping from element type name strings to their atomic number
struct element_type_ : boost::spirit::qi::symbols<char, int> {
  element_type_() {
    for (const auto& nameElementPair : Utils::ElementInfo::stringToElementType()) {
      /* Since we will be using these in case-insensitive matching, we need to
       * make the names lowercase
       */
      std::string lowercaseName = nameElementPair.first;
      std::transform(std::begin(lowercaseName), std::end(lowercaseName), std::begin(lowercaseName),
                     [](unsigned char c) { return std::tolower(c); });
      if (lowercaseName != "none") {
        add(lowercaseName, Utils::ElementInfo::Z(nameElementPair.second));
      }
    }
  };
} element_type;

boost::filesystem::path trySkfNames(const boost::filesystem::path& directory, const std::string& a, const std::string& b) {
  namespace fs = boost::filesystem;
  auto tryPattern = [&directory](std::string name) -> boost::optional<fs::path> {
    fs::path fullPath = directory / name;
    if (fs::exists(fullPath) && fs::is_regular_file(fullPath)) {
      return fullPath;
    }

    // Try an all-lowercase version too
    std::transform(std::begin(name), std::end(name), std::begin(name), [](unsigned char c) { return std::tolower(c); });

    fullPath = directory / name;
    if (fs::exists(fullPath) && fs::is_regular_file(fullPath)) {
      return fullPath;
    }

    return boost::none;
  };

  auto maybePath = tryPattern(a + "-" + b + ".skf");

  if (!maybePath) {
    maybePath = tryPattern(a + b + ".skf");
  }

  if (!maybePath) {
    throw std::runtime_error("Could not find a required element SKF file");
  }

  return maybePath.value();
}

} // namespace

ParameterSet ParameterSet::collect(const std::string& directory, const std::vector<int>& elements) {
  namespace fs = boost::filesystem;
  namespace qi = boost::spirit::qi;
  auto directoryPath = fs::path(directory);

  if (!fs::is_directory(directory)) {
    throw std::runtime_error("Argument to parameter set collection is not a directory");
  }

  ParameterSet parameters;

  if (elements.empty()) {
    // Read all parameter files we can find in the directory
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
      if (fs::is_regular_file(entry)) {
        const std::string extension = entry.path().extension().string();
        if (extension != ".skf") {
          continue;
        }

        const std::string stem = entry.path().stem().string();

        // NOTE: This parses BBr and BrB stems correctly
        std::pair<int, int> zPair;
        const bool successfulParse =
            qi::parse(std::begin(stem), std::end(stem), qi::no_case[element_type >> -qi::lit("-") >> element_type], zPair);

        if (!successfulParse) {
          throw std::runtime_error("Could not parse '" + stem + "' into element types!");
        }

        parameters.pairData.emplace(zPair, SkfData::read(entry.path().string()));
      }
    }
  }
  else {
    // Read only the required element parameter files
    const unsigned E = elements.size();
    std::vector<std::string> elementNames;
    elementNames.reserve(E);
    std::transform(
        std::begin(elements), std::end(elements), std::back_inserter(elementNames),
        [](const unsigned Z) -> std::string { return Utils::ElementInfo::symbol(Utils::ElementInfo::element(Z)); });

    for (unsigned i = 0; i < E; ++i) {
      for (unsigned j = 0; j < E; ++j) {
        const fs::path skf = trySkfNames(directory, elementNames.at(i), elementNames.at(j));
        parameters.pairData.emplace(std::make_pair(elements.at(i), elements.at(j)), SkfData::read(skf.string()));
      }
    }
  }

  const auto hubbardPath = directoryPath / "hubbard.dat";
  if (fs::is_regular_file(hubbardPath)) {
    parameters.hubbard = SkfHubbardDerivatives::read(hubbardPath.string());
  }

  const auto spinPath = directoryPath / "spin.dat";
  if (fs::is_regular_file(spinPath)) {
    parameters.spin = SkfSpinConstants::read(spinPath.string());
  }

  return parameters;
}

ParameterSet& ParameterSet::patch(ParameterSet patchSet) {
  // Patch diatomics
  for (auto& diatomicPair : patchSet.pairData) {
    auto findIter = pairData.find(diatomicPair.first);
    if (findIter == std::end(pairData)) {
      pairData.emplace_hint(findIter, diatomicPair);
    }
    else {
      findIter->second = std::move(diatomicPair.second);
    }
  }

  // Patch spins
  if (spin && patchSet.spin) {
    spin.value().patch(patchSet.spin.value());
  }
  else if (!spin && patchSet.spin) {
    spin = patchSet.spin;
  }
  // NOTE: All other cases require no work

  // Patch hubbard
  if (hubbard && patchSet.hubbard) {
    hubbard.value().patch(patchSet.hubbard.value());
  }
  else if (!hubbard && patchSet.hubbard) {
    hubbard = patchSet.hubbard;
  }
  // NOTE: All other cases require no work

  return *this;
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
