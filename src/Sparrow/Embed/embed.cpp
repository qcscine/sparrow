/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Basisfile.h"
#include "DftbParameterSet.h"
#include "NddoParameters.h"
#include <iostream>

const std::string usage = R"delim(
This binary is a serializer/meta-compiler of sorts.
It generates C++ .h and .cpp files that contain functions generating
instances for three runtime types:

- Sparrow::dftb::ParameterSet (a directory of SKF files and optionally a
  spin.dat and a hubbard.dat)
  - See dftbParameterSet
  - Generates a header and impl that construct only the parameters that
    are needed from a passed list of atomic numbers in a passed directory
- Sparrow::nddo::Parameters (runtime type directly mappable to custom JSON
  archive structure)
  - See nddoParameters
  - Generates a single parameters.cpp file in the same location as passed archive
- int -> Utils::GtoExpansion map (basically the result of parsing a
  turbomole minimal basisfile)
  - See basisfile
  - Generates a single cpp in the same location as passed file

Usage: Pass as arguments:

- directory: DFTB SKF parameters
- *.basis: Turbomole basisfile
- *.json: In-house NDDO parameter archive format
)delim";

int main(int argc, char* argv[]) {
  namespace fs = boost::filesystem;

  bool errors = false;

  for (int i = 1; i < argc; ++i) {
    const std::string argstr{argv[i]};

    if (argstr == "-h" || argstr == "--help") {
      std::cout << usage;
      return 0;
    }

    const fs::path arg{argstr};
    if (!fs::exists(arg)) {
      std::cout << "Argument '" << arg.string() << "' does not exist.\n";
      errors = true;
      continue;
    }

    if (fs::is_directory(arg)) {
      std::cout << "Passed a directory, assuming that it contains DFTB SKF parameters\n";
      dftbParameterSet(arg);
    }

    if (fs::is_regular_file(arg)) {
      if (arg.extension() == ".basis") {
        std::cout << "Passed a basisfile. Assuming it's turbomole format\n";
        basisfile(arg);
      }
      else if (arg.extension() == ".json") {
        std::cout << "Passed a JSON file. Assuming it's our nddo parameter archive format\n";
        nddoParameters(arg);
      }
    }
  }

  if (errors) {
    return 1;
  }

  return 0;
}
