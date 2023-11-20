/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "NddoParameters.h"
#include "Common.h"
#include "Sparrow/Implementations/Nddo/Parameters.h"

namespace {

using Parameters = Scine::Sparrow::nddo::Parameters;

std::ostream& operator<<(std::ostream& os, const Parameters::Atomic::Pack::Spd& spd) {
  if (spd.s != 0.0 || spd.p != 0.0 || spd.d != 0.0) {
    os << "Spd {" << spd.s << sep << spd.p << sep << spd.d << "}";
  }
  else {
    os << "Spd {}";
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const Parameters::Atomic::Pack& pack) {
  os << "Pack {";
  os << pack.oneCenterEnergy << sep;
  os << pack.beta << sep;
  os << pack.orbitalExponent << sep;
  os << pack.internalExponent << sep;
  const auto end = std::next(std::addressof(pack.alpha));
  for (auto it = std::addressof(pack.gss); it != end; ++it) {
    os << *it;
    if (it != end - 1) {
      os << sep;
    }
  }
  os << "}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Parameters::Atomic::GaussianRepulsion& gauss) {
  if (gauss.a != 0.0 || gauss.b != 0.0 || gauss.c != 0.0) {
    os << "Gauss {" << gauss.a << sep << gauss.b << sep << gauss.c << "}";
  }
  else {
    os << "Gauss {}";
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const std::pair<int, Parameters::Atomic>& intAtomicPair) {
  os << "{" << intAtomicPair.first << sep << "{" << intAtomicPair.second.pack;
  os << sep << "{{";
  const auto end = std::end(intAtomicPair.second.gaussianRepulsion);
  for (auto it = std::begin(intAtomicPair.second.gaussianRepulsion); it != end; ++it) {
    os << *it;

    if (it != end - 1) {
      os << sep;
    }
  }
  os << "}}";
  os << "}}";
  return os;
}

std::ostream& operator<<(std::ostream& os, const std::pair<Parameters::DiatomicKey, Parameters::Diatomic>& diatomicPair) {
  os << "{{" << diatomicPair.first.first << sep << diatomicPair.first.second << "}" << sep << "{"
     << diatomicPair.second.exponent << sep << diatomicPair.second.factor << "}}";
  return os;
}

void write(std::ostream& os, const Parameters& parameters) {
  CppFileStructure structure{{{"Sparrow/Implementations/Nddo/Parameters.h"}}, {{"Scine", "Sparrow", "nddo"}}};
  structure.writeHeader(os);
  os << "Parameters method() {" << nl;
  os << sp(2) << R"delim(// clang-format off)delim" << nl;
  os << sp(2) << "using Atomic = Parameters::Atomic;" << nl;
  os << sp(2) << "using Pack = Parameters::Atomic::Pack;" << nl;
  os << sp(2) << "using Spd = Parameters::Atomic::Pack::Spd;" << nl;
  os << sp(2) << "using Gauss = Parameters::Atomic::GaussianRepulsion;" << nl << nl;
  os << sp(2) << "using Diatomic = Parameters::Diatomic;" << nl;
  os << std::hexfloat;
  os << sp(2) << "return {" << nl;
  os << sp(4) << "{" << nl;
  {
    const unsigned size = parameters.atomic.size();
    unsigned count = 0;
    for (const auto& intAtomicPair : parameters.atomic) {
      os << sp(6) << intAtomicPair;

      if (count < size - 1) {
        os << sep;
      }

      os << nl;
      ++count;
    }
  }
  os << sp(4) << "}";
  if (!parameters.diatomic.empty()) {
    os << "," << nl;
    os << sp(4) << "{" << nl;
    {
      const unsigned size = parameters.diatomic.size();
      unsigned count = 0;
      for (const auto& diatomicPair : parameters.diatomic) {
        os << sp(6) << diatomicPair;

        if (count < size - 1) {
          os << sep;
        }

        os << nl;
        ++count;
      }
    }
    os << sp(4) << "}" << nl;
  }
  else {
    os << nl;
  }
  os << sp(2) << "};" << nl;
  os << sp(2) << R"delim(// clang-format on)delim" << nl;
  os << "}" << nl << nl;
  structure.writeFooter(os);
}

} // namespace

//! Writes a cpp file containing the parameters
void nddoParameters(boost::filesystem::path path) {
  auto parameters = Parameters::read(path.string());
  std::ofstream outfile(path.replace_extension(".cpp").string());
  write(outfile, parameters);
}
