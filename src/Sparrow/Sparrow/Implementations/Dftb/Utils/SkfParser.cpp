/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Sparrow/Implementations/Dftb/Utils/SkfParser.h"
#include "Sparrow/Implementations/Dftb/Utils/RepulsionParameters.h"
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Scf/MethodExceptions.h"
#include "boost/filesystem.hpp"
#include "boost/fusion/adapted/struct/adapt_struct.hpp"
#include "boost/fusion/include/adapt_struct.hpp"
#include "boost/phoenix/fusion/at.hpp"
#include "boost/phoenix/object/static_cast.hpp"
#include "boost/phoenix/stl/container.hpp"
#include "boost/spirit/include/phoenix.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/qi.hpp"
#include <fstream>
#include <iostream>

BOOST_FUSION_ADAPT_STRUCT(Scine::Sparrow::dftb::RepulsionParameters::Spline, (double, start), (double, end),
                          (double, c0), (double, c1), (double, c2), (double, c3))

BOOST_FUSION_ADAPT_STRUCT(Scine::Sparrow::dftb::RepulsionParameters, (int, nSplineInts), // 0
                          (double, cutoff), (double, a1),                                // 2
                          (double, a2), (double, a3),
                          (std::vector<Scine::Sparrow::dftb::RepulsionParameters::Spline>, splines), // 5
                          (double, c4),                                                              // 6
                          (double, c5))

BOOST_FUSION_ADAPT_STRUCT(Scine::Sparrow::dftb::SkfData::SameElementLine, (double, Ed), (double, Ep), (double, Es),
                          (double, SPE), (double, Ud), (double, Up), (double, Us), (unsigned, fd), (unsigned, fp),
                          (unsigned, fs))

namespace Scine {
namespace Sparrow {
namespace dftb {

namespace {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template<typename Iterator>
struct Skf : qi::grammar<Iterator> {
  using DoublesList = std::vector<double>;
  using NestedDoubles = std::vector<DoublesList>;

  std::string forwardName;
  double gridDistance;
  boost::optional<SkfData::SameElementLine> atomicParameters;
  NestedDoubles integralTable;
  RepulsionParameters repulsion;

  void clearState() {
    forwardName.clear();
    atomicParameters = boost::none;
    integralTable.clear();
    repulsion.splines.clear();
  }

  Skf() : Skf::base_type(start) {
    // clang-format off
    namespace bp = boost::phoenix;
    using bp::at_c;
    using bp::push_back;
    using bp::reserve;
    using bp::resize;
    using bp::size;
    using qi::_1;
    using qi::_val;
    using qi::char_;
    using qi::double_;
    using qi::eol;
    using qi::lit;
    using qi::uint_;

    splineToken %= qi::skip(ascii::blank | qi::lit(','))[
      double_ >> double_ >> double_ >> double_ >> double_ >> double_
    ];
    splineToken.name("spline token");

    /* Spline blocks have a pretty fixed format:
     * - header
     * - the number of spline integers and a cutoff
     * - three parameters
     * - lines with start, end, c0 - c3 (a spline)
     * - last line with start, end, c0 - c5 (spline + two parameters)
     */
    splineBlock = qi::skip(ascii::blank | qi::lit(','))[
      qi::lit("Spline") > +eol
      > uint_[at_c<0>(_val) = _1, reserve(at_c<5>(_val), _1)]
      > double_[at_c<1>(_val) = _1]
      > eol
      > double_[at_c<2>(_val) = _1]
      > double_[at_c<3>(_val) = _1]
      > double_[at_c<4>(_val) = _1]
      > eol
      > +(
        splineToken[push_back(at_c<5>(_val), _1)]
        >> -(
          double_[at_c<6>(_val) = _1]
          > double_[at_c<7>(_val) = _1]
        )
        > eol
      )
    ];
    splineBlock.name("spline block");
    qi::on_error<qi::fail>(splineBlock, std::cerr << bp::val("Expected spline block here:")
                                                  << bp::construct<std::string>(qi::_3, qi::_2) << bp::val("\"\n"));

    /* Integral lines can have either double tokens or something like this:
     * 8*1.4, which means 8 entries with value 1.4. We don't really know how many
     * tokens we're going to get on a line, so expect at least one.
     */
    integralLine = (
      qi::eps[reserve(_val, 20)]
      >> qi::skip(ascii::blank | qi::lit(','))[
        +(
          (uint_ >> qi::lit('*') >> double_)[resize(_val, size(_val) + _1, qi::_2)]
          | double_[push_back(_val, _1)]
        )
        > eol
      ]
    );
    integralLine.name("integral line");
    qi::on_error<qi::fail>(integralLine, std::cerr << bp::val("Expected integral line here:")
                                                   << bp::construct<std::string>(qi::_3, qi::_2) << bp::val("\"\n"));

    /* There can be up to three header lines, one optional with atomic
     * parameters if the skf file is for matching atoms.
     */
    headers = qi::skip(ascii::blank | qi::lit(','))[
      double_[bp::ref(gridDistance) = _1] > uint_ >> -uint_ > eol
        >> -sameElementLine[bp::ref(atomicParameters) = _1]
        // One line of either 20 doubles or fewer int*double tokens which we ignore
        >> +(qi::lexeme[uint_ >> qi::char_('*') >> double_] | double_) > eol
    ];

    headers.name("headers");
    qi::on_error<qi::fail>(headers, std::cerr << bp::val("Expected skf headers here:")
                                              << bp::construct<std::string>(qi::_3, qi::_2) << bp::val("\"\n"));

    // There can be a line with ten atomic parameters if the atoms match
    sameElementLine = qi::skip(ascii::blank | qi::lit(','))[
      double_[at_c<0>(_val) = _1]
      >> double_[at_c<1>(_val) = _1]
      >> double_[at_c<2>(_val) = _1]
      >> double_[at_c<3>(_val) = _1]
      >> double_[at_c<4>(_val) = _1]
      >> double_[at_c<5>(_val) = _1]
      >> double_[at_c<6>(_val) = _1]
      >> double_[at_c<7>(_val) = bp::static_cast_<unsigned>(_1)]
      >> double_[at_c<8>(_val) = bp::static_cast_<unsigned>(_1)]
      >> double_[at_c<9>(_val) = bp::static_cast_<unsigned>(_1)]
      >> eol
    ];
    sameElementLine.name("element parameters");

    // Instead of any data, the file can just consist of a forwarding filename
    skfFilename = +qi::char_(R"(A-Za-z\-)") >> qi::lit(".skf");
    skfFilename.name("skf filename");

    start = (
      skfFilename[bp::ref(forwardName) = _1]
      | (
        headers
        > +(integralLine[push_back(bp::ref(integralTable), _1)])
        > qi::skip(ascii::space)[*eol] > splineBlock[bp::ref(repulsion) = _1]
      )
    );
    start.name("start");
    qi::on_error<qi::fail>(start, std::cerr << bp::val("Expected basic SKF structure here:")
                                            << bp::construct<std::string>(qi::_3, qi::_2) << bp::val("\"\n"));

    // NOTE: If you need to debug this grammar, uncomment below:
    // qi::debug(splineToken);
    // qi::debug(splineBlock);
    // qi::debug(sameElementLine);
    // qi::debug(integralLine);
    // qi::debug(headers);
    // qi::debug(start);

    // clang-format on
  }

  qi::rule<Iterator, dftb::RepulsionParameters::Spline()> splineToken;
  qi::rule<Iterator, dftb::RepulsionParameters()> splineBlock;
  qi::rule<Iterator, std::vector<double>()> integralLine;
  qi::rule<Iterator, SkfData::SameElementLine()> sameElementLine;
  qi::rule<Iterator> headers;
  qi::rule<Iterator, std::string()> skfFilename;
  qi::rule<Iterator> start;
};

} // namespace

SkfData SkfData::read(const std::string& filename) {
  namespace fs = boost::filesystem;
  if (!fs::exists(filename)) {
    throw std::runtime_error("File to read does not exist");
  }

  if (!fs::is_regular_file(filename)) {
    throw std::runtime_error("Argument to SKF parse is not a file");
  }

  std::ifstream ifs(filename);
  ifs.unsetf(std::ios::skipws);

  using Iterator = boost::spirit::istream_iterator;
  using Parser = Skf<Iterator>;

  bool parsingSuccess;
  bool forwarding;
  Parser parser;
  do {
    Iterator iter(ifs);
    Iterator end;
    parsingSuccess = boost::spirit::qi::parse(iter, end, parser);
    forwarding = !parser.forwardName.empty();

    if (parsingSuccess && forwarding) {
      auto skfPath = boost::filesystem::path(filename).parent_path() / parser.forwardName;
      skfPath += ".skf";
      if (!fs::exists(skfPath) || !fs::is_regular_file(skfPath)) {
        throw std::runtime_error("SKF forwarding could not be resolved!");
      }

      ifs.close();
      ifs.open(skfPath.string());
      parser.clearState();
    }
  } while (parsingSuccess && forwarding);

  if (!parsingSuccess) {
    throw std::runtime_error("Failed to parse SKF file");
  }

  // Convert into SkfData
  SkfData data;
  data.gridDistance = parser.gridDistance;
  data.atomicParameters = std::move(parser.atomicParameters);
  data.repulsion = std::move(parser.repulsion);

  /* Transpose the integral table from Nx20 to 28xN (where the last 8 columns
   * are default-initialized).
   */
  const unsigned N = parser.integralTable.size();
  for (auto& column : data.integralTable) {
    column.resize(N);
  }

  for (unsigned i = 0; i < N; ++i) {
    for (unsigned j = 0; j < 20; ++j) {
      data.integralTable.at(j).at(i) = parser.integralTable.at(i).at(j);
    }
  }

  return data;
}

SkfSpinConstants SkfSpinConstants::read(const std::string& filename) {
  SkfSpinConstants result;

  std::ifstream fSpin(filename);
  fSpin.imbue(std::locale("C"));
  // Output warning if it doesn't exist
  if (!fSpin) {
    throw Utils::Methods::ParameterFileCannotBeOpenedException(filename);
  }

  std::string atom;
  std::stringstream buffer;
  while (getline(fSpin, atom)) {
    // skip comment lines denoted by #
    if (atom.find('#') == std::string::npos) {
      buffer = std::stringstream(atom);
      buffer >> atom;

      const int Z = Utils::ElementInfo::Z(Utils::ElementInfo::elementTypeForSymbol(atom));
      SkfSpinConstants::MatrixType sc;
      buffer >> sc[0][0] >> sc[0][1] >> sc[1][0] >> sc[1][1] >> sc[0][2] >> sc[1][2] >> sc[2][2] >> sc[2][0] >> sc[2][1];
      result.map.emplace(Z, std::move(sc));
    }
  }

  return result;
}

void SkfSpinConstants::patch(SkfSpinConstants other) {
  for (auto& entry : other.map) {
    auto findIter = map.find(entry.first);
    if (findIter == std::end(map)) {
      map.emplace_hint(findIter, entry);
    }
    else {
      findIter->second = std::move(entry.second);
    }
  }
}

SkfHubbardDerivatives SkfHubbardDerivatives::read(const std::string& filename) {
  SkfHubbardDerivatives result;

  std::ifstream fHubbard(filename);
  fHubbard.imbue(std::locale("C"));
  if (!fHubbard) {
    throw Utils::Methods::ParameterFileCannotBeOpenedException(filename);
  }

  std::string atom;
  fHubbard >> atom;
  while (!fHubbard.eof()) {
    double hubbard;
    fHubbard >> hubbard;

    const int Z = Utils::ElementInfo::Z(Utils::ElementInfo::elementTypeForSymbol(atom));
    result.map.emplace(Z, hubbard);
    fHubbard >> atom;
  }

  return result;
}

void SkfHubbardDerivatives::patch(SkfHubbardDerivatives other) {
  for (const auto& entry : other.map) {
    auto findIter = map.find(entry.first);
    if (findIter == std::end(map)) {
      map.emplace_hint(findIter, entry);
    }
    else {
      findIter->second = entry.second;
    }
  }
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
