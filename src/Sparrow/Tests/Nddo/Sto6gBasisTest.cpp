/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/Nddo/Parameters.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterTwoElectronIntegrals.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/AtomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/PM6DiatomicParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParameterProcessor.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/Gtf.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Typenames.h>
#include <gmock/gmock.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace Scine {
namespace Sparrow {

using namespace testing;
using namespace nddo;

class Sto6gBasisTest : public Test {
 public:
  bool coefficientsAreNonZero(const std::vector<Utils::Gtf>& gtfs) {
    return !std::all_of(gtfs.begin(), gtfs.end(), [](Utils::Gtf gtf) { return gtf.normalizedCoefficient == 0.0; });
  };

  void writeGtfLine(std::string& out, const Utils::Gtf& gtf) {
    double exp = gtf.exponent;
    double coeff = gtf.normalizedCoefficient;
    std::ostringstream streamObj;
    streamObj.precision(15);
    streamObj << std::scientific;
    streamObj << "    " << exp;
    if (coeff < 0) {
      streamObj << "   " << coeff << "\n";
    }
    else {
      streamObj << "    " << coeff << "\n";
    }
    out += streamObj.str();
  };

  std::string convertGtoToString(const Utils::ElementType& ele, const Utils::AtomicGtos& aos, std::string methodName) {
    std::string out = "*\n";
    std::string symbol = Utils::ElementInfo::symbol(ele);
    std::for_each(symbol.begin(), symbol.end(), [](char& c) { c = ::tolower(c); });
    std::for_each(methodName.begin(), methodName.end(), [](char& c) { c = ::toupper(c); });
    out += symbol;
    out += "   STO-6G-";
    out += methodName;
    out += "NOCORE\n*\n";
    if (aos.s) {
      out += "    6  s\n";
      for (const auto& gtf : aos.s->gtfs) {
        writeGtfLine(out, gtf);
      }
    }
    if (aos.p && coefficientsAreNonZero(aos.p->gtfs)) {
      out += "    6  p\n";
      for (const auto& gtf : aos.p->gtfs) {
        writeGtfLine(out, gtf);
      }
    }
    if (aos.d && coefficientsAreNonZero(aos.d->gtfs)) {
      out += "    6  d\n";
      for (const auto& gtf : aos.d->gtfs) {
        writeGtfLine(out, gtf);
      }
    }
    return out;
  };

  std::string writeBasisFile(const Parameters& param, std::string& methodName) {
    auto rawParameters = param;
    RawParameterProcessor processor(rawParameters, BasisFunctions::spd);

    std::string basisFileString = "$basis\n";
    Utils::ElementType e;
    for (int i = 1; i < 110; i++) {
      e = Utils::ElementInfo::element(i);
      if (rawParameters.atomic.count(Utils::ElementInfo::Z(e)) == 0) {
        continue;
      }
      auto par = processor.processAtomicParameters(e);
      basisFileString += convertGtoToString(e, par.first->GTOs(), methodName);
    }
    basisFileString += "*\n$end";
    std::ofstream outputFile;
    outputFile.open(methodName + ".basis");
    outputFile << basisFileString;
    outputFile.close();
    return methodName + ".basis";
  };

  bool filesAreIdentical(const std::string& p1, const std::string& p2) {
    std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
      return false; // file problem
    }

    if (f1.tellg() != f2.tellg()) {
      return false; // size mismatch
    }

    // seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()), std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(f2.rdbuf()));
  };
};

TEST_F(Sto6gBasisTest, Am1BasisFileIsCorrect) {
  std::string methodName = "am1";
  std::string fileName = Sto6gBasisTest::writeBasisFile(am1(), methodName);
  ASSERT_TRUE(Sto6gBasisTest::filesAreIdentical(fileName, "AM1-STO-6G.basis"));
  std::remove(fileName.c_str());
}

TEST_F(Sto6gBasisTest, MndoBasisFileIsCorrect) {
  std::string methodName = "mndo";
  std::string fileName = Sto6gBasisTest::writeBasisFile(mndo(), methodName);
  ASSERT_TRUE(Sto6gBasisTest::filesAreIdentical(fileName, "MNDO-STO-6G.basis"));
  std::remove(fileName.c_str());
}

TEST_F(Sto6gBasisTest, Pm3BasisFileIsCorrect) {
  std::string methodName = "pm3";
  std::string fileName = Sto6gBasisTest::writeBasisFile(pm3(), methodName);
  ASSERT_TRUE(Sto6gBasisTest::filesAreIdentical(fileName, "PM3-STO-6G.basis"));
  std::remove(fileName.c_str());
}

TEST_F(Sto6gBasisTest, Pm6BasisFileIsCorrect) {
  std::string methodName = "pm6";
  std::string fileName = Sto6gBasisTest::writeBasisFile(pm6(), methodName);
  ASSERT_TRUE(Sto6gBasisTest::filesAreIdentical(fileName, "PM6-STO-6G.basis"));
  std::remove(fileName.c_str());
}

TEST_F(Sto6gBasisTest, Rm1BasisFileIsCorrect) {
  std::string methodName = "rm1";
  std::string fileName = Sto6gBasisTest::writeBasisFile(rm1(), methodName);
  ASSERT_TRUE(Sto6gBasisTest::filesAreIdentical(fileName, "RM1-STO-6G.basis"));
  std::remove(fileName.c_str());
}

} // namespace Sparrow
} // namespace Scine
