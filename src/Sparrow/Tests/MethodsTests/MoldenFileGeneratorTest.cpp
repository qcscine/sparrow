/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "parameters_location.h"
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/MoldenFileGenerator.h>
#include <Sparrow/Implementations/Nddo/Pm6/Wrapper/PM6MethodWrapper.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <gmock/gmock.h>

using namespace testing;

namespace Scine {
namespace Sparrow {
namespace Tests {

class AMoldenFileGeneratorTest : public Test {
 public:
  std::shared_ptr<GenericMethodWrapper> calculator, calculatorFeO;
  Eigen::VectorXd mosFeO;

 protected:
  bool headerPresent(const std::string& filename) {
    std::ifstream molden("moldenTest.molden");
    std::string buffer;
    bool result = false;
    while (std::getline(molden, buffer)) {
      if (buffer.find("[Molden Format]") != std::string::npos) {
        std::getline(molden, buffer);
        result = buffer.find("Written by Sparrow, the semiempirical library of the SCINE software package") !=
                 std::string::npos;
      }
    }
    return result;
  }

 private:
  void SetUp() final {
    std::stringstream h2o{"3\n\n"
                          "O 0 0 -0.07\n"
                          "H 0 -0.5 0.45\n"
                          "H 0 0.5 0.45\n"};

    std::stringstream feo{"2\n\n"
                          "Fe      0.0 0.0 0.0\n"
                          "O      0.8 0.80 0.000000"};
    auto structure = Utils::XyzStreamHandler::read(h2o);

    calculator = std::make_shared<PM6MethodWrapper>();
    calculator->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, parameters_root);
    calculator->setStructure(structure);
    calculatorFeO = calculator->clone();
    calculatorFeO->setStructure(Utils::XyzStreamHandler::read(feo));
    mosFeO = Eigen::VectorXd(13);
    mosFeO << -0.1894410376, -0.0293471793, -0.0293471793, 0.0000000000, 0.2036461848, -0.0000000000, 0.0000000000,
        -0.0000000000, -0.3527255388, -0.8901491672, -0.0454805742, -0.0454805742, 0.0000000000;
  }

  void TearDown() final {
    std::remove("moldenTest.molden");
  }
};

TEST_F(AMoldenFileGeneratorTest, CanGenerateGTOAtomAndMOsForWavefunctionH2O) {
  calculator->calculate("");
  MoldenFileGenerator wfGen(*calculator);
  std::ofstream output("moldenTest.molden");
  if (output.is_open()) {
    wfGen.generateWavefunctionInformation(output);
  }
  {
    std::ifstream molden("moldenTest.molden");
    ASSERT_TRUE(molden.is_open());
    std::string buffer;
    bool fiveDPresent = false, sevenFPresent = false, nineGPresent = false;
    while (std::getline(molden, buffer)) {
      if (buffer.find("[Atoms]") != std::string::npos) {
        ASSERT_TRUE(buffer.find("(AU)"));
        std::string el;
        int index, atomicNumber;
        double x, y, z;
        std::getline(molden, buffer);
        std::stringstream line(buffer);
        line >> el >> index >> atomicNumber >> x >> y >> z;
        ASSERT_EQ(el, "O");
        ASSERT_EQ(index, 1);
        ASSERT_EQ(atomicNumber, 8);
        ASSERT_THAT(x, DoubleNear(0.0, 1e-9));
        ASSERT_THAT(y, DoubleNear(0.0, 1e-9));
        ASSERT_THAT(z, DoubleNear(-0.1322808288, 1e-9));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> el >> index >> atomicNumber >> x >> y >> z;
        ASSERT_EQ(el, "H");
        ASSERT_EQ(index, 2);
        ASSERT_EQ(atomicNumber, 1);
        ASSERT_THAT(x, DoubleNear(0.0, 1e-9));
        ASSERT_THAT(y, DoubleNear(-0.9448630627, 1e-9));
        ASSERT_THAT(z, DoubleNear(0.8503767565, 1e-9));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> el >> index >> atomicNumber >> x >> y >> z;
        ASSERT_EQ(el, "H");
        ASSERT_EQ(index, 3);
        ASSERT_EQ(atomicNumber, 1);
        ASSERT_THAT(x, DoubleNear(0.0, 1e-9));
        ASSERT_THAT(y, DoubleNear(0.9448630627, 1e-9));
        ASSERT_THAT(z, DoubleNear(0.8503767565, 1e-9));
      }
      if (buffer.find("[GTO]") != std::string::npos) {
        std::string index, zero, orbType, howMany, oneZeroZero;
        double exp, coeff;
        std::getline(molden, buffer);
        std::stringstream line(buffer);
        line >> index >> zero;
        ASSERT_EQ(index, "1");
        ASSERT_EQ(zero, "0");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> orbType >> howMany >> oneZeroZero;
        ASSERT_EQ(orbType, "s");
        ASSERT_EQ(howMany, "6");
        ASSERT_EQ(oneZeroZero, "1.00");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> exp >> coeff;
        ASSERT_THAT(exp, DoubleNear(813.8100984652, 1.0e-10));
        ASSERT_THAT(coeff, DoubleNear(-0.4508003008, 1.0e-10));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> exp >> coeff;
        ASSERT_THAT(exp, DoubleNear(149.2444978754, 1.0e-10));
        ASSERT_THAT(coeff, DoubleNear(-0.6290416320, 1.0e-10));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> exp >> coeff;
        ASSERT_THAT(exp, DoubleNear(41.9409236915, 1.0e-10));
        ASSERT_THAT(coeff, DoubleNear(-0.6049526644, 1.0e-10));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> exp >> coeff;
        ASSERT_THAT(exp, DoubleNear(5.9976452051, 1.0e-10));
        ASSERT_THAT(coeff, DoubleNear(0.9140223237, 1.0e-10));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> exp >> coeff;
        ASSERT_THAT(exp, DoubleNear(2.7221002652, 1.0e-10));
        ASSERT_THAT(coeff, DoubleNear(0.8489976775, 1.0e-10));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> exp >> coeff;
        ASSERT_THAT(exp, DoubleNear(1.2981542343, 1.0e-10));
        ASSERT_THAT(coeff, DoubleNear(0.1484775632, 1.0e-10));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> orbType >> howMany >> oneZeroZero;
        ASSERT_EQ(orbType, "p");
        ASSERT_EQ(howMany, "6");
        ASSERT_EQ(oneZeroZero, "1.00");
      }
      if (buffer.find("[5D]") != std::string::npos)
        fiveDPresent = true;
      if (buffer.find("[7F]") != std::string::npos)
        sevenFPresent = true;
      if (buffer.find("[9G]") != std::string::npos)
        nineGPresent = true;

      if (buffer.find("[MO]") != std::string::npos) {
        std::string description, stringEntry;
        double doubleEntry;
        std::getline(molden, buffer);
        std::stringstream line(buffer);
        line >> description >> stringEntry;
        ASSERT_EQ(description, "Sym=");
        ASSERT_EQ(stringEntry, "A1");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> doubleEntry;
        ASSERT_EQ(description, "Ene=");
        ASSERT_THAT(doubleEntry, DoubleNear(-1.3934000618, 1e-9));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> stringEntry;
        ASSERT_EQ(description, "Spin=");
        ASSERT_EQ(stringEntry, "Alpha");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> doubleEntry;
        ASSERT_EQ(description, "Occup=");
        ASSERT_THAT(doubleEntry, DoubleNear(2.0, 1e-9));
      }
    }
    ASSERT_TRUE(fiveDPresent && sevenFPresent && nineGPresent);
  }

  ASSERT_TRUE(headerPresent("moldenTest.molden"));
}

TEST_F(AMoldenFileGeneratorTest, CanGenerateWavefunctionThroughInterfaceH2O) {
  calculator->calculate("");

  auto& manager = Core::ModuleManager::getInstance();
  if (!manager.moduleLoaded("Sparrow"))
    manager.load("Sparrow");
  auto wfGen = manager.get<Core::WavefunctionOutputGenerator>("PM6");
  wfGen->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, parameters_root);
  wfGen->setStructure(*calculator->getStructure());
  wfGen->loadState(calculator->getState());
  wfGen->generateWavefunctionInformation("moldenTest.molden");
}

TEST_F(AMoldenFileGeneratorTest, CanGenerateMOsForWavefunctionFeO) {
  calculatorFeO->calculate("");
  MoldenFileGenerator wfGen(*calculatorFeO);
  std::ofstream output("moldenTest.molden");
  if (output.is_open()) {
    wfGen.generateWavefunctionInformation(output);
  }
  {
    std::ifstream molden("moldenTest.molden");
    ASSERT_TRUE(molden.is_open());
    std::string buffer;
    std::stringstream line;
    while (std::getline(molden, buffer)) {
      if (buffer.find("[MO]") != std::string::npos) {
        std::string description, stringEntry;
        double doubleEntry;
        std::getline(molden, buffer);
        std::stringstream line(buffer);
        line >> description >> stringEntry;
        ASSERT_EQ(description, "Sym=");
        ASSERT_EQ(stringEntry, "A1");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> doubleEntry;
        ASSERT_EQ(description, "Ene=");
        ASSERT_THAT(doubleEntry, DoubleNear(-1.1239734808, 1e-9));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> stringEntry;
        ASSERT_EQ(description, "Spin=");
        ASSERT_EQ(stringEntry, "Alpha");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> doubleEntry;
        ASSERT_EQ(description, "Occup=");
        ASSERT_THAT(doubleEntry, DoubleNear(2.0, 1e-9));
        for (int i = 1; i <= 13; ++i) {
          std::getline(molden, buffer);
          line = std::stringstream(buffer);
          line >> description >> doubleEntry;
          ASSERT_EQ(description, std::to_string(i));
          ASSERT_THAT(std::fabs(doubleEntry), DoubleNear(std::fabs(mosFeO(i - 1)), 1e-9));
        }
      }
    }
  }

  ASSERT_TRUE(headerPresent("moldenTest.molden"));
}

TEST_F(AMoldenFileGeneratorTest, CanGenerateMOsForUnrestrictedWavefunctionFeO) {
  calculatorFeO->settings().modifyBool(Utils::SettingsNames::unrestrictedCalculation, true);
  calculatorFeO->calculate("");
  MoldenFileGenerator wfGen(*calculatorFeO);
  std::ofstream output("moldenTest.molden");
  if (output.is_open()) {
    wfGen.generateWavefunctionInformation(output);
  }
  {
    std::ifstream molden("moldenTest.molden");
    ASSERT_TRUE(molden.is_open());
    std::string buffer;
    std::stringstream line;
    while (std::getline(molden, buffer)) {
      if (buffer.find("[MO]") != std::string::npos) {
        std::string description, stringEntry;
        double doubleEntry;
        std::getline(molden, buffer);
        std::stringstream line(buffer);
        line >> description >> stringEntry;
        ASSERT_EQ(description, "Sym=");
        ASSERT_EQ(stringEntry, "A1");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> doubleEntry;
        ASSERT_EQ(description, "Ene=");
        ASSERT_THAT(doubleEntry, DoubleNear(-1.1239734808, 1e-9));
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> stringEntry;
        ASSERT_EQ(description, "Spin=");
        ASSERT_EQ(stringEntry, "Alpha");
        std::getline(molden, buffer);
        line = std::stringstream(buffer);
        line >> description >> doubleEntry;
        ASSERT_EQ(description, "Occup=");
        ASSERT_THAT(doubleEntry, DoubleNear(1.0, 1e-9));
        for (int i = 1; i <= 13; ++i) {
          std::getline(molden, buffer);
          line = std::stringstream(buffer);
          line >> description >> doubleEntry;
          ASSERT_EQ(description, std::to_string(i));
          ASSERT_THAT(std::fabs(doubleEntry), DoubleNear(std::fabs(mosFeO(i - 1)), 1e-9));
        }
      }
    }
  }

  ASSERT_TRUE(headerPresent("moldenTest.molden"));
}

TEST_F(AMoldenFileGeneratorTest, CanGenerateWavefunctionThroughInterfaceFeO) {
  calculatorFeO->calculate("");

  auto& manager = Core::ModuleManager::getInstance();
  if (!manager.moduleLoaded("Sparrow"))
    manager.load("Sparrow");
  auto wfGen = manager.get<Core::WavefunctionOutputGenerator>("PM6");
  wfGen->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, parameters_root);
  wfGen->setStructure(*calculatorFeO->getStructure());
  wfGen->loadState(calculatorFeO->getState());
  wfGen->generateWavefunctionInformation("moldenTest.molden");
}

} // namespace Tests
} // namespace Sparrow
} // namespace Scine
