/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "parameters_location.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/GeometricDerivatives/NormalModeAnalyzer.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/FormattedIOUtils.h>
#include <Utils/Math/QuaternionFit.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>

using namespace Scine;

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "Usage: ./program struct1.xyz struct2.xyz" << std::endl;
  }

  std::string filename = argv[1];
  std::ifstream structureFile(filename);
  std::string filename2 = argv[2];
  std::ifstream structureFile2(filename2);

  if (!structureFile.is_open()) {
    throw std::runtime_error("File " + filename + " not accessible.");
  }
  auto structure = Utils::XyzStreamHandler::read(structureFile);
  auto structure2 = Utils::XyzStreamHandler::read(structureFile2);

  auto& manager = Core::ModuleManager::getInstance();

  if (!manager.moduleLoaded("Sparrow")) {
    manager.load("sparrow");
  }

  auto calculator = manager.get<Core::Calculator>("PM6");
  calculator->settings().modifyString(Utils::SettingsNames::parameterRootDirectory, "");
  calculator->settings().modifyString(Utils::SettingsNames::parameterFile, parameters_pm6);
  calculator->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-9);

  calculator->setStructure(structure);

  auto startOptimize1 = std::chrono::system_clock::now();
  Utils::GeometryOptimizer<Utils::Lbfgs> optimizer(*calculator);
  optimizer.optimize(*(calculator->getStructure()));
  auto endOptimize1 = std::chrono::system_clock::now();
  auto timeOptimize1 = std::chrono::duration_cast<std::chrono::milliseconds>(endOptimize1 - startOptimize1).count();

  std::cout << "Structure 1 optimization: " << timeOptimize1 << " ms." << std::endl;
  auto startCalculation = std::chrono::system_clock::now();
  calculator->setRequiredProperties(Utils::Property::Hessian);
  auto results = calculator->calculate("structure1");
  auto endCalculation = std::chrono::system_clock::now();

  auto timeCalculation1 = std::chrono::duration_cast<std::chrono::milliseconds>(endCalculation - startCalculation).count();
  std::cout << "Structure 1 calculation: " << timeCalculation1 << " ms." << std::endl;

  calculator->setStructure(structure2);

  auto startOptimize2 = std::chrono::system_clock::now();
  optimizer.optimize(*(calculator->getStructure()));
  auto endOptimize2 = std::chrono::system_clock::now();
  auto timeOptimize2 = std::chrono::duration_cast<std::chrono::milliseconds>(endOptimize2 - startOptimize2).count();
  std::cout << "Structure 2 optimization: " << timeOptimize2 << " ms." << std::endl;

  calculator->setRequiredProperties(Utils::Property::Hessian);

  auto startCalculation2 = std::chrono::system_clock::now();
  auto results2 = calculator->calculate("structure1");
  auto endCalculation2 = std::chrono::system_clock::now();
  auto timeCalculation2 = std::chrono::duration_cast<std::chrono::milliseconds>(endCalculation2 - startCalculation2).count();
  std::cout << "Structure 2 calculation: " << timeCalculation2 << " ms." << std::endl;

  auto startFrq1 = std::chrono::system_clock::now();
  Utils::NormalModeAnalyzer frequencyAnalyzer(results.get<Utils::Property::Hessian>(), structure.getElements(),
                                              structure.getPositions());
  auto normalModes = frequencyAnalyzer.calculateNormalModes();
  auto endFrq1 = std::chrono::system_clock::now();
  auto timeFrq1 = std::chrono::duration_cast<std::chrono::milliseconds>(endFrq1 - startFrq1).count();
  std::cout << "Structure 1 Diagonalization: " << timeFrq1 << " ms." << std::endl;
  auto startFrq2 = std::chrono::system_clock::now();
  Utils::NormalModeAnalyzer frequencyAnalyzer2(results2.get<Utils::Property::Hessian>(), structure2.getElements(),
                                               structure2.getPositions());
  auto normalModes2 = frequencyAnalyzer2.calculateNormalModes();
  auto endFrq2 = std::chrono::system_clock::now();
  auto timeFrq2 = std::chrono::duration_cast<std::chrono::milliseconds>(endFrq2 - startFrq2).count();
  std::cout << "Structure 2 Diagonalization: " << timeFrq2 << " ms." << std::endl;

  // align positions
  auto startQuaternionFit = std::chrono::system_clock::now();
  Utils::QuaternionFit qFit(structure.getPositions(), structure2.getPositions());
  auto fittedPositions2 = qFit.getFittedData();
  auto endQuaternionFit = std::chrono::system_clock::now();
  auto timeQuaternionFit = std::chrono::duration_cast<std::chrono::milliseconds>(endFrq2 - startFrq2).count();

  std::vector<int> atomIndices(structure.size());
  std::iota(atomIndices.begin(), atomIndices.end(), 0);

  Eigen::VectorXd euclideanDistance = (structure.getPositions() - fittedPositions2).rowwise().norm();
  //  std::vector<int> atomIndices {7, 8, 9, 10, 11, 12, 13, 38,37, 40,41,39, 42};
  std::vector<int> hessianIndices;
  atomIndices.erase(std::remove_if(atomIndices.begin(), atomIndices.end(),
                                   [&](std::size_t index) { return euclideanDistance(index) < 0.3; }),
                    atomIndices.end());

  for (int index : atomIndices) {
    hessianIndices.push_back(3 * index);
    hessianIndices.push_back(3 * index + 1);
    hessianIndices.push_back(3 * index + 2);
  }

  Eigen::MatrixXd hessian1 = results.get<Utils::Property::Hessian>();
  Eigen::MatrixXd hessian2 = results2.get<Utils::Property::Hessian>();
  Eigen::MatrixXd resultHessian = hessian1;

  for (int index : hessianIndices) {
    resultHessian.col(index) = hessian2.col(index);
  }
  Eigen::MatrixXd symmetrizedHessian = resultHessian.selfadjointView<Eigen::Upper>();

  auto startFrq3 = std::chrono::system_clock::now();
  Utils::NormalModeAnalyzer frequencyAnalyzer3(symmetrizedHessian, structure2.getElements(), structure2.getPositions());
  auto normalModes3 = frequencyAnalyzer3.calculateNormalModes();
  auto endFrq3 = std::chrono::system_clock::now();
  auto timeFrq3 = std::chrono::duration_cast<std::chrono::milliseconds>(endFrq3 - startFrq3).count();

  Utils::matrixPrettyPrint(std::cout, normalModes, structure.getElements());
  Utils::matrixPrettyPrint(std::cout, normalModes2, structure.getElements());
  Utils::matrixPrettyPrint(std::cout, normalModes3, structure.getElements());

  std::cout << (structure.getPositions() - fittedPositions2).rowwise().norm() << std::endl;
  std::cout << "Structure 1 geometry optimization: " << timeOptimize1 << " ms." << std::endl;
  std::cout << "Structure 1 Hessian calculation: " << timeCalculation1 << " ms." << std::endl;
  std::cout << "Structure 1 Hessian diagonalization: " << timeFrq1 << " ms." << std::endl;
  std::cout << "Structure 2 geometry optimization: " << timeOptimize2 << " ms." << std::endl;
  std::cout << "Structure 2 Hessian calculation: " << timeCalculation2 << " ms." << std::endl;
  std::cout << "Structure 2 Hessian diagonalization: " << timeFrq2 << " ms." << std::endl;
  std::cout << "Structure 2 quaternion fit: " << timeQuaternionFit << " ms." << std::endl;
  std::cout << "Structure 3 Hessian diagonalization: " << timeFrq3 << " ms." << std::endl;

  return 0;
}
