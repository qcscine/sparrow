/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/RealTimeSpectroscopy/IR/IRCalculator.h>
#include <Sparrow/Implementations/RealTimeSpectroscopy/Utils/Spectrum.h>
#include <Sparrow/Implementations/RealTimeSpectroscopy/UvVis/UvVisCalculator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/IO/Yaml.h>
#include <Utils/MolecularTrajectory.h>
#include <yaml-cpp/yaml.h>
#include <chrono>
#include <iomanip>
#include <iostream>

using namespace Scine;
using namespace Sparrow;
using namespace RealTimeSpectroscopy;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: ./program input.yaml" << std::endl;
    return -1;
  }

  // Load the input file
  std::string yamlInputName = argv[1];

  auto input = YAML::LoadFile(yamlInputName);
  if (!input["file"])
    std::cout << "Missing file field." << std::endl;
  std::string xyzName, methodName;
  xyzName = input["file"].as<std::string>();

  if (!input["method_name"])
    std::cout << "Missing method_name field." << std::endl;

  methodName = input["method_name"].as<std::string>();

  std::cout << "method name: " << methodName << std::endl;
  auto& manager = Core::ModuleManager::getInstance();

  std::ifstream structureFile(xyzName);
  if (!structureFile.is_open()) {
    throw std::runtime_error("File " + xyzName + " not accessible.");
  }

  // Parse the trajectory and set up the calculator
  auto trajectory = std::make_unique<Utils::MolecularTrajectory>(
      Utils::MolecularTrajectoryIO::read(Utils::MolecularTrajectoryIO::format::xyz, structureFile));
  auto irCalc = std::make_unique<IRCalculator>();
  auto uvVisCalc = std::make_unique<UvVisCalculator>();
  auto singlePointCalculator = manager.get<Core::Calculator>(methodName);

  // Parse and set the settings for the calculators
  if (auto a = input["single_point_settings"])
    Utils::nodeToSettings(singlePointCalculator->settings(), input["single_point_settings"]);
  else
    std::cout << "No single point settings" << std::endl;
  if (auto a = input["settings_uv_vis"])
    Utils::nodeToSettings(uvVisCalc->settings(), input["settings_uv_vis"]);
  else
    std::cout << "No UV/Vis settings" << std::endl;
  if (auto a = input["settings_ir"])
    Utils::nodeToSettings(irCalc->settings(), input["settings_ir"]);
  else
    std::cout << "No settings for IR" << std::endl;

  // Initialize the calculators
  irCalc->initialize(trajectory->getElementTypes());
  uvVisCalc->initialize(trajectory->getElementTypes());
  singlePointCalculator->setStructure({trajectory->getElementTypes(), *trajectory->begin()});
  singlePointCalculator->setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients);

  std::vector<double> irTimings, uvVisTimings;
  irTimings.reserve(trajectory->size());
  uvVisTimings.reserve(trajectory->size());
  int uvVisFrequency = input["uv_vis_frequency"].as<int>();
  if (uvVisFrequency <= 0)
    throw std::runtime_error("Uv/Vis frequency must be an integer > 0!");

  std::cout << "Setup complete!" << std::endl;
  int counter = 0;
  bool expectingNewMinimum = true;

  std::ofstream timingsUVVis("uvVisTimings.out"), timingsIR("irTimings.out");
  std::ofstream gsProperties("exploration.out");
  gsProperties << std::right << std::scientific << std::setw(15) << "Structure" << std::setw(15) << "Energies"
               << std::setw(15) << "Gradients" << std::endl;
  timingsIR << "Optimiziation  Hessian" << std::endl;

  for (const auto& structure : *trajectory) {
    std::cout << "Digesting structure " << ++counter << "/" << trajectory->size() << "." << std::endl;
    singlePointCalculator->modifyPositions(structure);
    singlePointCalculator->calculate();
    gsProperties << std::right << std::setw(15) << counter << std::setw(15)
                 << singlePointCalculator->results().get<Utils::Property::Energy>() << std::setw(15)
                 << singlePointCalculator->results().get<Utils::Property::Gradients>().rowwise().norm().sum() << std::endl;

    auto state = singlePointCalculator->getState();
    irCalc->updateState(state);
    uvVisCalc->updateState(std::move(state));
    if (irCalc->gradientAllowsIRCalculation(singlePointCalculator->results().get<Utils::Property::Gradients>())) {
      if (expectingNewMinimum) {
        std::ofstream spectraIR("irSpectra" + std::to_string(counter) + ".out");
        std::ofstream structureOut("structure_opt_" + std::to_string(counter) + ".xyz");
        auto irSpectrum = irCalc->calculate(structure, counter, timingsIR);
        for (int i = 0; i < irSpectrum.size(); ++i)
          spectraIR << irSpectrum.getXData(i) << " " << irSpectrum.getYData(i) << std::endl;
        Utils::XyzStreamHandler::write(structureOut, *irCalc->getOptimizedStructure());
        structureOut.close();
      }
      expectingNewMinimum = false;
    }
    else {
      expectingNewMinimum = true;
    }

    if (counter % uvVisFrequency == 0) {
      std::ofstream spectraUVVis("uvVisSpectra" + std::to_string(counter) + ".out");
      auto startUvVis = std::chrono::system_clock::now();
      auto uvVisSpectrum = uvVisCalc->calculate(structure);

      for (int i = 0; i < uvVisSpectrum.size(); ++i)
        spectraUVVis << uvVisSpectrum.getXData(i) << " " << uvVisSpectrum.getYData(i) << std::endl;

      auto endUvVis = std::chrono::system_clock::now();
      double uvVisTiming = std::chrono::duration_cast<std::chrono::milliseconds>(endUvVis - startUvVis).count();
      timingsUVVis << uvVisTiming << std::endl;
    }
  }

  gsProperties.close();

  return 0;
}
