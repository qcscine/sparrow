/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "parameters_location.h"
#include <App/SPARROWInitializer.h>
#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Settings.h>
#include <boost/dll/runtime_symbol_info.hpp>
#include <iostream>

int main(int argc, char* argv[]) {
  Scine::Sparrow::SPARROWInitializer initializer;
  initializer.initialize();

  auto& manager = Scine::Core::ModuleManager::getInstance();

  auto programPath = boost::dll::program_location().parent_path().parent_path();

  auto libraryPath = boost::filesystem::path{programPath / "Sparrow" / "sparrow"};

  std::cout << "Loading: " << libraryPath.string() << std::endl;

  manager.load(libraryPath);

  std::cout << "Loaded library: " << manager.getLoadedModuleNames()[0] << std::endl;

  // List all interfaces and models
  const auto interfaces = manager.getLoadedInterfaces();
  for (const auto& interfaceName : interfaces) {
    std::cout << interfaceName << "\n";
    for (const auto& modelName : manager.getLoadedModels(interfaceName)) {
      std::cout << "- " << modelName << "\n";
    }
  }

  auto stream = std::stringstream("5\n\n"
                                  "C      0.0000000000    0.0000000000    0.0000000000\n"
                                  "H      0.6287000000    0.6287000000    0.6287000000\n"
                                  "H     -0.6287000000   -0.6287000000    0.6287000000\n"
                                  "H     -0.6287000000    0.6287000000   -0.6287000000\n"
                                  "H      0.6287000000   -0.6287000000   -0.6287000000\n");
  /*  auto stream = std::stringstream("74\n\n"
                                    "C         16.90467       -2.36992      -15.89125\n"
                                    "C         17.73850       -1.98799      -14.67304\n"
                                    "C         19.21439       -2.33779      -14.98968\n"
                                    "C         17.62702       -0.49442      -14.27738\n"
                                    "C         16.23611        0.10697      -13.92977\n"
                                    "C         15.61241       -0.18736      -12.52137\n"
                                    "C         14.29698        0.66220      -12.22287\n"
                                    "C         13.25492        0.17585      -13.30259\n"
                                    "C         13.53121        0.55271      -10.81340\n"
                                    "C         13.33255       -0.95092      -10.51438\n"
                                    "C         13.77447       -1.21087       -9.10124\n"
                                    "C         14.30720        0.08896       -8.51975\n"
                                    "C         13.63805        1.21554       -9.30999\n"
                                    "C         14.39142        2.49784       -8.81981\n"
                                    "C         14.42315        2.66500       -7.19117\n"
                                    "C         14.81161        1.38124       -6.33899\n"
                                    "C         14.09619        0.16687       -6.99865\n"
                                    "C         14.59821       -1.15071       -6.38014\n"
                                    "C         14.26266       -1.18332       -4.95720\n"
                                    "C         14.27624       -0.11311       -4.16029\n"
                                    "C         14.61677        1.33208       -4.66821\n"
                                    "C         15.95442        1.66369       -3.89644\n"
                                    "C         15.88919        1.37681       -2.36536\n"
                                    "C         15.63792       -0.10847       -2.11160\n"
                                    "C         14.26542       -0.45416       -2.68056\n"
                                    "O         15.73105       -0.43082       -0.72916\n"
                                    "C         13.57880        2.37974       -4.18354\n"
                                    "C         12.13707        1.45210       -8.88864\n"
                                    "H         17.01846       -3.44140      -16.10804\n"
                                    "H         17.22662       -1.82944      -16.79353\n"
                                    "H         15.83161       -2.18010      -15.73179\n"
                                    "H         17.42327       -2.60619      -13.82530\n"
                                    "H         19.58093       -1.77083      -15.85370\n"
                                    "H         19.32058       -3.40102      -15.23299\n"
                                    "H         19.86582       -2.12127      -14.13375\n"
                                    "H         18.00635        0.08304      -15.13350\n"
                                    "H         18.31282       -0.29071      -13.44925\n"
                                    "H         15.54400       -0.17696      -14.72891\n"
                                    "H         16.34626        1.19636      -14.00025\n"
                                    "H         16.34912        0.01456      -11.74111\n"
                                    "H         15.39295       -1.26013      -12.50605\n"
                                    "H         14.50219        1.72020      -12.41675\n"
                                    "H         13.09268       -0.91006      -13.23635\n"
                                    "H         13.57544        0.38774      -14.32528\n"
                                    "H         12.27638        0.65387      -13.18441\n"
                                    "H         12.54354        0.92803      -11.12164\n"
                                    "H         13.88943       -1.61900      -11.17131\n"
                                    "H         12.27838       -1.22740      -10.63160\n"
                                    "H         12.92119       -1.59696       -8.53720\n"
                                    "H         14.56883       -1.96809       -9.07914\n"
                                    "H         15.38666        0.13314       -8.72595\n"
                                    "H         13.91573        3.39705       -9.23559\n"
                                    "H         15.42553        2.49767       -9.18063\n"
                                    "H         13.43533        3.02740       -6.88465\n"
                                    "H         15.11394        3.48361       -6.94709\n"
                                    "H         15.89085        1.25460       -6.50115\n"
                                    "H         13.02729        0.25492       -6.78129\n"
                                    "H         15.66604       -1.26970       -6.52831\n"
                                    "H         14.09287       -2.01790       -6.81809\n"
                                    "H         14.10039       -2.17532       -4.54740\n"
                                    "H         16.23526        2.71299       -4.03807\n"
                                    "H         16.77073        1.05968       -4.31224\n"
                                    "H         15.08896        1.95648       -1.88508\n"
                                    "H         16.81400        1.68277       -1.86510\n"
                                    "H         16.41469       -0.69890       -2.61365\n"
                                    "H         13.44929        0.03485       -2.14085\n"
                                    "H         14.10163       -1.52946       -2.51003\n"
                                    "H         14.88081       -0.23983       -0.30184\n"
                                    "H         12.77327        2.54740       -4.90032\n"
                                    "H         14.05774        3.35286       -4.00237\n"
                                    "H         13.10310        2.09735       -3.23695\n"
                                    "H         11.53095        0.54614       -8.99214\n"
                                    "H         11.68478        2.22951       -9.51992\n"
                                    "H         12.00940        1.79465       -7.85873");*/
  auto structure = Scine::Utils::XyzStreamHandler::read(stream);
  // Load PM3Method and adapt settings
  auto pm6 = manager.get<Scine::Core::Calculator>("PM3");
  pm6->settings->modifyString("parameter_root", "");
  pm6->settings->modifyString("parameter_file", parameters_pm3);
  pm6->setStructure(structure);
  auto res1 = pm6->calculate("DFTB3 first");
  std::cout << "Energy of calculation '" << res1.getDescription() << "' = " << res1.get<Utils::Property::Energy>()
            << std::endl;
  std::cout << "Changing Settings ... \n";
  pm6->settings->modifyInt("molecular_charge", 0);
  pm6->settings->modifyBool("unrestricted_calculation", false);
  pm6->settings->modifyInt("spin_multiplicity", 1);
  pm6->settings->modifyInt("max_iterations", 1000);
  pm6->settings->modifyDouble("self_consistence_criterion", 1e-8);
  //  pm6->settings->modifyString("SCF_mixer", "diis_mixer");
  pm6->setRequiredProperties(Scine::Utils::Property::Energy | Scine::Utils::Property::Gradients |
                             Scine::Utils::Property::Hessian);
  if (pm6->settings->check()) {
    auto result = pm6->calculate("DFTB3 Second");
    std::cout << "Energy of calculation '" << result.getDescription() << "' = " << result.get<Utils::Property::Energy>()
              << std::endl;
    std::cout << "And gradients = " << result.get<Utils::Property::Gradients>() << std::endl;
    std::cout << "And Hessian:" << std::endl;
    std::cout << result.get<Utils::Property::Hessian>() << std::endl;
  }
  else {
    std::cout << "Invalid settings!" << std::endl;
  }
  return 0;
}
