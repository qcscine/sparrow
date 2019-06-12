/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CalculationHandler.h"
#include "CommandLineOptions.h"
#include "SparrowInitializer.h"
#include <iostream>

using namespace Scine::Core;
using namespace Scine::Utils;
using namespace Scine::Sparrow;

int main(int argc, char* argv[]) {
  CommandLineOptions commandLineParser(argc, argv);

  if (commandLineParser.helpRequired() || argc == 1) {
    commandLineParser.printHelp(std::cout);
    return 0;
  }

  SparrowInitializer initializer;
  initializer.initialize();

  try {
    CalculationHandler calculationHandler(commandLineParser, initializer);
    calculationHandler.calculate(std::cout);
  }
  catch (MethodNotAvailableException& e) {
    std::cout << e.what() << std::endl;
  }

  return 0;
}
