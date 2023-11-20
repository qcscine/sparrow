/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
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

  SparrowInitializer::initialize();

  try {
    CalculationHandler calculationHandler(commandLineParser);
    calculationHandler.calculate(std::cout);
  }
  catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return 1;
  }

  return 0;
}
