/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_UVVISCALCULATOR_H
#define SPARROW_UVVISCALCULATOR_H

#include <Sparrow/Implementations/TimeDependent/LinearResponseCalculator.h>
#include <Utils/Typenames.h>
#include <fstream>
#include <memory>

namespace Scine {
namespace Core {
class CalculatorWithReference;
class State;
} // namespace Core
namespace Utils {
class Settings;
class Results;
} // namespace Utils
namespace Sparrow {
namespace RealTimeSpectroscopy {

class Spectrum;
class GuessPropagator;

class UvVisCalculator {
 public:
  UvVisCalculator();
  explicit UvVisCalculator(std::shared_ptr<Utils::Settings> settings);
  ~UvVisCalculator();

  void initialize(const Utils::ElementTypeCollection& elements);
  void modifyPositions(const Utils::PositionCollection& positions);
  /**
   * @brief Calculates an UV/VIS spectrum from a set of positions.
   * Care must be taken in calling code that the data are thread safe!!
   * @return A non convoluted spectrum to be processed with SpectrumProcessor.
   */
  Spectrum calculate(const Utils::PositionCollection& positions);
  void updateState(std::shared_ptr<Core::State> state);
  Utils::Settings& settings();
  const Utils::Settings& settings() const;

 private:
  void createReferenceCalculator();
  void createExcitedStatesCalculator();
  std::shared_ptr<LinearResponseCalculator> calculator_;
  std::unique_ptr<GuessPropagator> guessPropagator_;
  std::unique_ptr<std::ofstream> sizeOutput_;
  std::shared_ptr<LinearResponseCalculator::GuessSpecifier> lastExcitedStates_;
  std::shared_ptr<Utils::Settings> settings_;
};

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_UVVISCALCULATOR_H
