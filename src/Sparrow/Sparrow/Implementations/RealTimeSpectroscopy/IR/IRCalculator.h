/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_IRCALCULATOR_H
#define SPARROW_IRCALCULATOR_H

#include <Utils/Typenames.h>
#include <memory>
#include <vector>

namespace Scine {
namespace Core {
class Calculator;
class State;
} // namespace Core
namespace Utils {
class AtomCollection;
class Results;
class NormalModesContainer;
class Settings;
} // namespace Utils

namespace Sparrow {
namespace RealTimeSpectroscopy {

class Spectrum;

class IRCalculator {
 public:
  IRCalculator();

  void initialize(const Utils::ElementTypeCollection& elements);
  void modifyPositions(const Utils::PositionCollection& positions);
  /**
   * @brief Calculates an IR spectrum from a set of positions and a gradient.
   * Care must be taken in calling code that the data are thread safe!!
   * @return A non convoluted spectrum to be processed with SpectrumProcessor.
   */
  Spectrum calculate(const Utils::PositionCollection& position, int structureIndex, std::ostream& outs);
  const Utils::Settings& settings() const;
  Utils::Settings& settings();
  void updateState(std::shared_ptr<Core::State> state);
  bool gradientAllowsIRCalculation(const Utils::GradientCollection& gradient) const;
  std::unique_ptr<Utils::AtomCollection> getOptimizedStructure() const;

 private:
  Utils::Results calculateHessianAndDipoleGradient();
  Utils::NormalModesContainer calculateFrequencies(const Utils::HessianMatrix& hessian,
                                                   const Utils::GradientCollection& gradients = {}) const;
  Eigen::VectorXd calculateIntensities(const Utils::DipoleGradient& dipoleGradient, const Eigen::MatrixXd& normalModes) const;
  std::unique_ptr<Utils::PositionCollection> lastPositions_;
  std::unique_ptr<Utils::HessianMatrix> lastHessian_;
  std::unique_ptr<Utils::DipoleGradient> lastDipoleGradient_;
  std::shared_ptr<Core::Calculator> calculator_;
  std::unique_ptr<Utils::Settings> settings_;
};

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_IRCALCULATOR_H
