/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_LINEARRESPONSECALCULATOR_H
#define SPARROW_LINEARRESPONSECALCULATOR_H

#include <Core/Interfaces/CalculatorWithReference.h>

namespace Scine {
namespace Sparrow {

class LinearResponseCalculator : public Core::CalculatorWithReference {
 public:
  struct GuessSpecifier {
    Eigen::MatrixXd singlet, triplet, unrestricted;
  };
  ~LinearResponseCalculator() override = default;
  /**
   * @brief Sets the guess to be used in the next calculation. If empty, diagonal dominant guess will be used.
   */
  virtual void setGuess(std::shared_ptr<GuessSpecifier> guessVectorMatrix) = 0;
  /**
   * @brief Returns the guess in the full singles space (no pruning)
   */
  virtual auto getGuess() const -> std::shared_ptr<GuessSpecifier> = 0;
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_LINEARRESPONSECALCULATOR_H
