/**
 * @file GuessPropagation.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Sparrow/Implementations/TimeDependent/LinearResponseCalculator.h>
#include <Utils/Typenames.h>
#include <deque>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

/**
 * @class GuessPropagation
 * @brief calculates the expansion coefficients for a new set of atomic positions.
 * The coefficients are recovered by a DIIS-type procedure on error vectors
 * calculated by the difference in coordinate.
 * See
 * Accelerating Wave Function Convergence in Interactive QuantumChemical Reactivity Studies
 * Adrian  H.  Mühlbach,  Alain  C.  Vaucher,  and  Markus  Reiher
 * for implementation details.
 * Electronic states vectors are extrapolated instead of density matrices.
 *
 */
class GuessPropagator {
 public:
  /**
   * @brief Calculates the new guess for the excited states calculation.
   * If no previous structures are recorded, an empty vector is returned.
   * @return An Eigen::VectorXd containing the expansion coefficients of the previous positions.
   *         It is empty if no structures were previously recorded.
   */
  auto calculateGuessAtNewPosition(const Utils::PositionCollection& newPositions)
      -> std::shared_ptr<LinearResponseCalculator::GuessSpecifier>;

  /**
   * @brief Saves the excited states of the previous step to be used in the extrapolation.
   */
  void record(const LinearResponseCalculator::GuessSpecifier& excitedState);
  /**
   * @brief Saves the positions of the previous step to be used in the extrapolation.
   */
  void record(const Utils::PositionCollection& positions);

  void setDimension(int diisDimension);

  /**
   * @brief Empties the memory and resets the extrapolation.
   */
  void reset();

 private:
  enum class SpectralType { Singlet, Unrestricted };
  Eigen::MatrixXd calculateOverlapMatrix(const Utils::PositionCollection& positions) const;
  Eigen::MatrixXd calculateErrorVectors(const Utils::PositionCollection& newPosition) const;
  Eigen::VectorXd calculateRhsVector() const;

  std::deque<Eigen::MatrixXd> excitedStatesHistory_;
  std::deque<Utils::PositionCollection> positionsHistory_;
  int maxHistoryDimension_ = 4;
  SpectralType symmetry_;
};

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
