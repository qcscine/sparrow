/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_NDDOSTATE_H
#define SPARROW_NDDOSTATE_H

#include <Utils/CalculatorBasics/State.h>
#include <exception>
#include <utility>

namespace Scine {

namespace Utils {
class DensityMatrix;
class SpinAdaptedMatrix;
} // namespace Utils

namespace Sparrow {
class GenericMethodWrapper;
/**
 * @brief Exception for the case that a state is requested which is not available.
 */
class StateNotAvailableException : public std::exception {
 public:
  explicit StateNotAvailableException(std::string stateName)
    : stateName_("State " + std::move(stateName) + " is not available.") {
  }
  const char* what() const noexcept final {
    return stateName_.c_str();
  }

 private:
  std::string stateName_;
};
class EmptyStateException : public std::exception {
 public:
  const char* what() const noexcept final {
    return "Sparrow state is empty and cannot be loaded.";
  }
};

/**
 * @brief Definition of a calculation state for NDDO methods.
 * The calculation state is defined as being a collection of the density matrix, optionally the Fock matrix and
 * the one and two electron matrices. Only matrix states are saved here.
 */
struct SparrowState : public Utils::State {
  /**
   * @brief Constructor, calls the base class constructor to initialize the size of the state.
   */
  SparrowState(Utils::StateSize size, GenericMethodWrapper& methodWrapper);

  /**
   * @brief Getter for a Eigen::MatrixXd state identified with a std::string.
   * @param matrixState The key for the state's key-value pair.
   * @return The value of the key-value pair, an Eigen::MatrixXd.
   */
  const Eigen::MatrixXd& getMatrixState(const std::string& matrixState) const final;
  /**
   * @brief Getter for a std::string state identified with a std::string.
   * @param stringState The key for the state's key-value pair.
   * @return The value of the key-value pair, a std::string.
   */
  const std::string& getStringState(const std::string& stringState) const final;
  /**
   * @brief Getter for a integer state identified with a std::string.
   * @param intState The key for the state's key-value pair.
   * @return The value of the key-value pair, an integer.
   */
  int getIntState(const std::string& intState) const final;
  /**
   * @brief Getter for a double state identified with a std::string.
   * @param doubleState The key for the state's key-value pair.
   * @return The value of the key-value pair, a double.
   */
  double getDoubleState(const std::string& doubleState) const final;

  /**
   * @brief Initializer for the state, resets it to the initial state.
   * In this case the initial state is just the initial density matrix guess.
   */
  void initialize() final;

  bool hasState(const std::string& matrixState) const;

  void generateDensityMatrixState(const Utils::DensityMatrix& densityMatrix, bool isUnrestricted);
  void generateFockMatrixState(const Utils::SpinAdaptedMatrix& fockMatrix, bool isUnrestricted);

 private:
  /// @brief A std::map of std::string keys with Eigen::MatrixXd values.
  MatrixState matrixStates_;
  GenericMethodWrapper& methodWrapper_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_NDDOSTATE_H
