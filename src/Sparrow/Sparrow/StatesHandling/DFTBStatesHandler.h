/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTBSTATESHANDLER_H
#define SPARROW_DFTBSTATESHANDLER_H

#include <Utils/CalculatorBasics/StatesHandler.h>

namespace Scine {
namespace Utils {
class State;
class DensityMatrix;
} // namespace Utils
namespace Sparrow {
class DFTBMethodWrapper;
class SparrowState;

class DFTBStatesHandler : public Utils::StatesHandler {
 public:
  /**
   * @brief Constructor, it takes the embedding method wrapper as argument in order to read the state from it.
   * @param methodWrapper The embedding method wrapper.
   */
  explicit DFTBStatesHandler(DFTBMethodWrapper& methodWrapper);
  ~DFTBStatesHandler() final;
  /**
   * @brief The implementation of the store function. It stores a state with the specified size.
   * @param size The desired Utils::StateSize of the state to store.
   */
  void store(Utils::StateSize size) final;
  /**
   * @brief Loads a state in the embedding NDDOMethodWrapper
   * @param state The state to load. If it is not compatible with the NDDOMethodWrapper,
   *        a StateNotCompatibleException is thrown.
   */
  void load(std::shared_ptr<Utils::State> state) final;
  /**
   * @brief Returns the current state of the DFTB method.
   * @param The desired Utils::StateSize of the state to store.
   * This method does not need to store the current state in the stateHandler,
   * it is a const method and therefore suited to be used in a copy-constructor.
   */
  std::shared_ptr<Utils::State> getCurrentState(Utils::StateSize size) const final;

 private:
  DFTBMethodWrapper& methodWrapper_;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DFTBSTATESHANDLER_H
