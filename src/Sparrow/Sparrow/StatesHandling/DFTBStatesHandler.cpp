/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "DFTBStatesHandler.h"
#include "SparrowStatesHandlerUtils.h"
#include <Sparrow/Implementations/Dftb/DFTBMethodWrapper.h>
#include <Sparrow/StatesHandling/SparrowState.h>
#include <Utils/MethodEssentials/Methods/LCAOMethod.h>

namespace Scine {
namespace Sparrow {

DFTBStatesHandler::DFTBStatesHandler(DFTBMethodWrapper& methodWrapper) : methodWrapper_(methodWrapper) {
}

DFTBStatesHandler::~DFTBStatesHandler() = default;

void DFTBStatesHandler::store(Utils::StateSize size) {
  states_.emplace_back(std::dynamic_pointer_cast<SparrowState>(getCurrentState(size)));
}

void DFTBStatesHandler::load(std::shared_ptr<Utils::State> state) {
  if (!state) {
    throw EmptyStateException();
  }
  try {
    std::shared_ptr<SparrowState> sparrowState = std::dynamic_pointer_cast<SparrowState>(state);
    SparrowStatesHandlerUtils::loadDensityMatrix(*sparrowState, methodWrapper_.getLCAOMethod());
  }
  catch (std::bad_cast& e) {
    throw StateNotCompatibleException();
  }
}

std::shared_ptr<Utils::State> DFTBStatesHandler::getCurrentState(Utils::StateSize size) const {
  SparrowState state(size, methodWrapper_);
  bool unrestricted = methodWrapper_.getLCAOMethod().unrestrictedCalculationRunning();
  const auto& densityMatrix = methodWrapper_.getLCAOMethod().getDensityMatrix();

  if (size == Utils::StateSize::minimal) {
    state.generateDensityMatrixState(densityMatrix, unrestricted);
  }
  else if (size == Utils::StateSize::regular || size == Utils::StateSize::extensive) {
    state.generateDensityMatrixState(densityMatrix, unrestricted);
    state.generateFockMatrixState(methodWrapper_.getLCAOMethod().getFockMatrix(), unrestricted);
  }
  return std::make_shared<SparrowState>(state);
}

} // namespace Sparrow
} // namespace Scine
