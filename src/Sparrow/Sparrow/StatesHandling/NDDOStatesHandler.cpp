/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NDDOStatesHandler.h"
#include "SparrowStatesHandlerUtils.h"
#include <Sparrow/Implementations/Nddo/NDDOMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Utils/OneElectronMatrix.h>
#include <Sparrow/Implementations/Nddo/Utils/TwoElectronMatrix.h>
#include <Sparrow/StatesHandling/SparrowState.h>
#include <Utils/MethodEssentials/Methods/LCAOMethod.h>
#include <Utils/MethodEssentials/util/DensityMatrix.h>

namespace Scine {
namespace Sparrow {

NDDOStatesHandler::NDDOStatesHandler(NDDOMethodWrapper& methodWrapper)
  : StatesHandler(), methodWrapper_(methodWrapper) {
}

NDDOStatesHandler::~NDDOStatesHandler() = default;

void Scine::Sparrow::NDDOStatesHandler::store(Utils::StateSize size) {
  states_.emplace_back(std::dynamic_pointer_cast<SparrowState>(getCurrentState(size)));
}

void NDDOStatesHandler::load(std::shared_ptr<Utils::State> state) {
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

std::shared_ptr<Utils::State> NDDOStatesHandler::getCurrentState(Utils::StateSize size) const {
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
