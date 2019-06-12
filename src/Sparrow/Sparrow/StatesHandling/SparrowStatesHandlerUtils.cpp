/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SparrowStatesHandlerUtils.h"
#include "SparrowState.h"
#include <Utils/MethodEssentials/Methods/LCAOMethod.h>
#include <Utils/MethodEssentials/util/DensityMatrix.h>

namespace Scine {
namespace Sparrow {

void SparrowStatesHandlerUtils::loadDensityMatrix(SparrowState& state, Utils::LCAOMethod& method) {
  bool unrestricted = method.unrestrictedCalculationRunning();
  Utils::DensityMatrix densityMatrix;
  auto nElectrons = method.getNumberElectrons();

  // If the density matrix of the receiving state is restricted.
  if (!unrestricted) {
    // If the state to be loaded has a restricted density matrix
    if (state.hasState("Density Matrix")) {
      auto storedDensity = state.getMatrixState("Density Matrix");
      densityMatrix.setDensity(std::move(storedDensity), nElectrons);
    }
    else if (state.hasState("Alpha Density Matrix") && state.hasState("Beta Density Matrix")) {
      // If the state to be loaded was unrestricted (clone() issue)
      auto totDensity = state.getMatrixState("Alpha Density Matrix") + state.getMatrixState("Beta Density Matrix");
      densityMatrix.setDensity(std::move(totDensity), nElectrons);
    }
    else { // No density matrix available!
      throw EmptyStateException();
    }
  }

  // If it is unrestricted now and it had an unrestricted state, then add it to the density matrix.
  if (unrestricted) {
    if (state.hasState("Alpha Density Matrix") && state.hasState("Beta Density Matrix")) {
      auto alphaElectrons = method.getElectronicOccupation().numberAlphaElectrons();
      auto betaElectrons = method.getElectronicOccupation().numberBetaElectrons();
      auto alphaStoredDensity = state.getMatrixState("Alpha Density Matrix");
      auto betaStoredDensity = state.getMatrixState("Beta Density Matrix");
      densityMatrix.setDensity(std::move(alphaStoredDensity), std::move(betaStoredDensity), alphaElectrons, betaElectrons);
    }
    else if (state.hasState("Density Matrix")) {
      auto storedDensity = state.getMatrixState("Density Matrix");
      densityMatrix.setDensity(std::move(storedDensity), nElectrons);
    }
    else {
      throw EmptyStateException();
    }
  }
  method.setDensityMatrix(std::move(densityMatrix));
}
} // namespace Sparrow
} // namespace Scine
