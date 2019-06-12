/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "SparrowState.h"
#include <Sparrow/Implementations/GenericMethodWrapper.h>
#include <Utils/MethodEssentials/util/DensityMatrix.h>
#include <Utils/Settings.h>

namespace Scine {
namespace Sparrow {

SparrowState::SparrowState(Scine::Utils::StateSize size, GenericMethodWrapper& methodWrapper)
  : Utils::State(size), methodWrapper_(methodWrapper) {
}

const Eigen::MatrixXd& SparrowState::getMatrixState(const std::string& matrixState) const {
  auto it = matrixStates_.find(matrixState);
  if (it != matrixStates_.end()) {
    return it->second;
  }
  else {
    throw StateNotAvailableException(matrixState);
  }
}

const std::string& SparrowState::getStringState(const std::string& stringState) const {
  throw StateNotAvailableException(stringState);
}

int SparrowState::getIntState(const std::string& intState) const {
  throw StateNotAvailableException(intState);
}

double SparrowState::getDoubleState(const std::string& doubleState) const {
  throw StateNotAvailableException(doubleState);
}

void SparrowState::generateFockMatrixState(const Utils::SpinAdaptedMatrix& fockMatrix, bool isUnrestricted) {
  if (!isUnrestricted)
    matrixStates_.insert({"Fock Matrix", fockMatrix.restrictedMatrix()});
  else {
    matrixStates_.insert({"Alpha Fock Matrix", fockMatrix.alphaMatrix()});
    matrixStates_.insert({"Beta Fock Matrix", fockMatrix.betaMatrix()});
  }
}

void SparrowState::generateDensityMatrixState(const Utils::DensityMatrix& densityMatrix, bool isUnrestricted) {
  if (!isUnrestricted) {
    matrixStates_.insert({"Density Matrix", densityMatrix.restrictedMatrix()});
  }
  else {
    matrixStates_.insert({"Alpha Density Matrix", densityMatrix.alphaMatrix()});
    matrixStates_.insert({"Beta Density Matrix", densityMatrix.betaMatrix()});
  }
}

void SparrowState::initialize() {
  auto densityGuess = methodWrapper_.getDensityMatrixGuess();
  matrixStates_.clear();
  generateDensityMatrixState(densityGuess, false);
}

bool SparrowState::hasState(const std::string& matrixState) const {
  return matrixStates_.find(matrixState) != matrixStates_.end();
}

} // namespace Sparrow
} // namespace Scine
