/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "NDDODensityGuess.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/DataStructures/DensityMatrix.h>
#include <Utils/DataStructures/MatrixWithDerivatives.h>
#include <Utils/Math/DerivOrderEnum.h>
#include <Utils/Scf/MethodInterfaces/OverlapCalculator.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>

namespace Scine {
namespace Sparrow {

namespace nddo {

NDDODensityGuess::NDDODensityGuess(const Utils::ElementTypeCollection& elements, const ElementParameters& elementParameters,
                                   Utils::OverlapCalculator& overlapCalculator, const int& nElectrons, const int& nAOs)
  : elements_(elements),
    elementParameters_(elementParameters),
    overlapCalculator_(overlapCalculator),
    nElectrons_(nElectrons),
    nAOs_(nAOs) {
}

Utils::DensityMatrix NDDODensityGuess::calculateGuess() const {
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nAOs_, nAOs_);

  // stewart1990:
  // The guess is very crude: all off-diagonal matrix elements are set to zero, and
  // all on-diagonal terms on any atom are set equal to the core charge of that atom divided by the
  // number of atomic orbitals.

  // Alain:
  // Same as above, but, for off-diagonal elements,
  // use overlap matrix times factor for first guess of density matrix.
  // Division by two found by testing, seems to work better with it than without
  if (nAOs_ != 0) {
    overlapCalculator_.calculateOverlap(Utils::derivOrder::zero);
    P = overlapCalculator_.getOverlap().getMatrixXd() * nElectrons_ / (2 * nAOs_);
  }

  for (int i = 0; i < nAOs_; i++)
    for (int j = i + 1; j < nAOs_; j++)
      P(i, j) = P(j, i);

  // stewart 1990
  int index = 0;
  for (auto e : elements_) {
    double nEl = elementParameters_.get(e).coreCharge();
    auto nAOs = elementParameters_.get(e).nAOs();
    for (int j = 0; j < nAOs; j++) {
      P(index, index) = nEl / nAOs;
      index++;
    }
  }

  Utils::DensityMatrix d;
  d.setDensity(std::move(P), nElectrons_);
  return d;
}

void NDDODensityGuess::setNElectrons(int nElectrons) {
  nElectrons_ = nElectrons;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
