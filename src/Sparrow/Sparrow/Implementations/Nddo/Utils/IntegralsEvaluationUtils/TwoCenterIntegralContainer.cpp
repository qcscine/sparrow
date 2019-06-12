/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TwoCenterIntegralContainer.h"
#include "Global2c2eMatrix.h"
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Utils/Math/DerivOrderEnum.h>

namespace Scine {
namespace Sparrow {

namespace nddo {

TwoCenterIntegralContainer::TwoCenterIntegralContainer(const Utils::ElementTypeCollection& elements,
                                                       const Utils::PositionCollection& positions, const ElementParameters& ep)
  : elementParameters_(ep), elementTypes_(elements), positions_(positions) {
}

void TwoCenterIntegralContainer::initialize() {
  nAtoms_ = static_cast<unsigned int>(elementTypes_.size());

  matrices_.clear();
  matrices_.resize(nAtoms_, std::vector<integralMatrix_t>(nAtoms_));

  for (unsigned i = 0; i < nAtoms_; ++i) {
    for (unsigned int j = i + 1; j < nAtoms_; j++) {
      initializePair(i, j);
    }
  }
}

void TwoCenterIntegralContainer::initializePair(unsigned int i, unsigned int j) {
  Utils::ElementType e1 = elementTypes_[i];
  Utils::ElementType e2 = elementTypes_[j];
  const auto& p1 = elementParameters_.get(e1);
  const auto& p2 = elementParameters_.get(e2);
  int nAOs1 = p1.nAOs();
  int nAOs2 = p2.nAOs();
  unsigned int l1 = (nAOs1 == 1) ? 0 : (nAOs1 == 4) ? 1 : 2;
  unsigned int l2 = (nAOs2 == 1) ? 0 : (nAOs2 == 4) ? 1 : 2;

  matrices_[i][j] = std::make_shared<multipole::Global2c2eMatrix>(l1, l2, p1.chargeSeparations(), p2.chargeSeparations(),
                                                                  p1.klopmanParameters(), p2.klopmanParameters());
  if (e1 == e2)
    matrices_[i][j]->setSymmetric(true);
}

void TwoCenterIntegralContainer::update(Utils::derivOrder order) {
#pragma omp parallel for schedule(dynamic)
  for (unsigned i = 0; i < nAtoms_; ++i) {
    for (unsigned int j = i + 1; j < nAtoms_; j++) {
      updatePair(i, j, order);
    }
  }
}

void TwoCenterIntegralContainer::updatePair(unsigned int i, unsigned int j, Utils::derivOrder order) {
  auto pA = positions_.row(i);
  auto pB = positions_.row(j);
  Eigen::RowVector3d Rab = pB - pA;

  if (order == Utils::derivOrder::zero)
    matrices_[i][j]->calculate<Utils::derivOrder::zero>(Rab);
  else if (order == Utils::derivOrder::one)
    matrices_[i][j]->calculate<Utils::derivOrder::one>(Rab);
  else if (order == Utils::derivOrder::two)
    matrices_[i][j]->calculate<Utils::derivOrder::two>(Rab);
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine
