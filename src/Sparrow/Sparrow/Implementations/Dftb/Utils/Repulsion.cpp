/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Repulsion.h"
#include "PairwiseRepulsion.h"
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/Math/AutomaticDifferentiation/MethodsHelpers.h>
#include <Utils/Typenames.h>

namespace Scine {
namespace Sparrow {

using namespace Utils::AutomaticDifferentiation;

namespace dftb {

Repulsion::Repulsion(const Utils::ElementTypeCollection& elements, const Utils::PositionCollection& positions,
                     const DFTBCommon::DiatomicParameterContainer& diatomicParameters)
  : RepulsionCalculator(elements, positions), nAtoms_(0), diatomicParameters_(diatomicParameters) {
}

Repulsion::~Repulsion() = default;

void Repulsion::initialize() {
  nAtoms_ = elements_.size();

  // Create 2D-vector of empty unique_ptr's
  pairRepulsions_ = Container(nAtoms_);
  for (int i = 0; i < nAtoms_; ++i)
    pairRepulsions_[i] = std::vector<PairRepulsion>(nAtoms_);
  for (int i = 0; i < nAtoms_; i++) {
    for (int j = i + 1; j < nAtoms_; j++) {
      initializePair(i, j);
    }
  }
}

void Repulsion::initializePair(int i, int j) {
  Utils::ElementType e1 = elements_[i];
  Utils::ElementType e2 = elements_[j];

  SKPair* parameters = diatomicParameters_[Utils::ElementInfo::Z(e1)][Utils::ElementInfo::Z(e2)].get();

  pairRepulsions_[i][j] = std::make_unique<dftb::PairwiseRepulsion>(parameters->getRepulsionParameters());
}

void Repulsion::calculateRepulsion(Utils::derivOrder order) {
  for (int i = 0; i < nAtoms_; i++) {
    for (int j = i + 1; j < nAtoms_; j++) {
      calculatePairRepulsion(i, j, order);
    }
  }
}

void Repulsion::calculatePairRepulsion(int i, int j, Utils::derivOrder order) {
  const auto& pA = positions_.row(i);
  const auto& pB = positions_.row(j);
  Eigen::Vector3d Rab = pB - pA;

  pairRepulsions_[i][j]->calculate(Rab, order);
}

double Repulsion::getRepulsionEnergy() const {
  double repulsion = 0;

#pragma omp parallel for reduction(+ : repulsion)
  for (int a = 0; a < nAtoms_; ++a) {
    for (int b = a + 1; b < nAtoms_; ++b) {
      auto p = pairRepulsions_[a][b]->getRepulsionEnergy();
      repulsion += p;
    }
  }

  return repulsion;
}

void Repulsion::addRepulsionDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::first>& derivatives) const {
  addRepulsionDerivativesImpl<Utils::derivativeType::first>(derivatives);
}

void Repulsion::addRepulsionDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_atomic>& derivatives) const {
  addRepulsionDerivativesImpl<Utils::derivativeType::second_atomic>(derivatives);
}

void Repulsion::addRepulsionDerivatives(
    Utils::AutomaticDifferentiation::DerivativeContainerType<Utils::derivativeType::second_full>& derivatives) const {
  addRepulsionDerivativesImpl<Utils::derivativeType::second_full>(derivatives);
}

template<Utils::derivativeType O>
void Repulsion::addRepulsionDerivativesImpl(Utils::AutomaticDifferentiation::DerivativeContainerType<O>& derivatives) const {
#pragma omp parallel for
  for (int a = 0; a < nAtoms_; ++a) {
    for (int b = a + 1; b < nAtoms_; b++) {
      auto dRep = pairRepulsions_[a][b]->getDerivative<O>();
#pragma omp critical
      { addDerivativeToContainer<O>(derivatives, a, b, dRep); }
    }
  }
}

} // namespace dftb
} // namespace Sparrow
} // namespace Scine
