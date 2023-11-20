/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "CISPseudoDensityBuilder.h"

namespace Scine {
namespace Sparrow {

template<>
CISPseudoDensityBuilder<Utils::Reference::Restricted>::CISPseudoDensityBuilder(
    const Utils::MolecularOrbitals& molecularOrbitals, const Utils::LcaoUtils::ElectronicOccupation& occupation) {
  const Eigen::MatrixXd& allRestricted = molecularOrbitals.restrictedMatrix();
  int nRestrictedMOs = static_cast<int>(allRestricted.cols());
  const std::vector<int>& filledOrbitals = occupation.getFilledRestrictedOrbitals();

  Eigen::MatrixXd virtualOrbitals(allRestricted.rows(), nRestrictedMOs - filledOrbitals.size());
  Eigen::MatrixXd occupiedOrbitals(allRestricted.rows(), filledOrbitals.size());

  constructOrbitals(occupiedOrbitals, virtualOrbitals, allRestricted, filledOrbitals);

  virtualMolecularOrbitals_.restricted = virtualOrbitals;
  occupiedMolecularOrbitals_.restricted = occupiedOrbitals;
}

template<>
CISPseudoDensityBuilder<Utils::Reference::Unrestricted>::CISPseudoDensityBuilder(
    const Utils::MolecularOrbitals& molecularOrbitals, const Utils::LcaoUtils::ElectronicOccupation& occupation) {
  const Eigen::MatrixXd& allAlpha = molecularOrbitals.alphaMatrix();
  const Eigen::MatrixXd& allBeta = molecularOrbitals.betaMatrix();

  auto nAlphaMOs = static_cast<int>(allAlpha.cols());
  auto nBetaMOs = static_cast<int>(allBeta.cols());

  const auto& filledAlphaOrbitals = occupation.getFilledAlphaOrbitals();
  const auto& filledBetaOrbitals = occupation.getFilledBetaOrbitals();

  Eigen::MatrixXd virtualAlphaOrbitals(allAlpha.rows(), nAlphaMOs - filledAlphaOrbitals.size());
  Eigen::MatrixXd virtualBetaOrbitals(allBeta.rows(), nBetaMOs - filledBetaOrbitals.size());
  Eigen::MatrixXd occupiedAlphaOrbitals(allAlpha.rows(), filledAlphaOrbitals.size());
  Eigen::MatrixXd occupiedBetaOrbitals(allBeta.rows(), filledBetaOrbitals.size());

  constructOrbitals(occupiedAlphaOrbitals, virtualAlphaOrbitals, allAlpha, filledAlphaOrbitals);
  constructOrbitals(occupiedBetaOrbitals, virtualBetaOrbitals, allBeta, filledBetaOrbitals);

  virtualMolecularOrbitals_.alpha = virtualAlphaOrbitals;
  virtualMolecularOrbitals_.beta = virtualBetaOrbitals;
  occupiedMolecularOrbitals_.alpha = occupiedAlphaOrbitals;
  occupiedMolecularOrbitals_.beta = occupiedBetaOrbitals;
}

template<Utils::Reference restrictedness>
void CISPseudoDensityBuilder<restrictedness>::constructOrbitals(Eigen::MatrixXd& occupiedOrbitals,
                                                                Eigen::MatrixXd& virtualOrbitals, const Eigen::MatrixXd& allMOs,
                                                                const std::vector<int>& filledOrbitals) {
  auto nMOs = allMOs.cols();
  Eigen::Index iterEmpty = 0;
  Eigen::Index iterFilled = 0;

  for (Eigen::Index i = 0; i < nMOs; ++i) {
    if (iterFilled < static_cast<Eigen::Index>(filledOrbitals.size()) && i == filledOrbitals[iterFilled]) {
      occupiedOrbitals.col(iterFilled) = allMOs.col(i);
      ++iterFilled;
    }
    else {
      virtualOrbitals.col(iterEmpty) = allMOs.col(i);
      ++iterEmpty;
    }
  }
}

template<>
const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd>
CISPseudoDensityBuilder<Utils::Reference::Unrestricted>::getPseudoDensityMatrix(
    const Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::VectorXd>& guessVector) const {
  Utils::SpinAdaptedContainer<Utils::Reference::Unrestricted, Eigen::MatrixXd> pseudoDensity;

  pseudoDensity.alpha = mapAndMultiply(guessVector.alpha, virtualMolecularOrbitals_.alpha, occupiedMolecularOrbitals_.alpha);
  pseudoDensity.beta = mapAndMultiply(guessVector.beta, virtualMolecularOrbitals_.beta, occupiedMolecularOrbitals_.beta);

  return pseudoDensity;
}

template<>
const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd>
CISPseudoDensityBuilder<Utils::Reference::Restricted>::getPseudoDensityMatrix(
    const Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::VectorXd>& guessVector) const {
  Utils::SpinAdaptedContainer<Utils::Reference::Restricted, Eigen::MatrixXd> pseudoDensity;
  pseudoDensity.restricted =
      mapAndMultiply(guessVector.restricted, virtualMolecularOrbitals_.restricted, occupiedMolecularOrbitals_.restricted);

  return pseudoDensity;
}

template<Utils::Reference restrictedness>
Eigen::MatrixXd CISPseudoDensityBuilder<restrictedness>::mapAndMultiply(const Eigen::VectorXd& Vector,
                                                                        const Eigen::MatrixXd& virtualMOs,
                                                                        const Eigen::MatrixXd& occupiedMOs) const {
  Eigen::Map<const Eigen::MatrixXd> mappedVector(Vector.data(), virtualMOs.cols(), occupiedMOs.cols());

  return occupiedMOs * mappedVector.transpose() * virtualMOs.transpose();
}

} // namespace Sparrow
} // namespace Scine
