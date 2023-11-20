/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_TDDFTBDATA_H
#define SPARROW_TDDFTBDATA_H

#include <Sparrow/Implementations/TimeDependent/LinearResponseData.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Typenames.h>
#include <memory>

namespace Scine {
namespace Sparrow {

/**
 * @brief This class contains the infos needed to perform a TD-DFTB calculation.
 * This way excited states properties can be calculated in a LR-TD approach.
 * The data are given by reference to prevent useless copies.
 * Right now the data needed to perform TD-SCC-DFTB are present. In order to implement
 * TD-DFTB3 other data might be needed.
 */
class TDDFTBData : public LinearResponseData {
 public:
  /// @brief Gamma parameters size: nAtoms x nAtoms
  std::shared_ptr<Eigen::MatrixXd> gammaMatrix;
  /// @brief Magnetic Hubbard parameters (spin constants) size: nAtoms
  std::shared_ptr<Eigen::VectorXd> spinConstants;

  TDDFTBData(const Utils::MolecularOrbitals& MOs, const Utils::SingleParticleEnergies& orbitalEnergies,
             Utils::AtomsOrbitalsIndexes aoIndex, const Utils::ElementTypeCollection& elements,
             const Utils::LcaoUtils::ElectronicOccupation& occupation, const Eigen::MatrixXd& overlapMatrix,
             const Eigen::MatrixXd& gMatrix, std::shared_ptr<Eigen::VectorXd> spinConstantVector)
    : LinearResponseData(MOs, orbitalEnergies, std::move(aoIndex), elements, occupation, overlapMatrix) {
    gammaMatrix = std::make_shared<Eigen::MatrixXd>(gMatrix);
    spinConstants = std::move(spinConstantVector);
  }

  template<class DFTBMethod>
  static TDDFTBData constructTDDFTBDataFromDFTBMethod(const DFTBMethod& method) {
    return TDDFTBData(method.getMolecularOrbitals(), method.getSingleParticleEnergies(),
                      method.getInitializer()->getAtomsOrbitalsIndexes(), method.getElementTypes(),
                      method.getElectronicOccupation(), method.getOverlapMatrix(), method.calculateGammaMatrix(),
                      method.calculateSpinConstantVector());
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_TDDFTBDATA_H
