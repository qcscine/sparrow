/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_CISDATA_H
#define SPARROW_CISDATA_H

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/TwoCenterIntegralContainer.h>
#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterIntegralContainer.h>
#include <Sparrow/Implementations/TimeDependent/LinearResponseData.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
namespace Scine {
namespace Sparrow {

/**
 * @brief This class contains the infos needed to perform a CIS calculation.
 * This way excited states properties can be calculated in a LR-TD approach.
 * The data are given by value (i.e. copied) to allow the calculation to work also in an interactive setting.
 */
struct CISData : public LinearResponseData {
  const nddo::OneCenterIntegralContainer& oneCenterIntegrals;
  const nddo::TwoCenterIntegralContainer& twoCenterIntegrals;

  /**
   * @brief Constructor for an instance of the CISData struct.
   * @param oneCenterIntegralContainer
   * @param twoCenterIntegralContainer
   * @param MOs
   * @param orbitalEnergies
   * @param aoIndex
   * @param elements
   * @param occupation
   * @param overlapMatrix
   */
  CISData(const nddo::OneCenterIntegralContainer& oneCenterIntegralContainer,
          const nddo::TwoCenterIntegralContainer& twoCenterIntegralContainer, const Utils::MolecularOrbitals& MOs,
          const Utils::SingleParticleEnergies& orbitalEnergies, Utils::AtomsOrbitalsIndexes aoIndex,
          const Utils::ElementTypeCollection& elements, const Utils::LcaoUtils::ElectronicOccupation& occupation,
          const Eigen::MatrixXd& overlapMatrix)
    : LinearResponseData(MOs, orbitalEnergies, std::move(aoIndex), elements, occupation, overlapMatrix),
      oneCenterIntegrals(oneCenterIntegralContainer),
      twoCenterIntegrals(twoCenterIntegralContainer) {
  }

  template<class NDDOMethod>
  static CISData constructCISDataFromNDDOMethod(const NDDOMethod& method) {
    return CISData(method.getTwoElectronMatrix().getOneCenterIntegrals(),
                   method.getTwoElectronMatrix().getTwoCenterIntegrals(), method.getMolecularOrbitals(),
                   method.getSingleParticleEnergies(), method.getInitializer().getAtomsOrbitalsIndexes(),
                   method.getElementTypes(), method.getElectronicOccupation(), method.getOverlapMatrix());
  }
};

} // namespace Sparrow
} // namespace Scine
#endif // SPARROW_CISDATA_H
