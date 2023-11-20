/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_LINEARRESPONSEDATA_H
#define SPARROW_LINEARRESPONSEDATA_H

#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/Scf/LcaoUtils/ElectronicOccupation.h>
#include <Utils/Typenames.h>
#include <map>
#include <memory>
#include <vector>

namespace Scine {
namespace Sparrow {

/**
 * @brief This class contains the infos needed to perform a linear response calculation.
 * This way excited states properties can be calculated in a LR-TD approach.
 * It serves as a common base class for  TD-DFTB and CIS data.
 * The data are given by reference to prevent useless copies.
 */
struct LinearResponseData {
  const Utils::MolecularOrbitals& molecularOrbitals;
  const Utils::SingleParticleEnergies& MOEnergies;
  const Utils::AtomsOrbitalsIndexes AOInfo;
  const Utils::ElementTypeCollection& elements;
  const Eigen::MatrixXd& overlapMatrix;
  const Utils::LcaoUtils::ElectronicOccupation& occupation;

  LinearResponseData(const Utils::MolecularOrbitals& MOs, const Utils::SingleParticleEnergies& orbitalEnergies,
                     Utils::AtomsOrbitalsIndexes aoIndex, const Utils::ElementTypeCollection& elements,
                     const Utils::LcaoUtils::ElectronicOccupation& occupation, const Eigen::MatrixXd& overlapMatrix)
    : molecularOrbitals(MOs),
      MOEnergies(orbitalEnergies),
      AOInfo(std::move(aoIndex)),
      elements(elements),
      overlapMatrix(overlapMatrix),
      occupation(occupation) {
  }
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_LINEARRESPONSEDATA_H
