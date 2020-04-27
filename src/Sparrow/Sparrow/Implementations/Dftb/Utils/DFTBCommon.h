/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTBCOMMON_H
#define SPARROW_DFTBCOMMON_H

#include "SKAtom.h"
#include "SKPair.h"
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Scf/MethodInterfaces/StructureDependentInitializer.h>
#include <Utils/Typenames.h>
#include <Eigen/Core>
#include <array>
#include <memory>
#include <vector>

namespace Scine {
namespace Utils {
class MatrixWithDerivatives;
class SingleParticleEnergies;
class AtomsOrbitalsIndexes;
class DensityMatrix;
} // namespace Utils
namespace Sparrow {

namespace dftb {

class DFTBCommon : public Utils::StructureDependentInitializer {
 public:
  using AtomicParameterContainer = std::vector<std::unique_ptr<SKAtom>>;
  using SKPairPtr = std::unique_ptr<SKPair>;
  using SKPairLine = std::array<SKPairPtr, 110>;
  using DiatomicParameterContainer = std::array<SKPairLine, 110>;

  DFTBCommon(const Utils::ElementTypeCollection& elements, int& nEl, int& charge, AtomicParameterContainer& atomicPar,
             DiatomicParameterContainer& diatomicPar);
  ~DFTBCommon() override;
  void initialize(const std::string& path, unsigned dftbType); // TODO: delete.
  void setMethodDetails(const std::string& path, unsigned dftbType);
  void initialize(const Utils::ElementTypeCollection& elementTypes) override;

  void reinitializeParameters();

  unsigned int getnAOs() const {
    return nAOs;
  }
  bool hasHubbardDerivatives() const {
    return DFTB3Possible;
  }
  bool unrestrictedCalculationPossible() const override;
  const DiatomicParameterContainer& getPairParameters() const {
    return pairParameters;
  }
  std::vector<double> getCoreCharges() const override {
    return coreCharges_;
  }
  Utils::AtomsOrbitalsIndexes getAtomsOrbitalsIndexes() const override {
    return aoIndexes_;
  }
  unsigned getNumberElectronsForUnchargedSpecies() const override {
    return nInitialElectrons_;
  }

 private:
  void readSpinConstants(const std::string& path);
  void readHubbardDerivatives(const std::string& path);

  const int nElements_;
  std::vector<bool> atomTypePresent; // tells if atomType present in the structure
  std::vector<std::string> listAtomTypes;
  AtomicParameterContainer& atomParameters;   // parameters for atoms
  DiatomicParameterContainer& pairParameters; // List of pointers to parameters

  unsigned int nAOs;               // Total number of Atomic orbitals
  unsigned int nAtoms;             // Total number of atoms
  int& nElectrons_;                // Total number of electrons
  unsigned int nInitialElectrons_; // initial number of electrons
  int& molecularCharge_;
  double noInteractionEnergy; // sum of energies of isolated atoms
  bool spinPolarizedPossible; // tells if parameters for SDFTB are present
  bool DFTB3Possible;         // tells if parameters for DFTB3 are present
  std::vector<double> coreCharges_;
  Utils::AtomsOrbitalsIndexes aoIndexes_;
  const Utils::ElementTypeCollection& elementTypes_;

  std::string path_;
  unsigned dftbType_;
};

} // namespace dftb

} // namespace Sparrow
} // namespace Scine
#endif // DFTBCOMMON_H
