/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_MOLDENFILEGENERATOR_H
#define SPARROW_MOLDENFILEGENERATOR_H

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/GeneralTypes.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Eigen/Core>
#include <map>
#include <memory>
#include <ostream>
namespace Scine {
namespace Sparrow {
class GenericMethodWrapper;
/**
 * @class MoldenFileGenerator @file MoldenFileGenerator.h
 * @brief Class to create the wavefunction information needed for outputting densities,...
 * Note that NDDO methods have their own STO-6G expansion, fine tuned according to their parameters,
 * while for DFTB methods the STO-6G expansion of the PM6 method was used.
 * Since STO-6G expansions are very similar one-another, the implementation ease of this was
 * deemed a satisfactory compromise.
 * If the underlying calculator has not been initialized, then care must be taken that outside of this
 * function a calculation is performed.
 * This class handles only s, p and d orbitals.
 */
class MoldenFileGenerator {
 public:
  /**
   * @brief Sets the calculator from which to generate a molden file.
   * @param calculator a reference to a GenericMethodWrapper instance.
   */
  explicit MoldenFileGenerator(const GenericMethodWrapper& calculator);
  ~MoldenFileGenerator() = default;

  /**
   * @brief Method to generate the molden input and print it somewhere.
   * @param out The stream the molden input is written to.
   */
  void generateWavefunctionInformation(std::ostream& out) const;

 private:
  virtual std::vector<Utils::AtomicGtos> getStoNGExpansion() const;
  void generateAtomBlock(std::ostream& out) const;
  void generateGTOBlock(std::ostream& out) const;
  void generateMolecularOrbitalsBlock(std::ostream& out) const;
  void writeMOBlock(std::ostream& out, Eigen::MatrixXd moMatrix, const std::vector<int>& filledOrb,
                    const std::vector<double>& moEnergies, const std::string& spin) const;
  const GenericMethodWrapper& calculator_;
  static std::map<int, nddo::GeneralTypes::orb_t> indexMap_;
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_MOLDENFILEGENERATOR_H
