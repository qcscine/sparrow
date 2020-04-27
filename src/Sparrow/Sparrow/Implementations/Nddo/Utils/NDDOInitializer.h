/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_NDDOINITIALIZER_H
#define SPARROW_NDDOINITIALIZER_H

#include <Sparrow/Implementations/Nddo/Utils/IntegralsEvaluationUtils/oneCenterIntegralContainer.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementPairParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/ElementParameters.h>
#include <Sparrow/Implementations/Nddo/Utils/ParameterUtils/RawParametersContainer.h>
#include <Utils/DataStructures/AtomsOrbitalsIndexes.h>
#include <Utils/Scf/MethodInterfaces/StructureDependentInitializer.h>
#include <string>
#include <vector>

namespace Scine {
namespace Sparrow {

namespace nddo {

/**
 * @brief Settings for generic NDDO methods.
 * Reads the parameters and applies them to the system of interest.
 * Depending on the method, the basis functions used can either be BasisFunction::sp (MNDO, AM1) of BasisFunction::spd
 * (i.e. PM6, AM1*, MNDO/d) and it can have just atomic parameters (i.e. AM1, MNDO) or also diatomic parameters (i.e.
 * PM6)
 */
class NDDOInitializer : public Utils::StructureDependentInitializer {
 public:
  explicit NDDOInitializer(BasisFunctions basisFunctions = BasisFunctions::spd, bool hasDiatomicParameters = true)
    : basisFunctions_(basisFunctions), hasDiatomicParameters_(hasDiatomicParameters){};

  /**
   * @brief (Re)generate values and run-time parameters from the current raw parameters.
   *         Only needed if the parameters are modified manually.
   * @param elements a vector containing the elements constituting the molecule.
   * @param basisFunctions Whether the method just accomodate s and p basis functions (i.e. AM1, MNDO), or if it can
   *        also activate d basis functions.
   * @param hasDiatomicParameters Whether the method also has diatomic parameters (i.e. PM6).
   */
  void applyRawParameters(const Utils::ElementTypeCollection& elements);
  /*! Load the parameters from a file. */
  void readParameters(const std::string& parameterPath);
  /*! Save the parameters to a file. */
  void saveParameters(const std::string& fileName);
  /*! Initialize the method <b>after</b> the parameters have been set or loaded. */
  void initialize(const Utils::ElementTypeCollection& elements) override;

  Utils::AtomsOrbitalsIndexes getAtomsOrbitalsIndexes() const override;
  unsigned getNumberElectronsForUnchargedSpecies() const override;
  std::vector<double> getCoreCharges() const override;
  bool unrestrictedCalculationPossible() const override;

  /*! Get reference to the class for raw NDDO parameters. */
  RawParametersContainer& getRawParameters();
  /*! Get const reference to the class for raw NDDO parameters. */
  const RawParametersContainer& getRawParameters() const;

  const ElementParameters& getElementParameters();
  const ElementPairParameters& getElementPairParameters();
  const OneCenterIntegralContainer& getOneCenterIntegrals();

  BasisFunctions getBasisFunctions() const;

 private:
  ElementParameters elementParameters_;
  ElementPairParameters elementPairParameters_;
  OneCenterIntegralContainer oneCenterIntegrals_;
  RawParametersContainer rawParameters_;

  unsigned int nElectronsForUnchargedSpecies_ = 0;
  std::vector<double> coreCharges_;
  Utils::AtomsOrbitalsIndexes aoIndexes_;

  BasisFunctions basisFunctions_;
  bool hasDiatomicParameters_;
};

} // namespace nddo
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_NDDOINITIALIZER_H
