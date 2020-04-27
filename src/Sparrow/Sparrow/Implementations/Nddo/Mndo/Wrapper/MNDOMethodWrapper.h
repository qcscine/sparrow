/**
 * @file MNDOMethodWrapper.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_MNDOMETHODWRAPPER_H
#define SPARROW_MNDOMETHODWRAPPER_H

/* Internal Includes */

#include <Sparrow/Implementations/Nddo/Mndo/MNDOMethod.h>
#include <Sparrow/Implementations/Nddo/NDDOMethodWrapper.h>
/* External Includes */
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/Technical/CloneInterface.h>
#include <string>

namespace Scine {
namespace Sparrow {
/**
 * @class MNDOMethodWrapper MNDOMethodWrapper.h
 * @brief A method wrapper handling MNDO calculations.
 */
class MNDOMethodWrapper : public Utils::CloneInterface<MNDOMethodWrapper, NDDOMethodWrapper> {
 public:
  static constexpr const char* model = "MNDO";

  /// @brief Default Constructor.
  MNDOMethodWrapper();
  // Rule of 5
  MNDOMethodWrapper(const MNDOMethodWrapper& rhs);
  MNDOMethodWrapper& operator=(const MNDOMethodWrapper& rhs);
  MNDOMethodWrapper(MNDOMethodWrapper&& rhs) = delete;
  MNDOMethodWrapper& operator=(MNDOMethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~MNDOMethodWrapper() final;
  /**
   * @brief Getter for the name of the underlying method.
   * @returns Returns the name of the underlying method.
   */
  std::string name() const final;
  /**
   * @brief Function to apply the settings to the underlying method.
   */
  void applySettings() final;

 private:
  bool successfulCalculation() const final;
  Eigen::MatrixXd getOneElectronMatrix() const final;
  Utils::SpinAdaptedMatrix getTwoElectronMatrix() const final;

  Utils::DensityMatrix getDensityMatrixGuess() const final;

  //! Initializes a method with the parameter file present in the settings.
  void initialize() final;

  Utils::LcaoMethod& getLcaoMethod() final;
  const Utils::LcaoMethod& getLcaoMethod() const final;

  void calculateImpl(Utils::derivativeType requiredDerivative) final;

  nddo::MNDOMethod method_;
};

} /* namespace Sparrow */
} /* namespace Scine */

#endif /* SPARROW_MNDOMETHODWRAPPER_H */
