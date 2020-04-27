/**
 * @file PM6MethodWrapper.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_PM6METHODWRAPPER_H
#define SPARROW_PM6METHODWRAPPER_H

/* Internal Includes */

#include <Sparrow/Implementations/Nddo/NDDOMethodWrapper.h>
#include <Sparrow/Implementations/Nddo/Pm6/PM6Method.h>
/* External Includes */
#include <Utils/Technical/CloneInterface.h>
#include <string>

namespace Scine {
namespace Sparrow {
/**
 * @class PM6MethodWrapper PM6MethodWrapper.h
 * @brief A method wrapper running PM6 calculations.
 */
class PM6MethodWrapper : public Utils::CloneInterface<PM6MethodWrapper, NDDOMethodWrapper> {
 public:
  static constexpr const char* model = "PM6";

  /// @brief Default Constructor.
  PM6MethodWrapper();
  // Rule of 5
  PM6MethodWrapper(const PM6MethodWrapper& rhs);
  PM6MethodWrapper& operator=(const PM6MethodWrapper& rhs);
  PM6MethodWrapper(PM6MethodWrapper&& rhs) = delete;
  PM6MethodWrapper& operator=(PM6MethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~PM6MethodWrapper() final;
  /**
   * @brief Getter for the name of the method.
   * @returns Returns the name of the method.
   */
  std::string name() const final;
  /**
   * @brief Function to apply the settings to the underlying method
   */
  void applySettings() final;

 private:
  bool successfulCalculation() const final;
  Eigen::MatrixXd getOneElectronMatrix() const final;
  Utils::SpinAdaptedMatrix getTwoElectronMatrix() const final;

  Utils::DensityMatrix getDensityMatrixGuess() const final;
  //! Initializes a method with the parameter file present in the settings.
  void initialize() final;

  //! Returns the underlying method
  Utils::LcaoMethod& getLcaoMethod() final;
  const Utils::LcaoMethod& getLcaoMethod() const final;
  void calculateImpl(Utils::derivativeType requiredDerivative) final;
  nddo::PM6Method method_;
};

} /* namespace Sparrow */
} /* namespace Scine */

#endif /* SPARROW_PM6METHODWRAPPER_H */
