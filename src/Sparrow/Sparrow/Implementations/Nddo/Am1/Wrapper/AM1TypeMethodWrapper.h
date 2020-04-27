/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_AM1TYPEMETHODWRAPPER_H
#define SPARROW_AM1TYPEMETHODWRAPPER_H

/* Internal Includes */

#include <Sparrow/Implementations/Nddo/Am1/AM1Method.h>
#include <Sparrow/Implementations/Nddo/NDDOMethodWrapper.h>
/* External Includes */
#include <Utils/CalculatorBasics/PropertyList.h>
#include <string>

namespace Scine {
namespace Sparrow {

/**
 * @class AM1TypeMethodWrapper AM1TypeMethodWrapper.h
 * @brief A parent class for all AM1-type methods using the same formalism but different parameters (AM1, RM1, PM3).
 * This class uses the curiously recurrent template pattern to have compile-time resolution of the
 * model names, as well as the use of the right constructor. The constructors are method-dependant, having
 * independent classes would lead to unwanted code duplication.
 */
template<class AM1Type>
class AM1TypeMethodWrapper
  : public Utils::CloneInterface<Utils::Abstract<AM1TypeMethodWrapper<AM1Type>>, NDDOMethodWrapper> {
 public:
  /// @brief Constructor initializing the NDDOMethodWrapper with this instance of the class.
  AM1TypeMethodWrapper();
  /// @brief Default Destructor.
  ~AM1TypeMethodWrapper() override;
  /**
   * @brief Getter for the name of the method.
   * @returns Returns the name of the method.
   */
  std::string name() const final;
  /**
   * @brief Function to apply the settings to the underlying method
   */
  void applySettings() final;

 protected:
  bool successfulCalculation() const final;
  Eigen::MatrixXd getOneElectronMatrix() const final;
  Utils::SpinAdaptedMatrix getTwoElectronMatrix() const final;

  Utils::DensityMatrix getDensityMatrixGuess() const final;
  //! Initializes a method with the parameter file present in the settings.
  void initialize() override;
  Utils::LcaoMethod& getLcaoMethod() final;
  const Utils::LcaoMethod& getLcaoMethod() const final;

  void calculateImpl(Utils::derivativeType requiredDerivative) final;

  nddo::AM1Method method_;
};

class AM1MethodWrapper : public Utils::CloneInterface<AM1MethodWrapper, AM1TypeMethodWrapper<AM1MethodWrapper>> {
 public:
  static constexpr const char* model = "AM1";
  AM1MethodWrapper();
  // Rule of 5
  AM1MethodWrapper(const AM1MethodWrapper& rhs);
  AM1MethodWrapper& operator=(const AM1MethodWrapper& rhs);
  AM1MethodWrapper(AM1MethodWrapper&& rhs) = delete;
  AM1MethodWrapper& operator=(AM1MethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~AM1MethodWrapper() final;
};

class RM1MethodWrapper : public Utils::CloneInterface<RM1MethodWrapper, AM1TypeMethodWrapper<RM1MethodWrapper>> {
 public:
  static constexpr const char* model = "RM1";
  RM1MethodWrapper();
  // Rule of 5
  RM1MethodWrapper(const RM1MethodWrapper& rhs);
  RM1MethodWrapper& operator=(const RM1MethodWrapper& rhs);
  RM1MethodWrapper(RM1MethodWrapper&& rhs) = delete;
  RM1MethodWrapper& operator=(RM1MethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~RM1MethodWrapper() final;
};

class PM3MethodWrapper : public Utils::CloneInterface<PM3MethodWrapper, AM1TypeMethodWrapper<PM3MethodWrapper>> {
 public:
  static constexpr const char* model = "PM3";
  PM3MethodWrapper();
  // Rule of 5
  PM3MethodWrapper(const PM3MethodWrapper& rhs);
  PM3MethodWrapper& operator=(const PM3MethodWrapper& rhs);
  PM3MethodWrapper(PM3MethodWrapper&& rhs) = delete;
  PM3MethodWrapper& operator=(PM3MethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~PM3MethodWrapper() final;
};

} /* namespace Sparrow */
} /* namespace Scine */

#endif /* SPARROW_AM1TYPEMETHODWRAPPER_H */
