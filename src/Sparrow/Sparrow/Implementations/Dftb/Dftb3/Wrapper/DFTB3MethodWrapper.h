/**
 * @file DFTB3MethodWrapper.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_DFTB3METHODWRAPPER_H
#define SPARROW_DFTB3METHODWRAPPER_H

/* Internal Includes */
#include <Sparrow/Implementations/Dftb/DFTBMethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb3/DFTB3.h>
/* External Includes */
#include <Utils/Technical/CloneInterface.h>
#include <string>

namespace Scine {
namespace Sparrow {
/**
 * @class DFTB3MethodWrapper DFTB3MethodWrapper.h
 * @brief A method wrapper handling DFTB3 calculations.
 */
class DFTB3MethodWrapper : public Utils::CloneInterface<DFTB3MethodWrapper, DFTBMethodWrapper> {
 public:
  static constexpr const char* model = "DFTB3";

  /// @brief Default Constructor.
  DFTB3MethodWrapper();
  // Rule of 5
  DFTB3MethodWrapper(const DFTB3MethodWrapper& rhs);
  DFTB3MethodWrapper& operator=(const DFTB3MethodWrapper& rhs);
  DFTB3MethodWrapper(DFTB3MethodWrapper&& rhs) = delete;
  DFTB3MethodWrapper& operator=(DFTB3MethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~DFTB3MethodWrapper() final;
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
  Utils::DensityMatrix getDensityMatrixGuess() const final;
  //! Initializes a method with the parameter file present in the settings.
  void initialize() final;
  //! Calls the underlying method's calculate() function.
  void calculateImpl(Utils::derivativeType requiredDerivative) final;
  //! Get the underlying method as a LCAO method.
  Utils::LcaoMethod& getLcaoMethod() final;
  const Utils::LcaoMethod& getLcaoMethod() const final;
  dftb::DFTB3 method_;
};

} /* namespace Sparrow */
} /* namespace Scine */

#endif /* SPARROW_DFTB3METHODWRAPPER_H */
