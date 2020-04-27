/**
 * @file DFTB2MethodWrapper.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_DFTB2METHODWRAPPER_H
#define SPARROW_DFTB2METHODWRAPPER_H

/* Internal Includes */
#include <Sparrow/Implementations/Dftb/DFTBMethodWrapper.h>
#include <Sparrow/Implementations/Dftb/Dftb2/DFTB2.h>
/* External Includes */
#include <string>

namespace Scine {
namespace Sparrow {
/**
 * @class DFTB2MethodWrapper DFTB2MethodWrapper.h
 * @brief A method wrapper handling DFTB2 calculations,
 * also known as SCC-DFTB, self-consistent charge DFTB.
 */
class DFTB2MethodWrapper : public Utils::CloneInterface<DFTB2MethodWrapper, DFTBMethodWrapper> {
 public:
  static constexpr const char* model = "DFTB2";

  /// @brief Default Constructor.
  DFTB2MethodWrapper();
  // Rule of 5
  DFTB2MethodWrapper(const DFTB2MethodWrapper& rhs);
  DFTB2MethodWrapper& operator=(const DFTB2MethodWrapper& rhs);
  DFTB2MethodWrapper(DFTB2MethodWrapper&& rhs) = delete;
  DFTB2MethodWrapper& operator=(DFTB2MethodWrapper&& rhs) = delete;
  /// @brief Default Destructor.
  ~DFTB2MethodWrapper() final;
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
  const Utils::LcaoMethod& getLcaoMethod() const final;
  Utils::LcaoMethod& getLcaoMethod() final;

  dftb::DFTB2 method_;
};

} /* namespace Sparrow */
} /* namespace Scine */

#endif /* SPARROW_DFTB2METHODWRAPPER_H */
