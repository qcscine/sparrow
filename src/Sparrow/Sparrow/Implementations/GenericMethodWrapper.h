/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_GENERICMETHODWRAPPER_H
#define SPARROW_GENERICMETHODWRAPPER_H

/* External Includes */

#include <Core/Interfaces/Calculator.h>
#include <Core/Interfaces/WavefunctionOutputGenerator.h>
#include <Utils/CalculatorBasics.h>
#include <Utils/DataStructures/AtomicGtos.h>
#include <Utils/Settings.h>
#include <Utils/Technical/CloneInterface.h>
#include <string>

namespace Scine {
namespace Utils {
class LcaoMethod;
class DensityMatrix;
class DipoleMatrix;
enum class derivativeType;
} // namespace Utils
namespace Sparrow {

class DipoleMatrixCalculator;
class DipoleMomentCalculator;

/**
 * @class GenericMethodWrapper GenericMethodWrapper.h
 * @brief A MethodWrapper running Generic calculations.
 */
class GenericMethodWrapper : public Utils::CloneInterface<Utils::Abstract<GenericMethodWrapper>, Core::Calculator>,
                             public Core::WavefunctionOutputGenerator {
 public:
  /// @brief Default Constructor
  GenericMethodWrapper();
  /// @brief Default Destructor.
  ~GenericMethodWrapper() override;
  /**
   * @brief Sets new structure and initializes the underlying method with the parameter given in the settings.
   * @param structure The structure to be assigned.
   */
  void setStructure(const Utils::AtomCollection& structure) final;
  /**
   * @brief Getter for the underlying element types and positions.
   */
  std::unique_ptr<Utils::AtomCollection> getStructure() const final;
  /**
   * @brief Allows to modify the positions of the underlying Utils::AtomCollection
   * @param newPositions the new positions to be assigned to the underlying Utils::AtomCollection
   */
  void modifyPositions(Utils::PositionCollection newPositions) final;
  /**
   * @brief Getter for the coordinates of the underlying Utils::AtomCollection
   */
  const Utils::PositionCollection& getPositions() const final;
  /**
   * @brief Sets the properties to calculate.
   * @param requiredProperties a Utils::PropertyList object, contains an enum class that work as
   *        a bitset, switching on and off the bits representing a property.
   */
  void setRequiredProperties(const Utils::PropertyList& requiredProperties) final;
  /**
   * @brief Gets the properties to calculate.
   * @return requiredProperties a Utils::PropertyList object, contains an enum class that work as
   *        a bitset, switching on and off the bits representing a property.
   */
  Utils::PropertyList getRequiredProperties() const final;
  Utils::PropertyList possibleProperties() const override;
  /**
   * @brief The main function running calculations (dummy).
   *
   * @param dummy   A dummy parameter.
   * @return double Return the result of the calculation.
   */
  const Utils::Results& calculate(std::string description) final;
  /**
   * @brief Getter for the initial density matrix guess.
   * @return The initial density matrix guess.
   */
  virtual Utils::DensityMatrix getDensityMatrixGuess() const = 0;
  /**
   * @brief Accessor for the Settings used in this method wrapper.
   * @returns Utils::Settings& The Settings.
   */
  Utils::Settings& settings() final;
  /**
   * @brief Const accessor for the Settings used in this method wrapper.
   * @returns const Utils::Settings& The Settings.
   */
  const Utils::Settings& settings() const final;
  /**
   * @brief Accessor for the Results stored in this method wrapper.
   * @returns Utils::Results& The results of the previous calculation.
   */
  Utils::Results& results() final;
  /**
   * @brief Const accessor for the Results used in this method wrapper.
   * @returns const Utils::Results& The results of the previous calculation.
   */
  const Utils::Results& results() const final;

  /**
   * @brief Whether the calculator supports a method family
   * In this library, every calculator supports just one "method family": itself.
   * @param methodFamily identifier for the method family
   * @return whether the calculator supports a method family
   */
  bool supportsMethodFamily(const std::string& methodFamily) const final;

  //! Get the underlying LcaoMethod
  virtual Utils::LcaoMethod& getLcaoMethod() = 0;
  virtual const Utils::LcaoMethod& getLcaoMethod() const = 0;
  /// @brief Returns a Sto-6G expansion
  std::string getStoNGExpansionPath() const;
  void generateWavefunctionInformation(const std::string& filename) final;
  void generateWavefunctionInformation(std::ostream& out) final;

  void loadState(std::shared_ptr<Core::State> state) override;
  std::shared_ptr<Core::State> getState() const final;

 protected:
  std::unique_ptr<Utils::Settings> settings_;
  Utils::Results results_;
  //! Initializes a method with the parameter file present in the settings.
  virtual void initialize() = 0;
  //! Determines the highest derivative type needed based on the property needed for calculation.
  Utils::derivativeType highestDerivativeRequired() const;
  virtual void assembleResults(const std::string& description);

  virtual bool successfulCalculation() const = 0;
  virtual bool canCalculateAnalyticalHessian() const;

  /**
   * @brief Function to apply the settings to the actual calculation method.
   * This method is needed in every calculator as the modification of the
   * settings object does not cause a direct propagation to the underlying
   * calculation method. It is not part of the public interface as an external
   * user should not have to decide when/where to use this. In the implementation
   * of a calculator care must be taken that at calculation time the settings
   * have correctly been propagated to the underlying calculator method.
   */
  virtual void applySettings();
  virtual bool getZPVEInclusion() const = 0;

  //! Method-dependent implementation of the calculate member function
  virtual void calculateImpl(Utils::derivativeType requiredDerivative) = 0;

  std::unique_ptr<DipoleMomentCalculator> dipoleCalculator_;
  std::unique_ptr<DipoleMatrixCalculator> dipoleMatrixCalculator_;
  Utils::PropertyList requiredProperties_;
};

} /* namespace Sparrow */
} /* namespace Scine */

#endif /* SPARROW_GENERICMETHODWRAPPER_H */
