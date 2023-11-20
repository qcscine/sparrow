/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_NDDOMETHODWRAPPER_H
#define SPARROW_NDDOMETHODWRAPPER_H

#include <Sparrow/Implementations/GenericMethodWrapper.h>
#include <Utils/Technical/CloneInterface.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace Scine {

namespace Utils {
class Settings;
class ScfMethod;
class SpinAdaptedMatrix;
} // namespace Utils

namespace Sparrow {
class CISData;
class DipoleMatrixCalculator;

/**
 * @class NDDOMethodWrapper
 * @brief Abstract class acting as a genericWrapper for NDDO methods.
 */
class NDDOMethodWrapper
  : public Utils::CloneInterface<Utils::Abstract<NDDOMethodWrapper>, GenericMethodWrapper, Core::Calculator> {
 public:
  /**
   * @brief Constructor.
   * It plays with the befriended states handler by giving it a *this reference.
   */
  NDDOMethodWrapper();
  ~NDDOMethodWrapper() override;
  /**
   * @brief Get the list of the possible properties to calculate analytically.
   * Since it is the same for all NDDO, have it stated here.
   * If they diverge, override this in each method wrapper.
   */
  Utils::PropertyList possibleProperties() const override;

  friend class NDDOStatesHandler;

  /**
   * @brief This function is needed in the calcualtion of the CIS matrix in
   *        linear response method.
   */
  CISData getCISData() const;

  /**
   * @brief Whether the calculator has no underlying Python code and can therefore
   * release the global interpreter lock in Python bindings
   */
  bool allowsPythonGILRelease() const override {
    return true;
  };

 protected:
  // Extracted method from all copy constructors and copy assignment operators.
  template<class NDDOMethod>
  void copyInto(NDDOMethod& instance, const NDDOMethod& classToCopy) {
    auto results = classToCopy.results();
    instance.settings() = classToCopy.settings();
    // Concurrent calling of the logger introduces race conditions
    // that eventually trigger a segfault
    instance.setStructure(*classToCopy.getStructure());
    instance.results() = std::move(results);
    instance.loadState(classToCopy.getState());
    instance.setLog(classToCopy.getLog());
  }
  virtual Eigen::MatrixXd getOneElectronMatrix() const = 0;
  virtual Utils::SpinAdaptedMatrix getTwoElectronMatrix() const = 0;
  virtual CISData getCISDataImpl() const = 0;
  void assembleResults(const std::string& description) final;

  void applySettings(std::unique_ptr<Utils::Settings>& settings, Utils::ScfMethod& method);
  bool getZPVEInclusion() const final;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_NDDOMETHODWRAPPER_H
