/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
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

class DipoleMatrixCalculator;

/**
 * @class NDDOMethodWrapper
 * @brief Abstract class acting as a generic Wrapper for NDDO methods.
 */
class NDDOMethodWrapper : public Utils::CloneInterface<Utils::Abstract<NDDOMethodWrapper>, GenericMethodWrapper> {
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

 protected:
  // Extracted method from all copy constructors and copy assignment operators.
  template<class NDDOMethod>
  void copyInto(NDDOMethod& instance, const NDDOMethod& classToCopy) {
    instance.results() = classToCopy.results();
    instance.settings() = classToCopy.settings();
    // Concurrent calling of the logger introduces race conditions
    // that eventually trigger a segfault
    instance.setStructure(*classToCopy.getStructure());
    instance.loadState(classToCopy.getState());
  }
  virtual Eigen::MatrixXd getOneElectronMatrix() const = 0;
  virtual Utils::SpinAdaptedMatrix getTwoElectronMatrix() const = 0;
  void assembleResults(const std::string& description) final;

  void applySettings(std::unique_ptr<Utils::Settings>& settings, Utils::ScfMethod& method);
  bool getZPVEInclusion() const final;
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_NDDOMETHODWRAPPER_H
