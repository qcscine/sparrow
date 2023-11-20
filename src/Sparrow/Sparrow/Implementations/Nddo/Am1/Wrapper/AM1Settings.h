/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_AM1SETTINGS_H
#define SPARROW_AM1SETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <utility>
namespace Scine {
namespace Sparrow {

/**
 * @class AM1TypeSettings AM1TypeSettings.h
 * @brief The Settings specific to all the AM1-type methods, i.e. AM1, RM1, PM3.
 * This class uses the curiously recurrent template pattern to have compile-time resolution of the
 * parameters position and settings descriptions, which both are method-dependant. Having independent classes
 * would lead to unwanted code duplication.
 */
class AM1TypeSettings : public Scine::Utils::Settings {
 public:
  //! Take the names from the derived class and use them
  AM1TypeSettings(std::string settingsDescription) : Settings(std::move(settingsDescription)) {
    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "");
    Utils::UniversalSettings::BoolDescriptor useNDDODipoleApprox("Sets use of NDDO dipole approximation.");
    useNDDODipoleApprox.setDefaultValue(true);
    _fields.push_back(Utils::SettingsNames::NDDODipoleApproximation, std::move(useNDDODipoleApprox));

    resetToDefaults();
  }
};

class AM1Settings : public AM1TypeSettings {
 public:
  static constexpr const char* settingsDescription = "AM1Settings";
  AM1Settings() : AM1TypeSettings(settingsDescription) {
    // Method
    Utils::UniversalSettings::StringDescriptor method("The method to be used.");
    method.setDefaultValue("am1");
    _fields.push_back(Utils::SettingsNames::method, method);

    resetToDefaults();
  }
};

class RM1Settings : public AM1TypeSettings {
 public:
  static constexpr const char* settingsDescription = "RM1Settings";
  RM1Settings() : AM1TypeSettings(settingsDescription) {
    // Method
    Utils::UniversalSettings::StringDescriptor method("The method to be used.");
    method.setDefaultValue("rm1");
    _fields.push_back(Utils::SettingsNames::method, method);

    resetToDefaults();
  }
};

class PM3Settings : public AM1TypeSettings {
 public:
  static constexpr const char* settingsDescription = "PM3Settings";
  PM3Settings() : AM1TypeSettings(settingsDescription) {
    // Method
    Utils::UniversalSettings::StringDescriptor method("The method to be used.");
    method.setDefaultValue("pm3");
    _fields.push_back(Utils::SettingsNames::method, method);

    resetToDefaults();
  }
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_AM1SETTINGS_H
