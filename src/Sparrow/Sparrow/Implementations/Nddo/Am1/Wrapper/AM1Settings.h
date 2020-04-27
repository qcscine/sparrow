/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_AM1SETTINGS_H
#define SPARROW_AM1SETTINGS_H

#include <Sparrow/ParametersHeader.h>
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
  AM1TypeSettings(std::string settingsDescription, std::string parameterFolderName)
    : Settings(std::move(settingsDescription)) {
    auto fullParameterPath = std::move(parameterFolderName) + "/parameters.xml";

    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, fullParameterPath);
    Utils::UniversalSettings::BoolDescriptor useNDDODipoleApprox("Sets use of NDDO dipole approximation.");
    useNDDODipoleApprox.setDefaultValue(true);
    _fields.push_back(Utils::SettingsNames::NDDODipoleApproximation, std::move(useNDDODipoleApprox));

    resetToDefaults();
    modifyString(Utils::SettingsNames::parameterRootDirectory, std::string(parametersRootDir));
  }
};

class AM1Settings : public AM1TypeSettings {
 public:
  static constexpr const char* settingsDescription = "AM1Settings";
  static constexpr const char* parameterFolderName = "Am1";
  AM1Settings() : AM1TypeSettings(settingsDescription, parameterFolderName) {
  }
};

class RM1Settings : public AM1TypeSettings {
 public:
  static constexpr const char* settingsDescription = "RM1Settings";
  static constexpr const char* parameterFolderName = "Rm1";
  RM1Settings() : AM1TypeSettings(settingsDescription, parameterFolderName) {
  }
};

class PM3Settings : public AM1TypeSettings {
 public:
  static constexpr const char* settingsDescription = "PM3Settings";
  static constexpr const char* parameterFolderName = "Pm3";
  PM3Settings() : AM1TypeSettings(settingsDescription, parameterFolderName) {
  }
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_AM1SETTINGS_H
