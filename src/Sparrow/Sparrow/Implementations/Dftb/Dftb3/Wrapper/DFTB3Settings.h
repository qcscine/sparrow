/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB3SETTINGS_H
#define SPARROW_DFTB3SETTINGS_H

#include <Sparrow/ParametersHeader.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Sparrow {

/**
 * @class DFTB3Settings DFTB3Settings.h
 * @brief The Settings specific to the DFTB3 method.
 */
class DFTB3Settings : public Scine::Utils::Settings {
 public:
  DFTB3Settings() : Settings("DFTB3Settings") {
    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "Dftb/3ob-3-1/");

    resetToDefaults();
    modifyString(Utils::SettingsNames::parameterRootDirectory, std::string(parametersRootDir));
  }
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DFTB3CALCULATORSETTINGS_H
