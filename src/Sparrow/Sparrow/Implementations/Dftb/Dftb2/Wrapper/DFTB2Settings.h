/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB2SETTINGS_H
#define SPARROW_DFTB2SETTINGS_H

#include <Sparrow/ParametersHeader.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Sparrow {

/**
 * @class DFTB2Settings DFTB2Settings.h
 * @brief The Settings specific to the DFTB2 (SCC-DFTB) method.
 */
class DFTB2Settings : public Scine::Utils::Settings {
 public:
  DFTB2Settings() : Settings("DFTB2Settings") {
    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "Dftb/mio-1-1/");

    resetToDefaults();
    modifyString(Utils::SettingsNames::parameterRootDirectory, std::string(parametersRootDir));
  };
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DFTB2SETTINGS_H
