/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_DFTB0SETTINGS_H
#define SPARROW_DFTB0SETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Sparrow {

/**
 * @class DFTB0Settings DFTB0Settings.h
 * @brief The Settings specific to the DFTB0 method, a non SCF method.
 * Please note that since DFTB0 is not an SCF method, it will not contain SCF options. Furthermore, it is not capable
 * of performing calculation in unrestricted formalism.
 */
class DFTB0Settings : public Scine::Utils::Settings {
 public:
  DFTB0Settings() : Settings("DFTB0Settings") {
    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "3ob-3-1");

    // Method
    Utils::UniversalSettings::StringDescriptor method("The method to be used.");
    method.setDefaultValue("dftb0");
    _fields.push_back(Utils::SettingsNames::method, method);

    resetToDefaults();
  };
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_DFTB0SETTINGS_H
