/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_PM6SETTINGS_H
#define SPARROW_PM6SETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Sparrow {

/**
 * @class PM6Settings PM6Settings.h
 * @brief The Settings specific to the PM6 method.
 */
class PM6Settings : public Scine::Utils::Settings {
 public:
  PM6Settings() : Settings("PM6MethodWrapper") {
    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "");

    Utils::UniversalSettings::BoolDescriptor useNDDODipoleApprox("Sets use of NDDO dipole approximation.");
    useNDDODipoleApprox.setDefaultValue(true);
    _fields.push_back(Utils::SettingsNames::NDDODipoleApproximation, std::move(useNDDODipoleApprox));

    // Method
    Utils::UniversalSettings::StringDescriptor method("The method to be used.");
    method.setDefaultValue("pm6");
    _fields.push_back(Utils::SettingsNames::method, method);

    resetToDefaults();
  };
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_PM6SETTINGS_H
