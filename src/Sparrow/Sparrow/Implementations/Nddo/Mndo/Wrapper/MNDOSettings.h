/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MNDOCALCULATORSETTINGS_H
#define SPARROW_MNDOCALCULATORSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Sparrow {

/**
 * @class MNDOCalculatorSettings MNDOCalculatorSettings.h
 * @brief The Settings specific to the MNDOCalculator.
 */
class MNDOSettings : public Scine::Utils::Settings {
 public:
  MNDOSettings() : Settings("MNDOSettings") {
    Utils::UniversalSettings::SettingPopulator::populateLcaoSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "");

    Utils::UniversalSettings::BoolDescriptor useNDDODipoleApprox("Sets use of NDDO dipole approximation.");
    useNDDODipoleApprox.setDefaultValue(true);
    _fields.push_back(Utils::SettingsNames::NDDODipoleApproximation, std::move(useNDDODipoleApprox));

    // Method
    Utils::UniversalSettings::StringDescriptor method("The method to be used.");
    method.setDefaultValue("mndo");
    _fields.push_back(Utils::SettingsNames::method, method);

    resetToDefaults();
  }
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_MNDOCALCULATORSETTINGS_H
