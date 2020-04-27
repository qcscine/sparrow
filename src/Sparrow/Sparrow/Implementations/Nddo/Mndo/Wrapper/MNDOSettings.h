/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_MNDOCALCULATORSETTINGS_H
#define SPARROW_MNDOCALCULATORSETTINGS_H

#include <Sparrow/ParametersHeader.h>
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
    Utils::UniversalSettings::SettingPopulator::populateSemiEmpiricalSettings(_fields, "Mndo/parameters.xml");

    Utils::UniversalSettings::BoolDescriptor useNDDODipoleApprox("Sets use of NDDO dipole approximation.");
    useNDDODipoleApprox.setDefaultValue(true);
    _fields.push_back(Utils::SettingsNames::NDDODipoleApproximation, std::move(useNDDODipoleApprox));

    resetToDefaults();
    modifyString(Utils::SettingsNames::parameterRootDirectory, std::string(parametersRootDir));
  }
};

} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_MNDOCALCULATORSETTINGS_H
