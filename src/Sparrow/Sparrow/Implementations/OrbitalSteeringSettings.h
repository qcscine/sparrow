/**
 * @file Calculator.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_ORBITALSTEERINGSETTINGS_H
#define SPARROW_ORBITALSTEERINGSETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>

namespace Scine {
namespace Sparrow {

static constexpr const char* numberOrbitalsToMixKey = "number_orbitals_to_mix";
static constexpr const char* mixingFrequencyKey = "mixing_frequency";
static constexpr const char* minimalAngleKey = "minimal_mixing_angle";
static constexpr const char* maximalAngleKey = "maximal_mixing_angle";
static constexpr const char* numberOrbitalsToConsiderKey = "number_orbitals_to_consider";

class OrbitalSteeringSettings : public Utils::Settings {
 public:
  static constexpr const int SameNumberOfIterations = 0;
  static constexpr const int AllOrbitals = 0;
  static constexpr const char* SameMixer = "same_mixer";
  OrbitalSteeringSettings() : Utils::Settings() {
    Utils::UniversalSettings::IntDescriptor numberOrbitalsToMix("Sets the number of orbitals that will be mixed.");
    numberOrbitalsToMix.setDefaultValue(10);
    numberOrbitalsToMix.setMinimum(0);
    _fields.push_back(numberOrbitalsToMixKey, std::move(numberOrbitalsToMix));

    Utils::UniversalSettings::IntDescriptor mixingFrequency(
        "Sets after how many single point calculations a mixing occurs.");
    mixingFrequency.setDefaultValue(5);
    mixingFrequency.setMinimum(1);
    _fields.push_back(mixingFrequencyKey, std::move(mixingFrequency));

    Utils::UniversalSettings::DoubleDescriptor minimalAngle("Sets the minimal angle for the mixing.");
    minimalAngle.setDefaultValue(0);
    minimalAngle.setMinimum(0);
    minimalAngle.setMaximum(90);
    _fields.push_back(minimalAngleKey, std::move(minimalAngle));

    Utils::UniversalSettings::DoubleDescriptor maximalAngle("Sets the maximal angle for the mixing.");
    maximalAngle.setDefaultValue(90);
    maximalAngle.setMinimum(0);
    maximalAngle.setMaximum(90);
    _fields.push_back(maximalAngleKey, std::move(maximalAngle));

    Utils::UniversalSettings::IntDescriptor numberOrbitalsToConsider(
        "Sets the number of orbitals from which to sample the mixing pairs.");
    numberOrbitalsToConsider.setDefaultValue(AllOrbitals);
    numberOrbitalsToConsider.setMinimum(AllOrbitals);
    _fields.push_back(numberOrbitalsToConsiderKey, std::move(numberOrbitalsToConsider));

    Utils::UniversalSettings::IntDescriptor maxScfIterations("Maximal number of iterations to reach self consistence.");
    maxScfIterations.setMinimum(SameNumberOfIterations);
    maxScfIterations.setDefaultValue(SameNumberOfIterations);

    _fields.push_back(Utils::SettingsNames::maxScfIterations, std::move(maxScfIterations));

    Utils::UniversalSettings::OptionListDescriptor mixer("Convergence acceleration to use.");
    mixer.addOption(SameMixer);
    mixer.addOption(Utils::SettingsNames::ScfMixers::noMixer);
    mixer.addOption(Utils::SettingsNames::ScfMixers::diis);
    mixer.addOption(Utils::SettingsNames::ScfMixers::ediis);
    mixer.addOption(Utils::SettingsNames::ScfMixers::ediisDiis);
    mixer.setDefaultOption(SameMixer);

    _fields.push_back(Utils::SettingsNames::mixer, mixer);

    resetToDefaults();
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_ORBITALSTEERINGSETTINGS_H
