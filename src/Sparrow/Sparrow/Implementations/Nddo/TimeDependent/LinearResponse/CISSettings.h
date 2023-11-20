/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_CISSETTINGS_H
#define SPARROW_CISSETTINGS_H

#include <Sparrow/Implementations/TimeDependent/LinearResponseSettings.h>

namespace Scine {
namespace Sparrow {
class CISSettings : public LinearResponseSettings {
 public:
  CISSettings() : LinearResponseSettings() {
    // TODO: set default path to parameter file !!!!
    Utils::UniversalSettings::FileDescriptor excitedStatesParamFile("Sets path to excited states parameter files");
    excitedStatesParamFile.setDefaultValue("");

    Utils::UniversalSettings::DoubleDescriptor distanceThreshold(
        "Set the distance threshold after which no interaction is calculated (in Angstrom).");
    distanceThreshold.setMinimum(0.0);
    distanceThreshold.setDefaultValue(std::numeric_limits<double>::max());

    _fields.push_back(Utils::SettingsNames::excitedStatesParamFile, std::move(excitedStatesParamFile));
    _fields.push_back("distance_threshold", std::move(distanceThreshold));
    resetToDefaults();
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_CISSETTINGS_H
