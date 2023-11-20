/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_TDDFTBSETTINGS_H
#define SPARROW_TDDFTBSETTINGS_H

#include <Sparrow/Implementations/TimeDependent/LinearResponseSettings.h>

namespace Scine {
namespace Sparrow {
class TDDFTBSettings : public LinearResponseSettings {
 public:
  TDDFTBSettings() : LinearResponseSettings() {
    Utils::UniversalSettings::OptionListDescriptor prunedBasisCalculation(
        "Sets whether the basis of singly excited determinants should be pruned and with which method.");
    prunedBasisCalculation.addOption(Utils::SettingsNames::PruningOptions::none);
    prunedBasisCalculation.addOption(Utils::SettingsNames::PruningOptions::energy);
    prunedBasisCalculation.setDefaultOption(Utils::SettingsNames::PruningOptions::none);

    Utils::UniversalSettings::DoubleDescriptor energyThresholdForPruning(
        "Sets the threshold for pruning with an energy criterion in au.");
    energyThresholdForPruning.setMinimum(0.0);
    energyThresholdForPruning.setDefaultValue(0.0);

    Utils::UniversalSettings::DoubleDescriptor perturbationTheoryThresholdForPruning(
        "Sets the threshold for pruning with an intensity criterion in au.");
    perturbationTheoryThresholdForPruning.setMinimum(0.);
    perturbationTheoryThresholdForPruning.setDefaultValue(1e-4);

    Utils::UniversalSettings::BoolDescriptor TDAApproximation(
        "Switches on the TDA for the excited states calculation.");
    TDAApproximation.setDefaultValue(false);

    _fields.push_back(Utils::SettingsNames::pruneBasis, std::move(prunedBasisCalculation));
    _fields.push_back(Utils::SettingsNames::energyThreshold, std::move(energyThresholdForPruning));
    _fields.push_back(Utils::SettingsNames::perturbativeThreshold, std::move(perturbationTheoryThresholdForPruning));
    _fields.push_back("tda", std::move(TDAApproximation));
    resetToDefaults();
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_TDDFTBSETTINGS_H
