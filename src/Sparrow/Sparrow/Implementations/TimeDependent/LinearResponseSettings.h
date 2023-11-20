/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef SPARROW_LINEARRESPONSESETTINGS_H
#define SPARROW_LINEARRESPONSESETTINGS_H

#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>

namespace {
static constexpr const char* convergence = "convergence";
}

namespace Scine {
namespace Sparrow {
/**
 * @brief Generic settings for linear response methods, i.e. CIS and TD-DFTB.
 * @class LinearResponseSettings @file LinearResponseSettings.h
 */
class LinearResponseSettings : public Utils::Settings {
 public:
  LinearResponseSettings() : Settings("LinearResponseSettings") {
    Utils::UniversalSettings::IntDescriptor numberOfEigenstates("Sets the number of desired Eigenstates.");
    numberOfEigenstates.setMinimum(0);
    numberOfEigenstates.setDefaultValue(1);

    Utils::UniversalSettings::IntDescriptor initialSubspaceDimension(
        "Sets the initial space dimension to use in the Davidson diagonalizer.");
    initialSubspaceDimension.setMinimum(1);
    initialSubspaceDimension.setDefaultValue(1);

    Utils::UniversalSettings::IntDescriptor maxDavidsonIterations(
        "Sets the maximal iteration number for the iterative diagonalizer.");
    initialSubspaceDimension.setMinimum(0);
    initialSubspaceDimension.setDefaultValue(0);

    Utils::UniversalSettings::OptionListDescriptor spinBlock(
        "Determines the spin block to calculate in the TD-DFTB Configuration State Functions implementation.");
    spinBlock.addOption(Utils::SettingsNames::SpinBlocks::singlet);
    spinBlock.addOption(Utils::SettingsNames::SpinBlocks::triplet);
    spinBlock.addOption(Utils::SettingsNames::SpinBlocks::singletAndTriplet);
    spinBlock.setDefaultOption(Utils::SettingsNames::SpinBlocks::singlet);

    Utils::UniversalSettings::DoubleDescriptor maxMemory(
        "Sets maximum amount of memory that is allowed for the calculation.");
    maxMemory.setDefaultValue(4.0);

    Utils::UniversalSettings::DoubleDescriptor convergenceCriterionDavidson(
        "The convergence threshold for the Davidson solver.");
    convergenceCriterionDavidson.setDefaultValue(1.0e-5);
    convergenceCriterionDavidson.setMinimum(std::numeric_limits<double>::min());

    Utils::UniversalSettings::OptionListDescriptor gepAlgorithm("Algorithm to compute the stable Generalize"
                                                                "Eigenvalue Problem Ax=lBx when B is almost singular.");
    gepAlgorithm.addOption("standard");
    gepAlgorithm.addOption("cholesky");
    gepAlgorithm.addOption("simultaneous_diag");
    gepAlgorithm.setDefaultOption("standard");

    _fields.push_back(Utils::SettingsNames::numberOfEigenstates, std::move(numberOfEigenstates));
    _fields.push_back(Utils::SettingsNames::initialSubspaceDimension, std::move(initialSubspaceDimension));
    _fields.push_back(Utils::SettingsNames::spinBlock, std::move(spinBlock));
    _fields.push_back(Utils::SettingsNames::maxMemory, std::move(maxMemory));
    _fields.push_back(Utils::SettingsNames::maxDavidsonIterations, std::move(maxDavidsonIterations));
    _fields.push_back(convergence, std::move(convergenceCriterionDavidson));
    _fields.push_back("gep_algo", std::move(gepAlgorithm));

    resetToDefaults();
  }
};
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_LINEARRESPONSESETTINGS_H
