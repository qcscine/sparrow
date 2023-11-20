/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SPARROW_SPECTROSCOPYSETTINGS_H
#define SPARROW_SPECTROSCOPYSETTINGS_H

#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/Optimizer/GradientBased/Bfgs.h>
#include <Utils/Optimizer/GradientBased/GradientBasedCheck.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/DoubleDescriptor.h>
#include <Utils/UniversalSettings/ParametrizedOptionValue.h>
#include <Utils/UniversalSettings/SettingPopulator.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <cassert>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

constexpr const char* resolutionOption = "spectrum_resolution";
constexpr const char* fwhmOption = "spectrum_fwhm";
constexpr const char* gradientThresholdOption = "gradient_threshold";
constexpr const char* partialHessianOption = "partial_hessian";
constexpr const char* projectionOption = "gradient_projection";
constexpr const char* partialHessianRMSDDeviationOption = "partial_hessian_RMSD_deviation";
constexpr const char* uvVisGuessPropagatorDiisDimension = "uv_vis_diis_max_dimension";
constexpr const char* optimizationProfileOption = "optimization_profile";

class SpectroscopySettings : public Utils::Settings {
 public:
  SpectroscopySettings() : Utils::Settings("Spectroscopy Settings") {
  }
  explicit SpectroscopySettings(std::string settingName) : Utils::Settings(std::move(settingName)) {
    resetToDefaults();
  }

 protected:
  void addResolutionOption(double minValue, double maxValue, double defaultValue) {
    assert(minValue < maxValue && "Min value must be smaller than max value.");
    assert(minValue < defaultValue && "Min value must be smaller than default value.");
    assert(defaultValue < maxValue && "default value must be smaller than max value.");
    Utils::UniversalSettings::DoubleDescriptor resolution("The resolution the spectrum should have."
                                                          "For IR, this is in cm^-1, for UV/Vis in eV.");
    resolution.setMinimum(minValue);
    resolution.setMaximum(maxValue);
    resolution.setDefaultValue(defaultValue);
    _fields.push_back(resolutionOption, std::move(resolution));
  }
  void addFWHMOption(double defaultFwhm) {
    assert(defaultFwhm < 100. && "Min value must be smaller than default value.");
    assert(defaultFwhm > 0. && "default value must be smaller than max value.");
    Utils::UniversalSettings::DoubleDescriptor fwhm("The fwhm of the peaks in the spectrum."
                                                    "For IR, this is in cm^-1, for UV/Vis in eV.");
    fwhm.setMinimum(0.);
    fwhm.setMaximum(100.);
    fwhm.setDefaultValue(defaultFwhm);
    _fields.push_back(fwhmOption, std::move(fwhm));
  }
  void addMethodOption() {
    Utils::UniversalSettings::StringDescriptor parameterFile("Path to the parameter file");
    parameterFile.setDefaultValue("");
    _fields.push_back(Utils::SettingsNames::methodParameters, std::move(parameterFile));

    Utils::UniversalSettings::StringDescriptor method("The method used in the digester.");
    method.setDefaultValue("PM6");
    _fields.push_back(Utils::SettingsNames::method, std::move(method));

    Utils::UniversalSettings::SettingPopulator::populateScfSettings(_fields);
  }
  void addPruningOption() {
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
  }
};

class IRSettings : public SpectroscopySettings {
 public:
  IRSettings() : SpectroscopySettings("IR Settings") {
    addResolutionOption(0.1, 10, 1);
    addFWHMOption(30.);
    addMethodOption();

    Utils::UniversalSettings::BoolDescriptor gradientProjection("Flags the use of the gradient projection technique.");
    _fields.push_back(projectionOption, std::move(gradientProjection));

    Utils::UniversalSettings::OptionListDescriptor optimizationProfile(
        "How tight the structure is optimized before projection. VeryTight, Tight, Medium, Loose, VeryLoose.");
    optimizationProfile.addOption("very_tight");
    optimizationProfile.addOption("tight");
    optimizationProfile.addOption("medium");
    optimizationProfile.addOption("loose");
    optimizationProfile.addOption("very_loose");
    optimizationProfile.addOption("none");
    optimizationProfile.setDefaultOption("tight");

    _fields.push_back(optimizationProfileOption, std::move(optimizationProfile));

    Utils::UniversalSettings::DoubleDescriptor gradientThreshold("The gradient threshold for minima detection (au).");
    gradientThreshold.setMinimum(1.0e-4);
    gradientThreshold.setMaximum(1.0);
    gradientThreshold.setDefaultValue(1.0e-2);

    _fields.push_back(gradientThresholdOption, std::move(gradientThreshold));

    Utils::UniversalSettings::BoolDescriptor partialHessianCalculation(
        "Calculate only parts of the Hessian matrix corresponding"
        "to atoms that moved consistently wrt last calculation.");
    partialHessianCalculation.setDefaultValue(false);
    _fields.push_back(partialHessianOption, std::move(partialHessianCalculation));

    Utils::UniversalSettings::DoubleDescriptor partialHessianRMSD(
        "Deviation atoms in aligned structures must have to be"
        "considered diverging.");
    partialHessianRMSD.setMinimum(1.0e-2);
    partialHessianRMSD.setMaximum(1.0);
    partialHessianRMSD.setDefaultValue(0.2);

    _fields.push_back(partialHessianRMSDDeviationOption, std::move(partialHessianRMSD));

    resetToDefaults();
  }
};
class UvVisSettings : public SpectroscopySettings {
 public:
  UvVisSettings() : SpectroscopySettings("Uv/Vis Settings") {
    addResolutionOption(0.0001, 0.3, 0.002);
    addFWHMOption(0.3);
    addMethodOption();
    addPruningOption();

    Utils::UniversalSettings::IntDescriptor numberEigenvalues("Number of roots to be calculated");
    numberEigenvalues.setMinimum(1);
    numberEigenvalues.setDefaultValue(20);
    _fields.push_back(Utils::SettingsNames::numberOfEigenstates, std::move(numberEigenvalues));

    Utils::UniversalSettings::IntDescriptor initialSubspaceDimension("Dimension of the guess subspace");
    numberEigenvalues.setMinimum(1);
    numberEigenvalues.setDefaultValue(1);
    _fields.push_back(Utils::SettingsNames::initialSubspaceDimension, std::move(initialSubspaceDimension));

    Utils::UniversalSettings::OptionListDescriptor gepAlgorithm(
        "Algorithm to compute the stable generalized eigenvalue problem.");
    gepAlgorithm.addOption("standard");
    gepAlgorithm.addOption("cholesky");
    gepAlgorithm.addOption("simultaneous_diag");
    gepAlgorithm.setDefaultOption("simultaneous_diag");
    _fields.push_back("gep_algo", std::move(gepAlgorithm));

    Utils::UniversalSettings::IntDescriptor maxDiisDimension("Maximal dimension of the diis subspace");
    maxDiisDimension.setMinimum(1);
    maxDiisDimension.setDefaultValue(4);
    _fields.push_back(uvVisGuessPropagatorDiisDimension, std::move(maxDiisDimension));

    resetToDefaults();
  }
};

struct GeometryOptimizationProfile {
  GeometryOptimizationProfile(Utils::Settings settings) : settings_(std::move(settings)) {
  }
  virtual ~GeometryOptimizationProfile() = default;
  Utils::Settings toSettings() {
    applyProfile(settings_);
    return settings_;
  }
  virtual void applyProfile(Utils::Settings& settings) = 0;

 private:
  Utils::Settings settings_;
};

struct VeryTightOptimizationProfile : public GeometryOptimizationProfile {
  VeryTightOptimizationProfile(const Utils::Settings& settings) : GeometryOptimizationProfile(settings){};
  ~VeryTightOptimizationProfile() final = default;
  void applyProfile(Utils::Settings& settings) final {
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepMaxCoeff, 2.0e-5);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepRMS, 1.0e-5);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradMaxCoeff, 2.0e-5);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradRMS, 1.0e-5);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::deltaValue, 1.0e-7);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::maxIter, 1000);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::requirement, 4);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Bfgs::gdiisMaxStore, 10);
  }
};
struct TightOptimizationProfile : public GeometryOptimizationProfile {
  TightOptimizationProfile(const Utils::Settings& settings) : GeometryOptimizationProfile(settings){};
  ~TightOptimizationProfile() final = default;
  void applyProfile(Utils::Settings& settings) final {
    settings.modifyInt(Utils::SettingsNames::Optimizations::Bfgs::gdiisMaxStore, 8);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::maxIter, 1000);
  }
};
struct MediumOptimizationProfile : public GeometryOptimizationProfile {
  MediumOptimizationProfile(const Utils::Settings& settings) : GeometryOptimizationProfile(settings){};
  ~MediumOptimizationProfile() final = default;
  void applyProfile(Utils::Settings& settings) final {
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepMaxCoeff, 5.0e-4);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepRMS, 5.0e-3);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradMaxCoeff, 5.0e-4);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradRMS, 5.0e-4);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::deltaValue, 5.0e-6);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::maxIter, 100);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::requirement, 2);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Bfgs::gdiisMaxStore, 8);
  }
};
struct LooseOptimizationProfile : public GeometryOptimizationProfile {
  LooseOptimizationProfile(const Utils::Settings& settings) : GeometryOptimizationProfile(settings){};
  ~LooseOptimizationProfile() final = default;
  void applyProfile(Utils::Settings& settings) final {
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepMaxCoeff, 5.0e-3);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepRMS, 1.0e-2);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradMaxCoeff, 1.0e-3);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradRMS, 5.0e-3);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::deltaValue, 1.0e-5);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::maxIter, 100);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::requirement, 2);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Bfgs::gdiisMaxStore, 8);
  }
};
struct VeryLooseOptimizationProfile : public GeometryOptimizationProfile {
  VeryLooseOptimizationProfile(const Utils::Settings& settings) : GeometryOptimizationProfile(settings){};
  ~VeryLooseOptimizationProfile() final = default;
  void applyProfile(Utils::Settings& settings) final {
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepMaxCoeff, 1.0e-2);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::stepRMS, 5.0e-2);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradMaxCoeff, 5.0e-3);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::gradRMS, 1.0e-2);
    settings.modifyDouble(Utils::SettingsNames::Optimizations::Convergence::deltaValue, 1.0e-4);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::maxIter, 50);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Convergence::requirement, 2);
    settings.modifyInt(Utils::SettingsNames::Optimizations::Bfgs::gdiisMaxStore, 8);
  }
};

inline std::unique_ptr<GeometryOptimizationProfile> profileFactory(const Utils::Settings& settings, const std::string& profile) {
  if (profile == "very_tight")
    return std::make_unique<VeryTightOptimizationProfile>(settings);
  else if (profile == "tight")
    return std::make_unique<TightOptimizationProfile>(settings);
  else if (profile == "medium")
    return std::make_unique<MediumOptimizationProfile>(settings);
  else if (profile == "loose")
    return std::make_unique<LooseOptimizationProfile>(settings);
  else if (profile == "very_loose")
    return std::make_unique<VeryLooseOptimizationProfile>(settings);
  else
    throw std::runtime_error("Geometry optimization profile '" + profile +
                             "' not found."
                             "Possible profiles: very_tight, tight, medium, loose, very_loose");
}

} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine

#endif // SPARROW_SPECTROSCOPYSETTINGS_H
