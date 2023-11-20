/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "IRCalculator.h"
#include "../SpectroscopySettings.h"
#include "../Utils/Spectrum.h"
#include "IntensitiesCalculator.h"
#include <Core/Interfaces/Calculator.h>
#include <Core/Log.h>
#include <Core/ModuleManager.h>
#include <Sparrow/Implementations/RealTimeSpectroscopy/Utils/LineWidthGenerator.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/GeometricDerivatives/NormalModeAnalysis.h>
#include <Utils/GeometricDerivatives/NumericalHessianCalculator.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/GeometryOptimization/GeometryOptimizer.h>
#include <Utils/IO/ChemicalFileFormats/XyzStreamHandler.h>
#include <Utils/Optimizer/GradientBased/Lbfgs.h>
#include <Utils/Optimizer/GradientBased/SteepestDescent.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <chrono>
#include <iomanip>

namespace Scine {
namespace Sparrow {
namespace RealTimeSpectroscopy {

IRCalculator::IRCalculator() {
  settings_ = std::make_unique<IRSettings>();
}

void IRCalculator::updateState(std::shared_ptr<Core::State> state) {
  if (calculator_) {
    calculator_->loadState(std::move(state));
  }
}

Spectrum IRCalculator::calculate(const Utils::PositionCollection& positions, int structureIndex, std::ostream& out) {
  calculator_->modifyPositions(positions);
  // Optimize or project after.
  auto startOpt = std::chrono::system_clock::now();
  if (settings_->getString(optimizationProfileOption) != "none") {
    std::shared_ptr<std::ofstream> optimizationOut =
        std::make_shared<std::ofstream>("optimization" + std::to_string(structureIndex) + ".out");
    *optimizationOut << std::setw(5) << "#" << std::string(44, '=') << "#" << std::endl;
    *optimizationOut << std::setw(5) << "#" << std::setw(20) << "Iteration" << std::setw(20) << "Energy" << std::setw(5)
                     << "#" << std::endl;
    *optimizationOut << std::setw(5) << "#" << std::string(44, '=') << "#" << std::endl;
    auto observer = [&](const int& it, const double& en, const Eigen::VectorXd& gr) -> void {
      *optimizationOut << std::setw(5) << "" << std::setw(20) << it << std::setw(20) << en << std::setw(5) << "" << std::endl;
    };
    calculator_->getLog().output.add("optimize_observer", optimizationOut);
    calculator_->modifyPositions(positions);
    Utils::GeometryOptimizer<Utils::Bfgs> optimizer(*calculator_);
    optimizer.addObserver(std::move(observer));
    auto profile = profileFactory(optimizer.getSettings(), settings_->getString(optimizationProfileOption));
    optimizer.setSettings(profile->toSettings());
    // optimizer.transformCoordinates = false;
    optimizer.optimize(*(calculator_->getStructure()), calculator_->getLog());
    *optimizationOut << "END!" << std::endl;
    optimizationOut->close();
    calculator_->getLog().output.remove("optimize_observer");
  }
  auto endOpt = std::chrono::system_clock::now();
  out << std::chrono::duration_cast<std::chrono::milliseconds>(endOpt - startOpt).count() << " ";

  auto startHes = std::chrono::system_clock::now();
  Utils::Results results = calculateHessianAndDipoleGradient();

  Utils::NormalModesContainer normalModesContainer;

  if (settings_->getBool(projectionOption)) {
    const Utils::GradientCollection& gradient = calculator_->calculate().get<Utils::Property::Gradients>().normalized();
    normalModesContainer = calculateFrequencies(results.get<Utils::Property::Hessian>(), gradient);
  }
  else {
    auto startAlign = std::chrono::system_clock::now();
    normalModesContainer = calculateFrequencies(results.get<Utils::Property::Hessian>());
    auto endAlign = std::chrono::system_clock::now();
    calculator_->getLog().output << "Time to diagonalize Hessian: "
                                 << std::chrono::duration_cast<std::chrono::milliseconds>(endAlign - startAlign).count()
                                 << Core::Log::nl;
  }

  auto w = normalModesContainer.getWaveNumbers();
  Eigen::VectorXd wavenumbers = Eigen::Map<const Eigen::VectorXd>(w.data(), w.size());

  Eigen::VectorXd intensities =
      calculateIntensities(results.get<Utils::Property::DipoleGradient>(), normalModesContainer.getNormalModes());

  std::pair<std::string, std::string> labels = {"Frequency / cm^-1", "Integral Absorption Coefficient / km mol^-1"};
  Spectrum spectrum{wavenumbers, intensities, labels};

  std::ofstream spectraIR("irLineSpectra" + std::to_string(structureIndex) + ".out");
  std::ofstream states("irStates" + std::to_string(structureIndex) + ".out");
  std::ofstream hessian("hessian" + std::to_string(structureIndex) + ".out");

  states << normalModesContainer.getNormalModes() << std::endl;
  hessian << *lastHessian_ << std::endl;
  for (int i = 0; i < spectrum.size(); ++i)
    spectraIR << spectrum.getXData(i) << " " << spectrum.getYData(i) << std::endl;
  LineWidthGenerator processor(spectrum);
  auto endHes = std::chrono::system_clock::now();
  out << std::chrono::duration_cast<std::chrono::milliseconds>(endHes - startHes).count() << std::endl;

  spectraIR.close();
  states.close();
  hessian.close();
  return processor.generateLorentzianProfile(settings_->getDouble(resolutionOption), settings_->getDouble(fwhmOption));
}

Utils::Results IRCalculator::calculateHessianAndDipoleGradient() {
  Utils::NumericalHessianCalculator hessianCalculator(*calculator_);
  hessianCalculator.requiredDipoleGradient(true);
  Utils::Results results;

  if (settings_->getBool(partialHessianOption) && lastPositions_ && lastHessian_) {
    auto startAlign = std::chrono::system_clock::now();
    auto indices = Utils::Geometry::Distances::getListOfDivergingAtomsRobust(
        *lastPositions_, calculator_->getPositions(), settings_->getDouble(partialHessianRMSDDeviationOption), 1e-2, 20,
        calculator_->getStructure()->getElements());
    auto endAlign = std::chrono::system_clock::now();
    calculator_->getLog().output << "Number of diverging atoms for partial Hessian calculation: " << indices.size()
                                 << Core::Log::nl;
    calculator_->getLog().output << "Time to align iteratively the structure: "
                                 << std::chrono::duration_cast<std::chrono::milliseconds>(endAlign - startAlign).count()
                                 << Core::Log::nl;
    startAlign = std::chrono::system_clock::now();
    auto partialResults = hessianCalculator.calculate(indices);
    endAlign = std::chrono::system_clock::now();
    calculator_->getLog().output << "Time to calculate Hessian: "
                                 << std::chrono::duration_cast<std::chrono::milliseconds>(endAlign - startAlign).count()
                                 << Core::Log::nl;
    for (int index : indices) {
      lastHessian_->middleCols(index * 3, 3) = partialResults.get<Utils::Property::Hessian>().middleCols(index * 3, 3);
      // This doubles some of the work, but really any other way right now would be premature opt. Well readable
      lastHessian_->middleRows(index * 3, 3) = partialResults.get<Utils::Property::Hessian>().middleRows(index * 3, 3);
      lastDipoleGradient_->row(index) = partialResults.get<Utils::Property::DipoleGradient>().row(index);
    }
    results.set<Utils::Property::Hessian>(*lastHessian_);
    results.set<Utils::Property::DipoleGradient>(*lastDipoleGradient_);
    lastPositions_ = std::make_unique<Utils::PositionCollection>(calculator_->getPositions());
  }
  else {
    results = hessianCalculator.calculate();
    lastPositions_ = std::make_unique<Utils::PositionCollection>(calculator_->getPositions());
    lastHessian_ = std::make_unique<Utils::HessianMatrix>(results.get<Utils::Property::Hessian>());
    lastDipoleGradient_ = std::make_unique<Utils::DipoleGradient>(results.get<Utils::Property::DipoleGradient>());
  }

  return results;
}

Utils::NormalModesContainer IRCalculator::calculateFrequencies(const Utils::HessianMatrix& hessian,
                                                               const Utils::GradientCollection& gradients) const {
  auto structure = calculator_->getStructure();

  if (settings_->getBool(projectionOption)) {
    // Project with gradient if not optimized
    return Utils::NormalModeAnalysis::calculateOrthogonalNormalModes(hessian, structure->getElements(),
                                                                     structure->getPositions(), gradients);
  }
  return Utils::NormalModeAnalysis::calculateNormalModes(hessian, structure->getElements(), structure->getPositions());
}

Eigen::VectorXd IRCalculator::calculateIntensities(const Utils::DipoleGradient& dipoleGradient,
                                                   const Eigen::MatrixXd& normalModes) const {
  Eigen::VectorXd squaredNormalDipoleGradient = IntensitiesCalculator::transformCartesianToSquaredNormalDipoleGradient(
      normalModes, dipoleGradient, Utils::Geometry::Properties::getMasses(calculator_->getStructure()->getElements()));

  IntensitiesCalculator intensityCalculator;
  intensityCalculator.setSquaredNormalDipoleGradient(squaredNormalDipoleGradient);

  return intensityCalculator.getAdsorptionCoefficients();
}

// Utils::HessianMatrix IRCalculator::projectHessianMatrix(Utils::HessianMatrix hessian, Utils::GradientCollection
// gradient) const {
//  Eigen::MatrixXd gradientMatrix = Eigen::Map<const Eigen::VectorXd>(gradient.data(), hessian.rows()) *
//                                   Eigen::Map<const Eigen::VectorXd>(gradient.data(), hessian.rows()).transpose();
//  Eigen::MatrixXd projector = Eigen::MatrixXd::Identity(hessian.rows(), hessian.cols()) - gradientMatrix;
//  return projector * std::move(hessian).selfadjointView<Eigen::Lower>() * projector;
//}

const Utils::Settings& IRCalculator::settings() const {
  return *settings_;
}

Utils::Settings& IRCalculator::settings() {
  return *settings_;
}

void IRCalculator::initialize(const Utils::ElementTypeCollection& elements) {
  auto& manager = Core::ModuleManager::getInstance();
  try {
    calculator_ = manager.get<Core::Calculator>(settings_->getString(Utils::SettingsNames::method));
  }
  catch (std::exception& e) {
    throw std::runtime_error("Method " + settings_->getString(Utils::SettingsNames::method) + " not found.");
  }
  calculator_->setRequiredProperties(Utils::Property::Energy | Utils::Property::Gradients | Utils::Property::Hessian |
                                     Utils::Property::DipoleGradient);
  if (calculator_->settings().valueExists(Utils::SettingsNames::selfConsistenceCriterion)) {
    calculator_->settings().modifyDouble(Utils::SettingsNames::selfConsistenceCriterion, 1e-9);
  }
  if (!settings_->getString(Utils::SettingsNames::methodParameters).empty()) {
    calculator_->settings().modifyString(Utils::SettingsNames::methodParameters,
                                         settings_->getString(Utils::SettingsNames::methodParameters));
  }
  Utils::PositionCollection positions = Utils::PositionCollection::Zero(elements.size(), 3);
  auto structure = Utils::AtomCollection{elements, positions};
  calculator_->setStructure(structure);
  lastPositions_.reset(nullptr);
  lastHessian_.reset(nullptr);
  lastDipoleGradient_.reset(nullptr);
  calculator_->getLog().output.add("cout", Core::Log::coutSink());
  calculator_->getLog().warning.add("warning", Core::Log::coutSink());
  calculator_->getLog().error.add("error", Core::Log::coutSink());
}

bool IRCalculator::gradientAllowsIRCalculation(const Utils::GradientCollection& gradient) const {
  double gradientThreshold = settings_->getDouble("gradient_threshold");
  return calculator_ && gradient.rowwise().norm().sum() < gradientThreshold;
}

void IRCalculator::modifyPositions(const Utils::PositionCollection& positions) {
  if (calculator_)
    calculator_->modifyPositions(positions);
}

std::unique_ptr<Utils::AtomCollection> IRCalculator::getOptimizedStructure() const {
  return calculator_->getStructure();
}
} // namespace RealTimeSpectroscopy
} // namespace Sparrow
} // namespace Scine
