/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Calculation.h"

namespace Sparrow {
namespace Python {

Calculation::Calculation(const std::string& method) {
  // Load module manager
  auto& manager = Core::ModuleManager::getInstance();

  // Set up calculator
  try {
    calculator_ = manager.get<Core::Calculator>(method, "Sparrow");
  }
  catch (...) {
    std::cout << "No SCINE module providing '" << method << "' is currently loaded.\n";
    std::cout << "Please add the module to the SCINE_MODULE_PATH in order for it to be accessible.\n";
    throw std::runtime_error("Failed to load method/program.");
  }
}

void Calculation::setSettings(pybind11::dict settings) {
  for (auto setting : settings) {
    std::string settingName = setting.first.cast<std::string>();

    if (settingName == Utils::SettingsNames::molecularCharge || settingName == Utils::SettingsNames::spinMultiplicity ||
        settingName == Utils::SettingsNames::maxIterations) {
      int value = setting.second.cast<int>();
      calculator_->settings().modifyInt(settingName, value);
    }

    else if (settingName == Utils::SettingsNames::unrestrictedCalculation) {
      bool value = setting.second.cast<bool>();
      calculator_->settings().modifyBool(settingName, value);
    }

    else if (settingName == Utils::SettingsNames::selfConsistanceCriterion) {
      double value = setting.second.cast<double>();
      calculator_->settings().modifyDouble(settingName, value);
    }

    else if (settingName == Utils::SettingsNames::mixer || settingName == Utils::SettingsNames::parameterFile ||
             settingName == Utils::SettingsNames::parameterRootDirectory) {
      std::string value = setting.second.cast<std::string>();
      calculator_->settings().modifyString(settingName, value);
    }

    else {
      throw std::runtime_error("The setting '" + settingName + "' cannot be set");
    }
  }
}

std::map<std::string, std::string> Calculation::getSettings() {
  std::map<std::string, std::string> settingsMap;
  auto settings = calculator_->settings();
  auto keys = settings.getKeys();

  for (const auto& key : keys) {
    auto value = settings.getValue(key);
    if (value.isBool()) {
      settingsMap[key] = std::to_string(value.toBool());
    }
    if (value.isDouble()) {
      settingsMap[key] = std::to_string(value.toDouble());
    }
    if (value.isInt()) {
      settingsMap[key] = std::to_string(value.toInt());
    }
    if (value.isString()) {
      settingsMap[key] = value.toString();
    }
  }

  return settingsMap;
}

void Calculation::setStructure(const std::string& structureFile) {
  auto structure = Utils::ChemicalFileHandler::read(structureFile);
  positionCollection_ = structure.first.getPositions();
  elementTypeCollection_ = structure.first.getElements();
}

void Calculation::setElements(pybind11::list elements) {
  std::string elementString;
  Utils::ElementType elementType;
  for (auto element : elements) {
    elementString = element.cast<std::string>();
    elementType = Utils::ElementInfo::elementTypeForSymbol(elementString);
    elementTypeCollection_.push_back(elementType);
  }
}

std::vector<std::string> Calculation::getElements() {
  std::vector<std::string> elements;
  for (auto elementType : elementTypeCollection_) {
    std::string element = Utils::ElementInfo::symbol(elementType);
    elements.push_back(element);
  }
  return elements;
}

void Calculation::setPositions(const Utils::PositionCollection& positions) {
  positionCollection_ = positions * Utils::Constants::bohr_per_angstrom;
}

Utils::PositionCollection Calculation::getPositions() {
  return positionCollection_;
}

double Calculation::calculateEnergy() {
  const Utils::Results& results = calculate(Utils::Property::Energy);
  return results.get<Utils::Property::Energy>();
}

Utils::GradientCollection Calculation::calculateGradients() {
  const Utils::Results& results = calculate(Utils::Property::Gradients);
  return results.get<Utils::Property::Gradients>();
}

Utils::HessianMatrix Calculation::calculateHessian() {
  const Utils::Results& results = calculate(Utils::Property::Hessian);
  return results.get<Utils::Property::Hessian>();
}

Eigen::SparseMatrix<double> Calculation::calculateBondOrders() {
  const Utils::Results& results = calculate(Utils::Property::BondOrderMatrix);
  const auto& bondOrderCollection = results.get<Utils::Property::BondOrderMatrix>();
  return bondOrderCollection.getMatrix();
}

const Utils::Results& Calculation::calculate(Utils::Property property) {
  Utils::PropertyList requiredProperties;
  requiredProperties.addProperty(property);
  calculator_->setRequiredProperties(requiredProperties);

  calculator_->setStructure({elementTypeCollection_, positionCollection_});

  return calculator_->calculate();
}

bool Calculation::isConverged() const {
  return calculator_->results().get<Utils::Property::SuccessfulCalculation>();
}

} // namespace Python
} // namespace Sparrow
