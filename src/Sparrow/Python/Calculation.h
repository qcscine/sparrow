/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Core/Interfaces/Calculator.h>
#include <Core/ModuleManager.h>
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/CalculatorBasics/Results.h>
#include <Utils/Constants.h>
#include <Utils/Geometry.h>
#include <Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h>
#include <Utils/Settings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

using namespace Scine;

namespace Sparrow {
namespace Python {

class Calculation {
 public:
  /**
   * @brief Set up a new calculation.
   * @param method The method to be used in the calculation.
   */
  explicit Calculation(const std::string& method = "PM6");

  /**
   * @brief Specify settings for the calculation, e.g., molecular charge or spin multiplicity.
   * @param settings The settings to be applied in the calculation
   */
  void setSettings(pybind11::dict settings);

  /**
   * @brief Get the current settings of the calculation.
   * @return The current settings of the calculation.
   */
  std::map<std::string, std::string> getSettings();

  /**
   * @brief Specify a molecule to be calculated from an XYZ file.
   * @param structureFile The path of the file containing the molecular structure (XYZ format).
   */
  void setStructure(const std::string& structureFile);

  /**
   * @brief Specify the atomic elements of the molecule to be calculated.
   * @param elements The atomic elements of the molecule to be calculated.
   */
  void setElements(pybind11::list elements);

  /**
   * @brief Returns the chemical elements of the molecule to be calculated.
   * @return The chemical elements of the system to be calculated.
   */
  std::vector<std::string> getElements();

  /**
   * @brief Specify the Cartesian coordinates of the nuclei of the molecule to be calculated.
   * @param positions The Cartesian coordinates of the nuclei of the molecule to be calculated (in Angstrom).
   */
  void setPositions(const Utils::PositionCollection& positions);

  /**
   * @brief Returns the Cartesian coordinates of the nuclei of the molecule to be calculated.
   * @return The Cartesian coordinates of the nuclei of the molecule to be calculated (in bohr).
   */
  Utils::PositionCollection getPositions();

  /**
   * @brief Calculate and return the energy.
   * @return The total electronic energy (in hartree).
   */
  double calculateEnergy();

  /**
   * @brief Calculate and return the nuclear gradients.
   * @return The nuclear gradients (in hartree/bohr).
   */
  Utils::GradientCollection calculateGradients();

  /**
   * @brief Calculate and return the Hessian matrix.
   * @return The Hessian matrix (in hartree/bohr^2).
   */
  Utils::HessianMatrix calculateHessian();

  /**
   * @brief Calculate and return the bond order matrix.
   * @return The bond order matrix (dimensionless).
   */
  Eigen::SparseMatrix<double> calculateBondOrders();

  /**
   * @brief Extracts the property 'SuccessfulCalculation' from the results.
   * @return Whether the calculation converged.
   */
  bool isConverged() const;

 private:
  std::shared_ptr<Core::Calculator> calculator_;
  Utils::ElementTypeCollection elementTypeCollection_;
  Utils::PositionCollection positionCollection_;

  const Utils::Results& calculate(Utils::Property property);
};

} // namespace Python
} // namespace Sparrow

PYBIND11_MODULE(scine_sparrow, m) {
  m.doc() = "Pybind11 Bindings for SCINE Sparrow";

  pybind11::class_<Sparrow::Python::Calculation>(m, "Calculation")
      .def(pybind11::init<const std::string&>(), "initializes a new calculation", pybind11::arg("method") = "PM6")
      .def("set_settings", &Sparrow::Python::Calculation::setSettings, "sets the specified settings for the calculation")
      .def("get_settings", &Sparrow::Python::Calculation::getSettings, "returns the current calculation settings")
      .def("set_structure", &Sparrow::Python::Calculation::setStructure, "specifies the system to be calculated from an XYZ file")
      .def("set_elements", &Sparrow::Python::Calculation::setElements,
           "specifies the chemical elements of the atoms present in the system to be calculated")
      .def("get_elements", &Sparrow::Python::Calculation::getElements,
           "get the chemical elements of the atoms present in the system to be calculated")
      .def("set_positions", &Sparrow::Python::Calculation::setPositions,
           "specifies the Cartesian coordinates of the atoms present in the system to be calculated (in Angstrom)")
      .def("get_positions", &Sparrow::Python::Calculation::getPositions,
           "get the Cartesian coordinates of the atoms present in the system to be calculated (in bohr)")
      .def("calculate_energy", &Sparrow::Python::Calculation::calculateEnergy, "calculates the total electronic energy (in hartree)")
      .def("calculate_gradients", &Sparrow::Python::Calculation::calculateGradients,
           "calculates the nuclear gradients (in hartree/bohr)")
      .def("calculate_hessian", &Sparrow::Python::Calculation::calculateHessian, "calculates the Hessian matrix (in hartree/bohr^2)")
      .def("calculate_bond_orders", &Sparrow::Python::Calculation::calculateBondOrders, "calculates the bond order matrix")
      .def("is_converged", &Sparrow::Python::Calculation::isConverged, "returns whether calculation converged");
}
