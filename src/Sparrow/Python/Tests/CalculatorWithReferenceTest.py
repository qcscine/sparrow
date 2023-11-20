__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest
import scine_utilities as su
import scine_sparrow


class TestCalculatorWithReferencePythonBindings(unittest.TestCase):
    def setUp(self):
        structure = su.AtomCollection(4)
        structure.elements = [su.ElementType.C, su.ElementType.O, su.ElementType.H, su.ElementType.H]
        # in angstrom
        structure.positions = [[ 0.000000,    0.000000,    0.000000],
                         [ 0.000000,    0.000000,    1.212200],
                         [ 0.937197,    0.000000,   -0.584262],
                         [-0.937197,    0.000000,   -0.584262]]
        methane = su.AtomCollection(5)
        methane.elements = [su.ElementType.C, su.ElementType.H, su.ElementType.H,
                          su.ElementType.H, su.ElementType.H]
        # in angstrom
        methane.positions = [[0.0000000000,    0.0000000000,    0.0000000000],
                            [ 0.7,             0.72,            0.64        ],
                            [-0.635,          -0.635,           0.635       ],
                            [-0.71,            0.7287000000,   -0.7287000000],
                            [ 0.6287000000,   -0.6287000000,   -0.6287000000]]

        # transform to bohr
        structure.positions = structure.positions * su.BOHR_PER_ANGSTROM
        methane.positions = methane.positions * su.BOHR_PER_ANGSTROM
        manager = su.core.ModuleManager.get_instance()

        self.calculator = manager.get('calculator', 'PM3')
        self.calculator.settings["self_consistence_criterion"] = 1e-9
        self.calculator.structure = structure
        self.dftb = manager.get('calculator', 'DFTB2')
        self.dftb.structure = methane

        self.excited_states_calculator = manager.get('calculator_with_reference','CIS-NDDO')
        self.excited_states_calculator.reference_calculator = self.calculator
        self.excited_states_calculator.log = su.core.Log.silent()
        self.tddftb = manager.get('calculator_with_reference','TD-DFTB')
        self.tddftb.reference_calculator = self.dftb
        self.tddftb.log = su.core.Log.silent()

    def test_reference_calculation(self):
        self.calculator.set_required_properties([su.Property.Energy])
        self.excited_states_calculator.reference_calculation()
        energy = self.excited_states_calculator.reference_calculator.get_results().energy
        ref_energy = self.calculator.calculate().energy

        self.assertAlmostEqual(energy, ref_energy, places=6)

    def test_singlet_excited_states(self):
        self.excited_states_calculator.reference_calculation()
        self.excited_states_calculator.settings["number_eigenstates"] = 10
        results = self.excited_states_calculator.calculate()
        singlet_energies = results.excited_states.singlet.eigenstates.eigenvalues
        expected_transition_energies = [2.716, 5.476, 6.506, 7.986, 8.526, 9.111, 9.114, 9.606, 10.719, 11.041]
        for index, ref_energy in enumerate(expected_transition_energies):
            self.assertAlmostEqual(ref_energy, su.EV_PER_HARTREE * singlet_energies[index], places=2)

    def test_triplet_excited_states(self):
        self.excited_states_calculator.reference_calculation()
        self.excited_states_calculator.settings["number_eigenstates"] = 10
        self.excited_states_calculator.settings["spin_block"] = "triplet"
        results = self.excited_states_calculator.calculate()
        triplet_energies = results.excited_states.triplet.eigenstates.eigenvalues
        expected_transition_energies = [2.341, 4.316, 4.868, 5.706, 7.094, 7.908, 8.533, 9.037, 10.271, 10.369]
        for index, ref_energy in enumerate(expected_transition_energies):
            self.assertAlmostEqual(ref_energy, su.EV_PER_HARTREE * triplet_energies[index], places=2)

    def test_unrestricted_excited_states(self):
        self.excited_states_calculator.reference_calculator.settings["spin_mode"] = "unrestricted"
        self.excited_states_calculator.reference_calculation()
        self.excited_states_calculator.settings["number_eigenstates"] = 6
        results = self.excited_states_calculator.calculate()
        unrestricted_energies = results.excited_states.unrestricted.eigenstates.eigenvalues
        expected_transition_energies = [2.341, 2.716, 4.316, 4.868, 5.476, 5.706]
        for index, ref_energy in enumerate(expected_transition_energies):
            self.assertAlmostEqual(ref_energy, su.EV_PER_HARTREE * unrestricted_energies[index], places=2)

    def test_singlet_tddftb(self):
        self.tddftb.reference_calculation()
        self.tddftb.settings["number_eigenstates"] = 10
        results = self.tddftb.calculate()
        singlet_energies = results.excited_states.singlet.eigenstates.eigenvalues
        # These are the eigenvalue problem solution, i.e. the squares of the transition energies from DFTB+
        expected_transition_energies = [13.355, 13.485, 14.582, 15.581, 15.594, 16.894, 18.039, 18.660, 19.323, 20.513]
        for index, ref_energy in enumerate(expected_transition_energies):
            self.assertAlmostEqual(ref_energy, su.EV_PER_HARTREE * singlet_energies[index], places=2)

    def test_triplet_tddftb(self):
        self.tddftb.reference_calculation()
        self.tddftb.settings["number_eigenstates"] = 10
        self.tddftb.settings["spin_block"] = "triplet"
        results = self.tddftb.calculate()
        triplet_energies = results.excited_states.triplet.eigenstates.eigenvalues
        expected_transition_energies = [11.381, 13.168, 13.467, 14.001, 15.117, 15.479, 17.111, 17.876, 18.063, 19.585]
        for index, ref_energy in enumerate(expected_transition_energies):
            self.assertAlmostEqual(ref_energy, su.EV_PER_HARTREE * triplet_energies[index], places=2)

    def test_unrestricted_tddftb(self):
        self.tddftb.reference_calculator.settings["spin_mode"] = "unrestricted"
        self.tddftb.reference_calculation()
        self.tddftb.settings["number_eigenstates"] = 6
        results = self.tddftb.calculate()
        unrestricted_energies = results.excited_states.unrestricted.eigenstates.eigenvalues
        expected_transition_energies = [11.381, 13.168, 13.355, 13.467, 13.485, 14.001]
        for index, ref_energy in enumerate(expected_transition_energies):
            self.assertAlmostEqual(ref_energy, su.EV_PER_HARTREE * unrestricted_energies[index], places=2)

if __name__ == '__main__':
    unittest.main()
