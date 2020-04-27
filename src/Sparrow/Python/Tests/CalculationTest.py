import os
import unittest

import scine_sparrow


class TestSparrowPythonBindings(unittest.TestCase):
    def setUp(self):
        self.calculation = scine_sparrow.Calculation()
        # TODO: Maybe there is a better solution to make sure the parameters are found during the tests (i.e., when
        # the stuff is not yet installed)
        self.test_script_directory = os.path.dirname(os.path.realpath(__file__))
        parameter_directory = os.path.join(self.test_script_directory, '../../Sparrow/Resources')
        self.calculation.set_settings({'parameter_root': parameter_directory})
        self.calculation.set_elements(['H', 'H'])
        self.calculation.set_positions([[0, 0, 0], [1, 0, 0]])

    def test_energy(self):

        energy = self.calculation.calculate_energy()
        self.assertAlmostEqual(energy, -0.9962796281, places=6)

    def test_gradient(self):
        gradients = self.calculation.calculate_gradients()
        expected_gradients = [[-0.13828252490352153, 0.0, 0.0], [0.13828252490352153, 0.0, 0.0]]
        for i in range(0, 2):
            for j in range(0, 3):
                self.assertAlmostEqual(gradients[i][j], expected_gradients[i][j], places=6)

    def test_hessian(self):
        hessian = self.calculation.calculate_hessian()
        expected_hessian = [[1.40270507e-01, -5.55111512e-15, -5.55111512e-15,
                             -1.40270507e-01, -5.55111512e-15, -5.55111512e-15],
                            [0.00000000e+00,  7.31761956e-02,  0.00000000e+00,
                             0.00000000e+00, -7.31761956e-02,  0.00000000e+00],
                            [0.00000000e+00,  0.00000000e+00,  7.31761956e-02,
                             0.00000000e+00,  0.00000000e+00, -7.31761956e-02],
                            [-1.40270507e-01,  5.55111512e-15,  5.55111512e-15,
                             1.40270507e-01,  5.55111512e-15,  5.55111512e-15],
                            [0.00000000e+00, -7.31761956e-02,  0.00000000e+00,
                             0.00000000e+00,  7.31761956e-02,  0.00000000e+00],
                            [0.00000000e+00,  0.00000000e+00, -7.31761956e-02,
                             0.00000000e+00,  0.00000000e+00,  7.31761956e-02]]
        for i in range(len(expected_hessian)):
            for j in range(len(expected_hessian)):
                self.assertAlmostEqual(hessian[i][j], expected_hessian[i][j], places=6)

    def test_bond_orders(self):
        bond_orders = self.calculation.calculate_bond_orders().todense()
        expected_bond_orders = [[0, 1], [1, 0]]
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(bond_orders[i, j], expected_bond_orders[i][j], places=6)

    def test_settings(self):
        settings = self.calculation.get_settings()
        self.assertEqual(settings['molecular_charge'], '0')
        self.assertEqual(settings['spin_multiplicity'], '1')
        self.assertEqual(settings['unrestricted_calculation'], '0')
        self.assertEqual(settings['scf_mixer'], 'diis')
        self.assertEqual(settings['max_scf_iterations'], '100')
        self.assertEqual(settings['self_consistence_criterion'], '0.000010')

        self.calculation.set_settings({'molecular_charge': 1, 'spin_multiplicity': 2, 'unrestricted_calculation': True,
                                       'scf_mixer': 'ediis', 'max_scf_iterations': 200,
                                       'self_consistence_criterion': 0.1})
        settings = self.calculation.get_settings()
        self.assertEqual(settings['molecular_charge'], '1')
        self.assertEqual(settings['spin_multiplicity'], '2')
        self.assertEqual(settings['unrestricted_calculation'], '1')
        self.assertEqual(settings['scf_mixer'], 'ediis')
        self.assertEqual(settings['max_scf_iterations'], '200')
        self.assertEqual(settings['self_consistence_criterion'], '0.100000')

    def test_set_structure(self):
        file_path = os.path.join(self.test_script_directory, 'f2.xyz')
        self.calculation.set_structure(file_path)
        elements = self.calculation.get_elements()
        positions = self.calculation.get_positions()
        expected_elements = ['F', 'F']
        expected_positions = [[1, 2, 3], [4, 5, 6]]
        for i in range(2):
            self.assertEqual(elements[i], expected_elements[i])
            for j in range(3):
                self.assertAlmostEqual(positions[i][j], expected_positions[i][j] * 1.88972613, places=6)


if __name__ == '__main__':
    unittest.main()
