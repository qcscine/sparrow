__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import unittest
import scine_utilities as su
import scine_sparrow


class TestCalculatorPythonBindings(unittest.TestCase):
    def setUp(self):
        structure = su.AtomCollection(2)
        structure.elements = [su.ElementType.H, su.ElementType.H]
        structure.positions = [[-0.7, 0, 0], [0.7, 0, 0]]
        manager = su.core.ModuleManager.get_instance()
        self.calculator = manager.get('calculator', 'PM6')
        self.calculator.structure = structure

    def test_energy(self):
        self.calculator.set_required_properties([su.Property.Energy])
        results = self.calculator.calculate()
        self.assertAlmostEqual(results.energy, -1.0331782672319838, places=6)

    def test_gradient(self):
        self.calculator.set_required_properties([su.Property.Gradients])
        results = self.calculator.calculate()
        expected_gradients = [[0.02468439607921747, 0.0, 0.0], [-0.02468439607921747, 0.0, 0.0]]
        for i in range(0, 2):
            for j in range(0, 3):
                self.assertAlmostEqual(results.gradients[i][j], expected_gradients[i][j], places=6)

    def test_hessian(self):
        self.calculator.set_required_properties([su.Property.Hessian])
        results = self.calculator.calculate()
        expected_hessian = [
            [ 0.60221448,  0.00000000,  0.00000000, -0.60221448,  0.00000000,  0.00000000],
            [ 0.00000000, -0.01761590,  0.00000000,  0.00000000,  0.01761590,  0.00000000],
            [ 0.00000000,  0.00000000, -0.01761590,  0.00000000,  0.00000000,  0.01761590],
            [-0.60221448,  0.00000000,  0.00000000,  0.60221448,  0.00000000,  0.00000000],
            [ 0.00000000,  0.01761590,  0.00000000,  0.00000000, -0.01761590,  0.00000000],
            [ 0.00000000,  0.00000000,  0.01761590,  0.00000000,  0.00000000, -0.01761590]
        ]
        for i in range(len(expected_hessian)):
            for j in range(len(expected_hessian)):
                self.assertAlmostEqual(results.hessian[i][j], expected_hessian[i][j], places=6)

    def test_atomic_hessians(self):
        self.calculator.set_required_properties([su.Property.AtomicHessians])
        results = self.calculator.calculate()
        calculated_atomic_hessians = results.atomic_hessian.get_atomic_hessians()
        # Reference data has been generated with original implementation by A. Vaucher
        expected_atomic_hessians = [
                                    [[0.60216,  0.00000000,  0.00000000],
                                     [0.0000000,  -0.0176317,   0.00000000],
                                     [0.00000000,  0.0000000,   -0.0176317]],
                                    [[0.60216,  0.00000000,  0.00000000],
                                     [0.0000000,  -0.0176317,   0.00000000],
                                     [0.00000000,  0.0000000,   -0.0176317]],
                                    ]
        for i in range(len(expected_atomic_hessians[0])):
            for j in range(len(expected_atomic_hessians[0])):
                self.assertAlmostEqual(calculated_atomic_hessians[0][i][j], expected_atomic_hessians[0][i][j], places=5)
                self.assertAlmostEqual(calculated_atomic_hessians[1][i][j], expected_atomic_hessians[1][i][j], places=5)

    def test_bond_orders(self):
        self.calculator.set_required_properties([su.Property.BondOrderMatrix])
        results = self.calculator.calculate()
        bond_orders = results.bond_orders.matrix.todense()
        expected_bond_orders = [[0, 1], [1, 0]]
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(bond_orders[i, j], expected_bond_orders[i][j], places=6)


if __name__ == '__main__':
    unittest.main()
