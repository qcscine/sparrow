import os
import subprocess
import unittest


class Calculation:
    def __init__(self, **settings):
        self.structure = os.path.join(os.path.dirname(os.path.realpath(__file__)), settings['structure'])
        self.method = settings['method']
        self.molecular_charge = settings['molecular_charge']
        self.spin_multiplicity = settings['spin_multiplicity']
        self.unrestricted = settings['unrestricted']
        self.calculate_gradients = settings['calculate_gradients']
        self.calculate_hessian = settings['calculate_hessian']
        self.energy = 0
        self.gradients = []
        self.first_hessian_element = 0
        self.hessian = []
        self._output = ''
        self._command_line = []

    def run(self):
        self.__prepare_input()
        self.__execute()
        self.__parse_output()

    def get_energy(self):
        return self.energy

    def get_gradients(self):
        return self.gradients

    def get_first_hessian_element(self):
        return self.first_hessian_element

    def get_hessian(self):
        return self.hessian

    def __prepare_input(self):
        self._command_line = ["sparrow", '--structure', self.structure, '--method', self.method,
                              '--molecular_charge', str(self.molecular_charge), '--spin_multiplicity',
                              str(self.spin_multiplicity), '--self_consistence_criterion', '1.e-9']

        if self.unrestricted:
            self._command_line.append('--unrestricted')

        if self.calculate_gradients:
            self._command_line.append('--gradients')

        if self.calculate_hessian:
            self._command_line.append('--hessian')
            self._command_line.append('--suppress_normal_modes')
            self._command_line.append('--output_to_file')

    def __execute(self):
        # If the return code is non-zero, this automatically raises an exception
        self._output = subprocess.check_output(self._command_line)

    def __parse_output(self):
        splitted_output = self._output.splitlines()
        for index, line in enumerate(splitted_output):
            if line == 'Energy [hartree]:':
                self.energy = float(splitted_output[index + 1])

            if line == 'Gradients [hartree/bohr]:':
                self.gradients = []
                counter = 2
                while True:
                    if splitted_output[index + counter].strip() == '':
                        break
                    splitted_line = splitted_output[index + counter].split()
                    gradient_x = float(splitted_line[1])
                    gradient_y = float(splitted_line[2])
                    gradient_z = float(splitted_line[3])
                    self.gradients.append([gradient_x, gradient_y, gradient_z])
                    counter += 1

            if line == 'Hessian [hartree/bohr^2]:':
                splitted_line = splitted_output[index + 2].split()
                self.first_hessian_element = float(splitted_line[1])

                self.hessian = []
                with open('hessian.dat') as fp:
                    for line in fp:
                        if line.startswith('#'):
                            continue
                        self.hessian.append([float(i) for i in line.split()])


class TestSparrowFast(unittest.TestCase):
    def test_energy_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -101.594627999, 1)

    def test_energy_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='AM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -101.235731061, 2)

    def test_energy_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -99.830098647, 2)

    def test_energy_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -94.616482772, 1)

    def test_energy_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -94.82507436, 2)

    def test_energy_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB0', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.63477, 2)

    def test_energy_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.71026, 2)

    def test_energy_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.57471, 2)

    # In order to convert the MOPAC gradients (given in kcal/mol/A) into atomic units, multiply them by 0.000843297.
    def test_gradients_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.011928714, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], 0.013450254, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], -0.012023743, 4)

    def test_gradients_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='AM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.004677906, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.010189013, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.007222255, 4)

    def test_gradients_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.000963553, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.00573894, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.002689228, 4)

    def test_gradients_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.004269321, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.008702266, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.006013156, 4)

    def test_gradients_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.005186126, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.005791405, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.004327844, 3)

    def test_gradients_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB0', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.01083, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.02137, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.01000, 3)

    def test_gradients_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.00325, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.01845, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.00801, 4)

    def test_gradients_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.01004, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.00976, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.00162, 3)

    # In order to convert the MOPAC Hessian into atomic units, remove their mass-weighting and then multiply them by
    # 0.0642.
    def test_hessian_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.710528114, 2)

        # Check that using a different number of CPUs yields the same Hessian
        parallel_hessian = calculation.get_hessian()
        os.environ["OMP_NUM_THREADS"] = "1"
        calculation.run()
        serial_hessian = calculation.get_hessian()
        for i in range(0, 84):
            for j in  range(0, 84):
                self.assertAlmostEqual(parallel_hessian[i][j], serial_hessian[i][j], places=5)

    def test_hessian_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB0', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.78975, 2)

    # TODO: Try to find a way to check that Sparrow fails due to missing parameters
    # def test_energy_sft_3_mndo_charge_0_multiplicity_1_unrestricted(self):
    #     calculation = Calculation(structure='sft_3.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
    #                               unrestricted=True, calculate_gradients=False)

    def test_relative_energy_differences_sft_5(self):
        mopac_energy_sft_5a = -51.394681808
        mopac_energy_sft_5b = -51.402473909

        calculation = Calculation(structure='sft_5a.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_5a = calculation.get_energy()

        calculation = Calculation(structure='sft_5b.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_5b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_5a, mopac_energy_sft_5a, 3)
        self.assertAlmostEqual(sparrow_energy_sft_5b, mopac_energy_sft_5b, 3)
        self.assertAlmostEqual(sparrow_energy_sft_5a - sparrow_energy_sft_5b,
                               mopac_energy_sft_5a - mopac_energy_sft_5b, 4)


class TestSparrowAll(TestSparrowFast):
    # In order to convert the MOPAC Hessian into atomic units, remove their mass-weighting and then multiply them by
    # 0.0642.
    def test_hessian_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.842550342, 2)

    def test_hessian_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='AM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.870422725, 3)

    def test_hessian_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.787611527, 3)

    def test_hessian_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.827046779, 3)

    def test_hessian_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.87420, 2)

    def test_hessian_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_hessian_element(), 0.80680, 2)

    def test_relative_energy_differences_sft_3(self):
        mopac_energy_sft_3a = -373.476133576
        mopac_energy_sft_3b = -373.243414082

        calculation = Calculation(structure='sft_3a.xyz', method='PM6', molecular_charge=0, spin_multiplicity=2,
                                  unrestricted=True, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_3a = calculation.get_energy()

        calculation = Calculation(structure='sft_3b.xyz', method='PM6', molecular_charge=0, spin_multiplicity=2,
                                  unrestricted=True, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_3b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_3a, mopac_energy_sft_3a, 2)
        self.assertAlmostEqual(sparrow_energy_sft_3b, mopac_energy_sft_3b, 2)
        self.assertAlmostEqual(sparrow_energy_sft_3a - sparrow_energy_sft_3b,
                               mopac_energy_sft_3a - mopac_energy_sft_3b, 4)

    def test_relative_energy_differences_sft_6_rm1(self):
        mopac_energy_sft_6a = -513.761856058
        mopac_energy_sft_6b = -513.764655824

        calculation = Calculation(structure='sft_6a.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_6a = calculation.get_energy()

        calculation = Calculation(structure='sft_6b.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_6b = calculation.get_energy()

        #self.assertAlmostEqual(sparrow_energy_sft_6a, mopac_energy_sft_6a, 3)
        #self.assertAlmostEqual(sparrow_energy_sft_6b, mopac_energy_sft_6b, 3)
        self.assertAlmostEqual(sparrow_energy_sft_6a - sparrow_energy_sft_6b,
                               mopac_energy_sft_6a - mopac_energy_sft_6b, 4)

    def test_relative_energy_differences_sft_6_mndo(self):
        mopac_energy_sft_6a = -519.135403913
        mopac_energy_sft_6b = -519.143483946

        calculation = Calculation(structure='sft_6a.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_6a = calculation.get_energy()

        calculation = Calculation(structure='sft_6b.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False)
        calculation.run()
        sparrow_energy_sft_6b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_6a, mopac_energy_sft_6a, 0)
        self.assertAlmostEqual(sparrow_energy_sft_6b, mopac_energy_sft_6b, 0)
        self.assertAlmostEqual(sparrow_energy_sft_6a - sparrow_energy_sft_6b,
                               mopac_energy_sft_6a - mopac_energy_sft_6b, 3)



if __name__ == '__main__':
    unittest.main()
