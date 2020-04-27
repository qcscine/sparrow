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
        self.calculate_thermochemistry = settings['thermochemistry']
        self.wavefunction = settings['wavefunction']
        self.requested_temperature = 300.0
        self.energy = 0
        self.gradients = []
        self.hessian = []
        self.frequencies = []
        self.actual_temperature = 0
        self.enthalpy = 0
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

    def get_hessian(self):
        return self.hessian

    def get_frequencies(self):
        return self.frequencies

    def get_actual_temperature(self):
        return self.actual_temperature

    def get_enthalpy(self):
        return self.enthalpy

    def __prepare_input(self):
        self._command_line = ["sparrow", '--structure', self.structure, '--method', self.method,
                              '--molecular_charge', str(self.molecular_charge), '--spin_multiplicity',
                              str(self.spin_multiplicity), '--self_consistence_criterion', '1.e-9',
                              '--temperature', str(self.requested_temperature)]

        if self.unrestricted:
            self._command_line.append('--unrestricted')

        if self.calculate_gradients:
            self._command_line.append('--gradients')

        if self.calculate_hessian:
            self._command_line.append('--hessian')
            self._command_line.append('--output_to_file')

        if self.calculate_thermochemistry:
            self._command_line.append('--thermochemistry')

        if self.wavefunction:
            self._command_line.append('--wavefunction')

    def __execute(self):
        # If the return code is non-zero, this automatically raises an exception
        self._output = subprocess.check_output(self._command_line)

    def __parse_output(self):
        splitted_output = self._output.decode().splitlines()
        for index, line in enumerate(splitted_output):
            if 'Energy [hartree]:' in line:
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

            if line.startswith('Frequency:'):
                splitted_line = splitted_output[index].split()
                self.frequencies += [float(freq) for freq in splitted_line[1:]]

                self.hessian = []
                with open('hessian.dat') as fp:
                    for line in fp:
                        if line.startswith('#'):
                            continue
                        self.hessian.append([float(i) for i in line.split()])

            if line == 'Thermochemistry:':
                temperature_splitted_line = splitted_output[index + 1].split()
                overall_splitted_line = splitted_output[index + 8].split()
                self.actual_temperature = float(temperature_splitted_line[1])
                self.enthalpy = float(overall_splitted_line[2])


class TestSparrowFast(unittest.TestCase):
    def test_energy_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -101.594627999, 1)

    def test_energy_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='AM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -101.235731061, 2)

    def test_energy_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -99.830098647, 2)

    def test_energy_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -94.616482772, 1)

    def test_energy_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -94.82507436, 2)

    def test_energy_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB0', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.63477, 2)

    def test_energy_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.71026, 2)

    def test_energy_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.57471, 2)

    # In order to convert the MOPAC gradients (given in kcal/mol/A) into atomic units, multiply them by 0.000843297.
    def test_gradients_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.011928714, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], 0.013450254, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], -0.012023743, 4)

    def test_gradients_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='AM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.004677906, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.010189013, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.007222255, 4)

    def test_gradients_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.000963553, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.00573894, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.002689228, 4)

    def test_gradients_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.004269321, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.008702266, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.006013156, 4)

    def test_gradients_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.005186126, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.005791405, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.004327844, 3)

    def test_gradients_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB0', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.01083, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.02137, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.01000, 3)

    def test_gradients_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.00325, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.01845, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.00801, 4)

    def test_gradients_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.01004, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.00976, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.00162, 3)

    # In order to convert the MOPAC Hessian into atomic units, remove their mass-weighting and then multiply them by
    # 0.0642.
    def test_hessian_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        # Reference value from MOPAC 2016 is -122.7; actual value from Sparrow should be -125.522
        self.assertAlmostEqual(calculation.get_frequencies()[0], -125.522, places=2)

        # Check that using a different number of CPUs yields the same Hessian
        parallel_hessian = calculation.get_hessian()
        os.environ["OMP_NUM_THREADS"] = "1"
        calculation.run()
        serial_hessian = calculation.get_hessian()
        for i in range(0, 84):
            for j in range(0, 84):
                self.assertAlmostEqual(parallel_hessian[i][j], serial_hessian[i][j], places=5)

    def test_hessian_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB0', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        # Reference value from DFTB+ 18.2 is -136.73; actual value from Sparrow should be -136.824
        self.assertAlmostEqual(calculation.get_frequencies()[0], -136.824, 2)

    # TODO: Try to find a way to check that Sparrow fails due to missing parameters
    # def test_energy_sft_3_mndo_charge_0_multiplicity_1_unrestricted(self):
    #     calculation = Calculation(structure='sft_3.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
    #                               unrestricted=True, calculate_gradients=False)

    def test_relative_energy_differences_sft_5(self):
        mopac_energy_sft_5a = -51.394681808
        mopac_energy_sft_5b = -51.402473909

        calculation = Calculation(structure='sft_5a.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        sparrow_energy_sft_5a = calculation.get_energy()

        calculation = Calculation(structure='sft_5b.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=True, thermochemistry=True)
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
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        # Reference value from MOPAC 2016 is -316.4; actual value from Sparrow should be -320.529
        self.assertAlmostEqual(calculation.get_frequencies()[0], -320.529, places=2)

    def test_hessian_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='AM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        # Reference value from MOPAC 2016 is -406.5; actual value from Sparrow should be -410.378
        self.assertAlmostEqual(calculation.get_frequencies()[0], -410.378, places=2)

    def test_hessian_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        self.assertAlmostEqual(calculation.get_frequencies()[0], -112.5, places=1)

    def test_hessian_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        # Reference value from MOPAC 2016 is -202.3; actual value from Sparrow should be -202.165
        self.assertAlmostEqual(calculation.get_frequencies()[0], -202.165, places=2)

    def test_hessian_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()

        # Reference value from DFTB+ 18.2 is -165.02
        self.assertAlmostEqual(calculation.get_frequencies()[0], -158.752, 2)

    def test_hessian_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=True,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        # Reference value from DFTB+ 18.2 is -71.33; actual value from Sparrow should be -69.3909
        self.assertAlmostEqual(calculation.get_frequencies()[0], -69.3909, 2)

    def test_relative_energy_differences_sft_3(self):
        mopac_energy_sft_3a = -373.476133576
        mopac_energy_sft_3b = -373.243414082

        calculation = Calculation(structure='sft_3a.xyz', method='PM6', molecular_charge=0, spin_multiplicity=2,
                                  unrestricted=True, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        sparrow_energy_sft_3a = calculation.get_energy()

        calculation = Calculation(structure='sft_3b.xyz', method='PM6', molecular_charge=0, spin_multiplicity=2,
                                  unrestricted=True, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
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
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        sparrow_energy_sft_6a = calculation.get_energy()

        calculation = Calculation(structure='sft_6b.xyz', method='RM1', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        sparrow_energy_sft_6b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_6a, mopac_energy_sft_6a, places=2)
        self.assertAlmostEqual(sparrow_energy_sft_6b, mopac_energy_sft_6b, places=2)
        self.assertAlmostEqual(sparrow_energy_sft_6a - sparrow_energy_sft_6b,
                               mopac_energy_sft_6a - mopac_energy_sft_6b, places=4)

    def test_relative_energy_differences_sft_6_mndo(self):
        mopac_energy_sft_6a = -519.135403913
        mopac_energy_sft_6b = -519.143483946

        calculation = Calculation(structure='sft_6a.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        sparrow_energy_sft_6a = calculation.get_energy()

        calculation = Calculation(structure='sft_6b.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False)
        calculation.run()
        sparrow_energy_sft_6b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_6a, mopac_energy_sft_6a, 0)
        self.assertAlmostEqual(sparrow_energy_sft_6b, mopac_energy_sft_6b, 0)
        self.assertAlmostEqual(sparrow_energy_sft_6a - sparrow_energy_sft_6b,
                               mopac_energy_sft_6a - mopac_energy_sft_6b, 3)

    def test_molden_input_file_generation_sft_1(self):
        calculation = Calculation(structure='sft_1.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=True, thermochemistry=False)
        calculation.run()
        with open('wavefunction.molden.input') as mol_output, open('molden_input_reference.molden') as mol_test:
            in_mo_block = False
            for line1, line2 in zip(mol_output, mol_test):
                if 'Sym=' in line1:
                    in_mo_block = False
                if 'Ene=' in line1:
                    energy1 = float(line1.split()[1])
                    energy2 = float(line1.split()[1])
                    self.assertAlmostEqual(energy1, energy2, 6)
                elif in_mo_block:
                    mo1 = float(line1.split()[1])
                    mo2 = float(line2.split()[1])
                    self.assertAlmostEqual(abs(mo1), abs(mo2), 6)
                else:
                    self.assertTrue(line1 == line2)
                if 'Occup=' in line1:
                    in_mo_block = True

    def test_thermochemistry_sft_1(self):
        calculation = Calculation(structure='sft_1_opt.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=True, thermochemistry=True)
        calculation.run()
        self.assertEqual(calculation.get_actual_temperature(), 300.0)
        self.assertAlmostEqual(calculation.get_enthalpy(), -248928.5, 1)


if __name__ == '__main__':
    unittest.main()
