__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import os
import subprocess
import unittest
import numpy as np


class Calculation:
    def __init__(self, **settings):
        self.structure = os.path.join(os.path.dirname(os.path.realpath(__file__)), settings['structure'])
        self.method = settings['method']
        self.molecular_charge = settings['molecular_charge']
        self.spin_multiplicity = settings['spin_multiplicity']
        self.unrestricted = settings['unrestricted']
        self.calculate_gradients = settings['calculate_gradients']
        self.calculate_hessian = settings['calculate_hessian']
        self.calculate_atomic_hessians = True if ('calculate_atomic_hessians' in settings and settings['calculate_atomic_hessians']) else False
        self.calculate_thermochemistry = settings['thermochemistry']
        self.requested_symmetry_number = settings['symmetry_number']
        self.wavefunction = settings['wavefunction']
        self.excited_states = True if ('excited_states' in settings and settings['excited_states']) else False
        self.orbital_mixing = settings['number_orbital_mixes'] if ('number_orbital_mixes' in settings) else 0
        self.spin_block = 'singlet' if 'spin_block' not in settings else settings['spin_block']
        self.requested_temperature = 300.0
        self.symmetry_number = 1
        self.energy = 0
        self.first_singlet_energy = 0
        self.first_triplet_energy = 0
        self.first_unrestricted_energy = 0
        self.gradients = []
        self.hessian = []
        self.atomic_hessians = []
        self.frequencies = []
        self.actual_temperature = 0
        self.enthalpy = 0
        self.mixer = settings['scf_mixer']
        self._output = ''
        self._command_line = []

    def run(self):
        self.__prepare_input()
        self.__execute()
        self.__parse_output()

    def get_energy(self):
        return self.energy

    def get_first_singlet_transition_energy(self):
        return self.first_singlet_energy

    def get_first_triplet_transition_energy(self):
        return self.first_triplet_energy

    def get_first_unrestricted_transition_energy(self):
        return self.first_unrestricted_energy

    def get_gradients(self):
        return self.gradients

    def get_hessian(self):
        return self.hessian

    def get_frequencies(self):
        return self.frequencies

    def get_actual_temperature(self):
        return self.actual_temperature

    def get_actual_symmetry_number(self):
        return self.symmetry_number

    def get_enthalpy(self):
        return self.enthalpy

    def get_atomic_hessians(self):
        return self.atomic_hessians

    def __prepare_input(self):
        self._command_line = [
            "sparrow",
            '--structure', self.structure,
            '--method', self.method,
            '--scf_mixer', self.mixer,
            '--molecular_charge', str(self.molecular_charge),
            '--spin_multiplicity', str(self.spin_multiplicity),
            '--self_consistence_criterion', '1.e-9',
            '--density_rmsd_criterion', '1.e-9',
            '--temperature', str(self.requested_temperature),
            "--symmetry_number", str(self.requested_symmetry_number),
            "--spin_block", self.spin_block,
            "--number_orbital_mixes", str(self.orbital_mixing)
            ]

        if self.unrestricted:
            self._command_line.append('--unrestricted')

        if self.calculate_gradients:
            self._command_line.append('--gradients')

        if self.calculate_hessian:
            self._command_line.append('--hessian')
            self._command_line.append('--output_to_file')

        if self.calculate_atomic_hessians:
            self._command_line.append('--atomic_hessians')

        if self.calculate_thermochemistry:
            self._command_line.append('--thermochemistry')

        if self.wavefunction:
            self._command_line.append('--wavefunction')

        if self.excited_states:
            self._command_line.append('--excited_states')

    def __execute(self):
        # If the return code is non-zero this automatically raises an exception
        self._output = subprocess.check_output(self._command_line)

    def __parse_output(self):
        splitted_output = self._output.decode().splitlines()
        for index, line in enumerate(splitted_output):
            if 'Energy [hartree]:' in line:
                self.energy = float(splitted_output[index + 1])

            if 'Temperature:' in line:
                self.actual_temperature = float(splitted_output[index].split()[1])

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
                        if line.startswith('#') or not line.strip():
                            continue
                        self.hessian.append([float(i) for i in line.split()])

            if line == 'Thermochemistry:':
                overall_splitted_line = splitted_output[index + 8].split()
                self.enthalpy = float(overall_splitted_line[2])
                symmetry_splitted_line = splitted_output[index + 1].split()
                self.symmetry_number = float(symmetry_splitted_line[3])

            if line == 'Atomic Hessians [hartree/bohr^2]:':
                self.atomic_hessians = []
                counter = 1
                while True:
                    if splitted_output[index + counter].startswith('='):
                        break
                    if splitted_output[index+counter].split()[0] == 'Atomic':
                        counter += 2
                    splitted_line = splitted_output[index + counter].split()
                    xx = float(splitted_line[1])
                    xy = float(splitted_line[2])
                    xz = float(splitted_line[3])

                    counter += 1
                    splitted_line = splitted_output[index + counter].split()
                    yx = float(splitted_line[1])
                    yy = float(splitted_line[2])
                    yz = float(splitted_line[3])

                    counter += 1
                    splitted_line = splitted_output[index + counter].split()
                    zx = float(splitted_line[1])
                    zy = float(splitted_line[2])
                    zz = float(splitted_line[3])

                    self.atomic_hessians.append([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])
                    counter += 2

            if line.startswith('Singlet electronic transitions'):
                self.first_singlet_energy = float(splitted_output[index + 3].split()[6])

            if line.startswith('Triplet electronic transitions'):
                self.first_triplet_energy = float(splitted_output[index + 3].split()[6])

            if line.startswith('Unrestricted electronic transitions'):
                self.first_unrestricted_energy = float(splitted_output[index + 3].split()[6])

class TestSparrowFast(unittest.TestCase):

    def test_energy_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='MNDO', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -101.594627999, 1)

    def test_energy_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='AM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -101.235731061, 2)

    def test_energy_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='RM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -99.830098647, 2)

    def test_energy_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -94.616482772, 1)

    def test_energy_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -94.82507436, 2)

    def test_energy_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB0', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.63477, 2)

    def test_energy_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB2', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.71026, 2)

    def test_energy_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_energy(), -35.57471, 2)

    # In order to convert the MOPAC gradients (given in kcal/mol/A) into atomic units, multiply them by 0.000843297.
    def test_gradients_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='MNDO', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.011928714, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], 0.013450254, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], -0.012023743, 4)

    def test_gradients_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='AM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.004677906, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.010189013, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.007222255, 4)

    def test_gradients_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='RM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.000963553, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.00573894, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.002689228, 4)

    def test_gradients_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.004269321, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.008702266, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.006013156, 4)

    def test_gradients_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], -0.005186126, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.005791405, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.004327844, 3)

    def test_gradients_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB0', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.01083, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.02137, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.01000, 3)

    def test_gradients_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB2', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.00325, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.01845, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.00801, 4)

    def test_gradients_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=True, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_gradients()[0][0], 0.01004, 4)
        self.assertAlmostEqual(calculation.get_gradients()[0][1], -0.00976, 3)
        self.assertAlmostEqual(calculation.get_gradients()[0][2], 0.00162, 3)

    def test_hessian_ethane_opt_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='ethane_opt.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from MOPAC 2016 is 2787.80; actual value from Sparrow should be 2788.55
        self.assertAlmostEqual(calculation.get_frequencies()[11], 2788.55, places=2)

        # Check that using a different number of CPUs yields the same Hessian
        parallel_hessian = calculation.get_hessian()
        os.environ["OMP_NUM_THREADS"] = "1"
        calculation.run()
        serial_hessian = calculation.get_hessian()
        for i in range(0, 18):
            for j in range(0, 18):
                self.assertAlmostEqual(parallel_hessian[i][j], serial_hessian[i][j], places=5)

    def test_steering_stretched_metane_pm6_charge_0_multiplicity_1_unrestricted(self):
        calculation_without_steering = Calculation(structure='stretched_methane.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=True, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', symmetry_number=1, number_orbital_mixes=0)
        calculation = Calculation(structure='stretched_methane.xyz', method='PM6', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=True, calculate_gradients=True, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', symmetry_number=1, number_orbital_mixes=3)
        calculation_without_steering.run()
        calculation.run()
        self.assertLessEqual(calculation.get_energy(), calculation_without_steering.get_energy())

    def test_hessian_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB0', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from DFTB+ 18.2 is -136.73; actual value from Sparrow should be -136.615
        self.assertAlmostEqual(calculation.get_frequencies()[0], -136.615, 2)

    # TODO: Try to find a way to check that Sparrow fails due to missing parameters
    # def test_energy_sft_3_mndo_charge_0_multiplicity_1_unrestricted(self):
    #     calculation = Calculation(structure='sft_3.xyz', method='MNDO', molecular_charge=0, spin_multiplicity=1,
    #                               unrestricted=True, calculate_gradients=False)

    def test_relative_energy_differences_sft_5(self):
        mopac_energy_sft_5a = -51.394681808
        mopac_energy_sft_5b = -51.402473909

        calculation = Calculation(
            structure='sft_5a.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_5a = calculation.get_energy()

        calculation = Calculation(
            structure='sft_5b.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=True,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_5b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_5a, mopac_energy_sft_5a, 3)
        self.assertAlmostEqual(sparrow_energy_sft_5b, mopac_energy_sft_5b, 3)
        self.assertAlmostEqual(sparrow_energy_sft_5a - sparrow_energy_sft_5b,
                               mopac_energy_sft_5a - mopac_energy_sft_5b, 4)

    def test_atomic_hessians_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='MNDO', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original re-implementation in Sparrow
        self.assertAlmostEqual(atomic_hessians[0][0][0], 1.42004, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][1], 0.00703889, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.157166, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][0], 0.00703889, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 1.2543, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.325494, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.157166, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.325494, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.886366, 5)

    def test_atomic_hessians_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='AM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original re-implementation in Sparrow
        self.assertAlmostEqual(atomic_hessians[0][0][0], 1.42634, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][1], 0.00326763, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.158518, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][0], 0.00326763, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 1.25643, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.326845, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.158518, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.326845, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.888875, 5)

    def test_atomic_hessians_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='RM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original re-implementation in Sparrow
        self.assertAlmostEqual(atomic_hessians[0][0][0], 1.34508, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][1], 0.01476, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.148802, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][0], 0.01476, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 1.19146, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.298299, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.148802, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.298299, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.858051, 5)

    def test_atomic_hessians_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original re-implementation in Sparrow
        self.assertAlmostEqual(atomic_hessians[0][0][0], 1.33617, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][1], 0.00174254, 5)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.134834, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][0], 0.00174254, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 1.16866, 5)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.271352, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.134834, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.271352, 5)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.86339, 5)

    def test_atomic_hessians_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original implementation by A. Vaucher
        self.assertAlmostEqual(atomic_hessians[0][0][0], 1.28073, 4)
        self.assertAlmostEqual(atomic_hessians[0][0][1], -0.0118772, 4)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.134253, 4)
        self.assertAlmostEqual(atomic_hessians[0][1][0], -0.0118772, 4)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 1.08065, 4)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.269756, 4)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.134253, 4)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.269756, 4)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.771995, 4)

    def test_atomic_hessians_sft_1_dftb0_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB0', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original implementation by A. Vaucher
        self.assertAlmostEqual(atomic_hessians[0][0][0], 0.807415, 2)
        self.assertAlmostEqual(atomic_hessians[0][0][1], -0.0374171, 3)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.0638167, 3)
        self.assertAlmostEqual(atomic_hessians[0][1][0], -0.0374171, 3)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 0.599173, 2)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.103017, 3)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.0638167, 3)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.103017, 3)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.476781, 2)

    def test_atomic_hessians_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB2', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original implementation by A. Vaucher
        self.assertAlmostEqual(atomic_hessians[0][0][0], 0.890185, 4)
        self.assertAlmostEqual(atomic_hessians[0][0][1], -0.0562036, 4)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.0637091, 4)
        self.assertAlmostEqual(atomic_hessians[0][1][0], -0.0562036, 4)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 0.622531, 4)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.0982637, 4)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.0637091, 4)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.0982637, 4)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.500765, 4)

    def test_atomic_hessians_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=True, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        atomic_hessians = calculation.get_atomic_hessians()
        self.assertEqual(len(atomic_hessians), 28)
        # Reference data has been generated with original implementation by A. Vaucher
        self.assertAlmostEqual(atomic_hessians[0][0][0], 0.818021, 2)
        self.assertAlmostEqual(atomic_hessians[0][0][1], -0.0396974, 3)
        self.assertAlmostEqual(atomic_hessians[0][0][2], -0.0623885, 2)
        self.assertAlmostEqual(atomic_hessians[0][1][0], -0.0396974, 3)
        self.assertAlmostEqual(atomic_hessians[0][1][1], 0.594924, 2)
        self.assertAlmostEqual(atomic_hessians[0][1][2], -0.0978562, 3)
        self.assertAlmostEqual(atomic_hessians[0][2][0], -0.0623885, 2)
        self.assertAlmostEqual(atomic_hessians[0][2][1], -0.0978562, 3)
        self.assertAlmostEqual(atomic_hessians[0][2][2], 0.478853, 2)

    @classmethod
    def tearDownClass(cls):
        if os.path.isfile('energy.dat'):
            os.remove('energy.dat')
        if os.path.isfile('gradients.dat'):
            os.remove('gradients.dat')
        if os.path.isfile('hessian.dat'):
            os.remove('hessian.dat')
        if os.path.isfile('wavefunction.molden.input'):
            os.remove('wavefunction.molden.input')


class TestSparrowSlow(unittest.TestCase):

    # In order to convert the MOPAC Hessian into atomic units, remove their mass-weighting and then multiply them by
    # 0.0642.
    def test_hessian_sft_1_mndo_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='MNDO', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from MOPAC 2016 is -316.4; actual value from Sparrow should be -320.713
        self.assertAlmostEqual(calculation.get_frequencies()[0], -320.713, places=2)

    def test_hessian_sft_1_am1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='AM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from MOPAC 2016 is -406.5; actual value from Sparrow should be -410.523
        self.assertAlmostEqual(calculation.get_frequencies()[0], -410.523, places=2)

    def test_hessian_sft_1_rm1_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='RM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_frequencies()[0], -112.562, places=1)

    def test_hessian_sft_1_pm3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from MOPAC 2016 is -202.3; actual value from Sparrow should be -202.18
        self.assertAlmostEqual(calculation.get_frequencies()[0], -202.18, places=2)

    def test_hessian_sft_1_dftb2_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB2', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()

        # Reference value from DFTB+ 18.2 is -165.02
        self.assertAlmostEqual(calculation.get_frequencies()[0], -158.922, 2)

    def test_hessian_sft_1_dftb3_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='DFTB3', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from DFTB+ 18.2 is -71.33; actual value from Sparrow should be -69.2154
        self.assertAlmostEqual(calculation.get_frequencies()[0], -69.2154, 2)

    # In order to convert the MOPAC Hessian into atomic units, remove their mass-weighting and then multiply them by
    # 0.0642.
    def test_hessian_sft_1_pm6_charge_0_multiplicity_1_restricted(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=True,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        # Reference value from MOPAC 2016 is -122.7; actual value from Sparrow should be -125.561
        self.assertAlmostEqual(calculation.get_frequencies()[0], -125.561, places=2)

        # Check that using a different number of CPUs yields the same Hessian
        parallel_hessian = calculation.get_hessian()
        os.environ["OMP_NUM_THREADS"] = "1"
        calculation.run()
        serial_hessian = calculation.get_hessian()
        for i in range(0, 84):
            for j in range(0, 84):
                self.assertAlmostEqual(parallel_hessian[i][j], serial_hessian[i][j], places=5)

    def test_relative_energy_differences_sft_3(self):
        mopac_energy_sft_3a = -373.476133576
        mopac_energy_sft_3b = -373.243414082

        calculation = Calculation(
            structure='sft_3a.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=2, unrestricted=True,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_3a = calculation.get_energy()

        calculation = Calculation(
            structure='sft_3b.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=2, unrestricted=True,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_3b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_3a, mopac_energy_sft_3a, 2)
        self.assertAlmostEqual(sparrow_energy_sft_3b, mopac_energy_sft_3b, 2)
        self.assertAlmostEqual(sparrow_energy_sft_3a - sparrow_energy_sft_3b,
                               mopac_energy_sft_3a - mopac_energy_sft_3b, 4)

    def test_relative_energy_differences_sft_6_rm1(self):
        mopac_energy_sft_6a = -513.761856058
        mopac_energy_sft_6b = -513.764655824

        calculation = Calculation(
            structure='sft_6a.xyz', method='RM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_6a = calculation.get_energy()

        calculation = Calculation(
            structure='sft_6b.xyz', method='RM1', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_6b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_6a, mopac_energy_sft_6a, places=2)
        self.assertAlmostEqual(sparrow_energy_sft_6b, mopac_energy_sft_6b, places=2)
        self.assertAlmostEqual(sparrow_energy_sft_6a - sparrow_energy_sft_6b,
                               mopac_energy_sft_6a - mopac_energy_sft_6b, places=4)

    def test_relative_energy_differences_sft_6_mndo(self):
        mopac_energy_sft_6a = -519.135403913
        mopac_energy_sft_6b = -519.143483946

        calculation = Calculation(
            structure='sft_6a.xyz', method='MNDO', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_6a = calculation.get_energy()

        calculation = Calculation(
            structure='sft_6b.xyz', method='MNDO', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        sparrow_energy_sft_6b = calculation.get_energy()

        self.assertAlmostEqual(sparrow_energy_sft_6a, mopac_energy_sft_6a, 0)
        self.assertAlmostEqual(sparrow_energy_sft_6b, mopac_energy_sft_6b, 0)
        self.assertAlmostEqual(sparrow_energy_sft_6a - sparrow_energy_sft_6b,
                               mopac_energy_sft_6a - mopac_energy_sft_6b, 3)

    def test_molden_input_file_generation_sft_1(self):
        calculation = Calculation(
            structure='sft_1.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=True,
            thermochemistry=False, scf_mixer='diis', symmetry_number=1)
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
                    normalizedMO1 = mo1 / np.linalg.norm(mo1)
                    normalizedMO2 = mo2 / np.linalg.norm(mo2)
                    # Check that the overlap between the two coefficients vectors is equal to 1.
                    self.assertAlmostEqual(np.abs(np.dot(normalizedMO1, normalizedMO2)), 1, 6)
                else:
                    self.assertTrue(line1 == line2)
                if 'Occup=' in line1:
                    in_mo_block = True

    def test_thermochemistry_sft_1(self):
        calculation = Calculation(
            structure='sft_1_opt.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=True, scf_mixer='diis', symmetry_number=1)
        calculation.run()
        conversionFactor = 4.359744650e-18 / 1000 * 6.022140857e23
        self.assertAlmostEqual(calculation.get_enthalpy(), 41.9 + calculation.get_energy() * conversionFactor, 1)
        self.assertEqual(calculation.get_actual_temperature(), 300.0)

    def test_symmetry_factor_sft_1(self):
        calculation = Calculation(
            structure='sft_1_opt.xyz', method='PM6', molecular_charge=0,
            spin_multiplicity=1, unrestricted=False,
            calculate_gradients=False, calculate_hessian=False,
            calculate_atomic_hessians=False, wavefunction=False,
            thermochemistry=True, scf_mixer='diis', symmetry_number=2)
        calculation.run()
        self.assertEqual(calculation.get_actual_symmetry_number(), 2)
        self.assertEqual(calculation.get_actual_temperature(), 300.0)

    def test_dftb2_first_singlet_excitation_energy_beta_carotene(self):
        calculation = Calculation(structure='beta_carotene_opt.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', excited_states=True, symmetry_number=1)
        calculation.run()
        self.assertAlmostEqual(calculation.get_first_singlet_transition_energy(), 1.807, 2)

    def test_dftb2_first_excitation_energies_sft_1(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', excited_states=True, symmetry_number=1,
                                  spin_block='both')
        calculation.run()
        # Reference transition energies calculated with DFTB+.
        self.assertAlmostEqual(calculation.get_first_triplet_transition_energy(), 4.086, 2)
        self.assertAlmostEqual(calculation.get_first_singlet_transition_energy(), 4.343, 2)

    def test_dftb3_first_excitation_energies_sft_1(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', excited_states=True, symmetry_number=1,
                                  spin_block='singlet')
        calculation.run()
        # Reference transition energies calculated with DFTB+.
        self.assertAlmostEqual(calculation.get_first_singlet_transition_energy(), 4.492, 2)

    def test_dftb2_first_excitation_energy_unrestricted_sft_1(self):
        calculation = Calculation(structure='sft_1.xyz', method='DFTB2', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=True, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', excited_states=True, symmetry_number=1)
        calculation.run()
        # Reference transition energy taken as the lowest energy from test before, as ground state is closed-shell.
        self.assertAlmostEqual(calculation.get_first_unrestricted_transition_energy(), 4.086, 2)

    def test_PM3_first_excitation_energies_formaldehyde(self):
        calculation = Calculation(structure='formaldehyde.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=False, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', excited_states=True, symmetry_number=1,
                                  spin_block='both')
        calculation.run()
        # Reference transition energies calculated with Orca 4.2.0
        self.assertAlmostEqual(calculation.get_first_triplet_transition_energy(), 2.341, 2)
        self.assertAlmostEqual(calculation.get_first_singlet_transition_energy(), 2.716, 2)

    def test_PM3_first_excitation_energy_unrestricted_formaldehyde(self):
        calculation = Calculation(structure='formaldehyde.xyz', method='PM3', molecular_charge=0, spin_multiplicity=1,
                                  unrestricted=True, calculate_gradients=False, calculate_hessian=False,
                                  wavefunction=False, thermochemistry=False, scf_mixer='diis', excited_states=True, symmetry_number=1)
        calculation.run()
        # Reference transition energy taken as the lowest energy from test before, as ground state is closed-shell.
        self.assertAlmostEqual(calculation.get_first_unrestricted_transition_energy(), 2.341, 2)

    @classmethod
    def tearDownClass(cls):
        if os.path.isfile('energy.dat'):
            os.remove('energy.dat')
        if os.path.isfile('gradients.dat'):
            os.remove('gradients.dat')
        if os.path.isfile('hessian.dat'):
            os.remove('hessian.dat')
        if os.path.isfile('wavefunction.molden.input'):
            os.remove('wavefunction.molden.input')

if __name__ == '__main__':
    unittest.main()
