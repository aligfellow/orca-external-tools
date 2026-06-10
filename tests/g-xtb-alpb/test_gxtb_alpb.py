import unittest

from oet import ROOT_DIR
from oet.core.test_utilities import (
    OH,
    WATER,
    get_filenames,
    read_result_file,
    run_wrapper,
    write_input_file,
    write_xyz_file,
)

# Path to the scripts, adjust if needed.
gxtb_alpb_script_path = ROOT_DIR / "../../bin/oet_gxtb_alpb"
# Leave xtb_executable_path empty, if the (g-xTB-enabled) xtb from the system path should be called.
# NOTE: this must be the modified xtb build that supports `--gxtb` (github.com/grimme-lab/g-xTB),
# and the .gxtb/.eeq/.basisq parameter files must be discoverable (e.g. in $HOME or $GXTBHOME).
xtb_executable_path = ""
# Solvent passed to xtb's --alpb for the Delta-solv correction.
solvent = "thf"


def run_gxtb_alpb(inputfile: str, output_file: str) -> None:
    arguments = ["-s", solvent]
    if xtb_executable_path:
        arguments += ["--exe", xtb_executable_path]
    run_wrapper(
        inputfile=inputfile,
        script_path=gxtb_alpb_script_path,
        outfile=output_file,
        args=arguments,
    )


class GxtbAlpbTests(unittest.TestCase):
    def test_H2O_engrad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("H2O")

        write_xyz_file(xyz_file, WATER)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=0,
            multiplicity=1,
            ncores=2,
            do_gradient=1,
        )
        run_gxtb_alpb(input_file, output_file)
        expected_num_atoms = 3
        # TODO: regenerate the reference values below on a machine with the g-xTB-enabled xtb
        # build + parameter files (solvent="thf"), then remove this note.
        expected_energy = 0.0
        expected_gradients = [
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=9)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=9)

    def test_OH_anion_eng_grad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_anion")
        write_xyz_file(xyz_file, OH)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=-1,
            multiplicity=1,
            ncores=2,
            do_gradient=1,
        )
        run_gxtb_alpb(input_file, output_file)
        expected_num_atoms = 2
        # TODO: regenerate the reference values below on a machine with the g-xTB-enabled xtb
        # build + parameter files (solvent="thf"), then remove this note.
        expected_energy = 0.0
        expected_gradients = [
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=9)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=9)

    def test_OH_rad_eng_grad(self):
        xyz_file, input_file, engrad_out, output_file = get_filenames("OH_rad")
        write_xyz_file(xyz_file, OH)
        write_input_file(
            filename=input_file,
            xyz_filename=xyz_file,
            charge=0,
            multiplicity=2,
            ncores=2,
            do_gradient=1,
        )
        run_gxtb_alpb(input_file, output_file)
        expected_num_atoms = 2
        # TODO: regenerate the reference values below on a machine with the g-xTB-enabled xtb
        # build + parameter files (solvent="thf"), then remove this note.
        expected_energy = 0.0
        expected_gradients = [
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        ]

        try:
            num_atoms, energy, gradients = read_result_file(engrad_out)
        except Exception as e:
            raise FileNotFoundError(
                f"Error wrapper outputfile not found. Check {output_file} for details"
            ) from e

        self.assertEqual(num_atoms, expected_num_atoms)
        self.assertAlmostEqual(energy, expected_energy, places=7)
        for g1, g2 in zip(gradients, expected_gradients):
            self.assertAlmostEqual(g1, g2, places=7)


if __name__ == "__main__":
    unittest.main()
