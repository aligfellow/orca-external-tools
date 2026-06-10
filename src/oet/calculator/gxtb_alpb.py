#!/usr/bin/env python3
"""
Wrapper that adds an ALPB(solvent) implicit-solvation correction (computed at the GFN2-xTB level)
on top of a g-xTB gas-phase energy and gradient, compatible with ORCA's ExtTool interface.

E_total = E(g-xTB, gas) + [E(GFN2, ALPB) - E(GFN2, gas)]
grad    = grad(g-xTB, gas) + [grad(GFN2, ALPB) - grad(GFN2, gas)]

Uses the modified xtb binary (github.com/grimme-lab/g-xTB), which supports both --gxtb
and standard GFN2-xTB. Three xtb single-points are performed per ORCA gradient call.

Provides
--------
class: GxtbAlpbCalc(BaseCalc)
    Class for performing a g-xTB calculation with a GFN2-ALPB Δsolv correction
main: function
    Main function
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import Any

from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import (
    check_path,
    mult_to_nue,
    nat_from_xyzfile,
    run_command,
)


class GxtbAlpbCalc(BaseCalc):
    @property
    def PROGRAM_NAMES(self) -> list[str]:
        """Program names to search for in PATH"""
        return ["xtb", "otool_xtb"]

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add gxtb_alpb parsing options.

        Parameters
        ----------
        parser: ArgumentParser
            Parser that should be extended
        """
        parser.add_argument("-e", "--exe", dest="prog", help="Path to the xtb executable")
        parser.add_argument(
            "-s",
            "--solvent",
            dest="solvent",
            default="thf",
            help="Solvent name passed to xtb's --alpb (default: thf)",
        )

    def read_xtbout(
        self,
        output_file: str | Path,
        grad_file: str | Path,
        natoms: int,
        dograd: bool,
    ) -> tuple[float, list[float]]:
        """
        Read energy and (optionally) gradient from one xtb run.

        Parameters
        ----------
        output_file: str | Path
            xtb stdout file containing the 'TOTAL ENERGY' line
        grad_file: str | Path
            Turbomole-style .gradient file (only read if dograd)
        natoms: int
            Number of atoms in the system
        dograd: bool
            Whether to read the gradient

        Returns
        -------
        float
            The computed energy (Eh)
        list[float]
            The gradient (X,Y,Z) for each atom (Eh/Bohr)
        """
        energy = None
        gradient = []
        # read the energy from the output file
        out_path = check_path(output_file)
        with out_path.open() as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    energy = float(line.split()[3])
                    break
        # read the gradient from the .gradient file
        if dograd:
            grad_path = check_path(grad_file)
            natoms_read = 0
            with grad_path.open() as f:
                for line in f:
                    if "$grad" in line:
                        break
                for line in f:
                    fields = line.split()
                    if len(fields) == 4:
                        natoms_read += 1
                    elif len(fields) == 3:
                        gradient += [float(i) for i in fields]
                    elif "$end" in line:
                        break
                if natoms_read != natoms:
                    print(
                        f"Number of atoms read: {natoms_read} does not match the expected: {natoms}"
                    )
                    exit(1)
                if len(gradient) != 3 * natoms:
                    print(
                        f"Number of gradient entries: {len(gradient)} does not match 3x number of atoms: {natoms}"
                    )
                    exit(1)
        if not energy:
            raise ValueError(f"Total energy not found in file {output_file}")
        return energy, gradient

    def run_xtb(
        self,
        calc_data: CalculationData,
        namespace: str,
        method_args: list[str],
        extra_args: list[str],
    ) -> Path:
        """
        Run the xtb program once with a dedicated namespace and redirect its
        STDOUT and STDERR to a per-run logfile.

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        namespace: str
            Value passed to xtb's --namespace; also used to name the logfile
        method_args: list[str]
            Method-specific xtb arguments (e.g., ['--gxtb'] or ['--gfn', '2', '--alpb', 'thf'])
        extra_args: list[str]
            Additional pass-through arguments from the ORCA input

        Returns
        -------
        Path
            Path to the per-run logfile that was written
        """
        # For the xtb package, not all additional command-line arguments (e.g., specifying a solvation model)
        # can be passed before the structure input, so pass the additional args at the end.
        args = [
            str(i)
            for i in [
                calc_data.xyzfile,
                "-c",
                calc_data.charge,
                "-P",
                calc_data.ncores,
                "--namespace",
                namespace,
            ]
        ] + method_args + extra_args
        nue = mult_to_nue(calc_data.mult)
        if nue:
            args += ["-u", str(nue)]
        if calc_data.dograd:
            args += ["--grad"]
        if not calc_data.prog_path:
            raise RuntimeError("Path to program is None.")
        run_log = calc_data.tmp_dir / f"{namespace}.out"
        run_command(calc_data.prog_path, run_log, args)
        return run_log

    def calc(
        self,
        calc_data: CalculationData,
        args_parsed: dict[str, Any],
        args_not_parsed: list[str],
    ) -> tuple[float, list[float]]:
        """
        Routine for calculating energy and optional gradient.
        Runs three xtb jobs (g-xTB gas, GFN2-ALPB, GFN2 gas) and combines them
        into a g-xTB + Δsolv(GFN2-ALPB) energy and gradient.
        Writes ORCA output

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        args_parsed: dict[str, Any]
            Arguments parsed as defined in extend_parser
        args_not_parsed: list[str]
            Arguments not parsed so far. Forwarded to all three xtb runs.

        Returns
        -------
        float
            The combined energy (Eh): E(g-xTB) + [E(GFN2,ALPB) - E(GFN2,gas)]
        list[float]
            Combined flattened gradient vector (Eh/Bohr), if computed, otherwise empty
        """
        # Get parsed options
        prog = args_parsed.get("prog")
        solvent = args_parsed.get("solvent")
        # Set and check the program path if its executable
        calc_data.set_program_path(prog)
        if calc_data.prog_path:
            print(f"Using executable {calc_data.prog_path}")
        else:
            raise FileNotFoundError(
                f"Could not find a valid executable from standard program names: {self.PROGRAM_NAMES}"
            )

        # get the number of atoms from the xyz file (used for parsing each gradient)
        natoms = nat_from_xyzfile(xyz_file=calc_data.xyzfile)

        # Per-component xtb namespaces (also used to name the gradient and logfiles)
        ns_gxtb = f"{calc_data.basename}_gxtb"
        ns_alpb = f"{calc_data.basename}_alpb"
        ns_gas = f"{calc_data.basename}_gas"

        # 1) g-xTB gas-phase
        out_gxtb = self.run_xtb(
            calc_data=calc_data,
            namespace=ns_gxtb,
            method_args=["--gxtb"],
            extra_args=args_not_parsed,
        )
        e_gxtb, g_gxtb = self.read_xtbout(
            output_file=out_gxtb,
            grad_file=f"{ns_gxtb}.gradient",
            natoms=natoms,
            dograd=calc_data.dograd,
        )

        # 2) GFN2-xTB with ALPB(solvent)
        out_alpb = self.run_xtb(
            calc_data=calc_data,
            namespace=ns_alpb,
            method_args=["--gfn", "2", "--alpb", solvent],
            extra_args=args_not_parsed,
        )
        e_alpb, g_alpb = self.read_xtbout(
            output_file=out_alpb,
            grad_file=f"{ns_alpb}.gradient",
            natoms=natoms,
            dograd=calc_data.dograd,
        )

        # 3) GFN2-xTB gas-phase reference
        out_gas = self.run_xtb(
            calc_data=calc_data,
            namespace=ns_gas,
            method_args=["--gfn", "2"],
            extra_args=args_not_parsed,
        )
        e_gas, g_gas = self.read_xtbout(
            output_file=out_gas,
            grad_file=f"{ns_gas}.gradient",
            natoms=natoms,
            dograd=calc_data.dograd,
        )

        # Δsolv combination
        energy = e_gxtb + (e_alpb - e_gas)
        if calc_data.dograd:
            gradient = [g_gxtb[i] + (g_alpb[i] - g_gas[i]) for i in range(3 * natoms)]
        else:
            gradient = []

        # Concatenate the three sub-outputs into the file the framework prints to STDOUT
        with open(calc_data.output_file, "w") as fout:
            for label, fpath in (
                ("g-xTB (gas)", out_gxtb),
                (f"GFN2 ALPB({solvent})", out_alpb),
                ("GFN2 (gas)", out_gas),
            ):
                fout.write(f"\n========== {label}: {fpath.name} ==========\n")
                with open(fpath) as fin:
                    fout.write(fin.read())
            fout.write("\n========== Δsolv combination ==========\n")
            fout.write(f"E(g-xTB, gas)         = {e_gxtb:.10f} Eh\n")
            fout.write(f"E(GFN2 ALPB,{solvent}) = {e_alpb:.10f} Eh\n")
            fout.write(f"E(GFN2, gas)          = {e_gas:.10f} Eh\n")
            fout.write(f"E(total)              = {energy:.10f} Eh\n")

        return energy, gradient


def main() -> None:
    """
    Main routine for execution
    """
    calculator = GxtbAlpbCalc()
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
