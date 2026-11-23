"""
GFN2-ALPB implicit-solvation correction for wrappers whose gas-phase method has no ALPB of its own.

E_total = E(method, gas) + [E(GFN2, ALPB) - E(GFN2, gas)]
grad    = grad(method, gas) + [grad(GFN2, ALPB) - grad(GFN2, gas)]

The two GFN2-xTB single points run through the xtb binary (github.com/grimme-lab/xtb).
Without -s/--solvent (or with -s none), no correction is applied and no xtb job is run.

Provides
--------
class: AlpbCorrection(BaseCalc)
    Base for a calculator that adds the correction: its -e/--exe and -s/--solvent options,
    the xtb runs and the combined output
"""

from argparse import ArgumentParser
from pathlib import Path
from typing import Any

from oet.core.base_calc import BaseCalc, CalculationData
from oet.core.misc import check_path, mult_to_nue, run_command

# Solvent names that request the gas phase, i.e. no correction
GAS_PHASE = ("", "none", "gas", "vacuum")


class AlpbCorrection(BaseCalc):
    """
    Subclass it before the gas-phase calculator, e.g. `class UmaAlpbCalc(AlpbCorrection, UmaCalc)`,
    and finish `calc` with `add_alpb_correction`.
    """

    @property
    def PROGRAM_NAMES(self) -> list[str]:
        """Program names to search for in PATH"""
        return ["xtb", "otool_xtb"]

    @classmethod
    def extend_parser(cls, parser: ArgumentParser) -> None:
        """Add the gas-phase calculator's options, then the xtb executable and the solvent.

        Parameters
        ----------
        parser: ArgumentParser
            Parser that should be extended
        """
        super().extend_parser(parser)
        parser.add_argument("-e", "--exe", dest="prog", help="Path to the xtb executable")
        parser.add_argument(
            "-s",
            "--solvent",
            dest="solvent",
            default=None,
            help=(
                "Solvent name passed to xtb's --alpb. If omitted (or set to 'none'/'gas'/'vacuum'), "
                "no solvation correction is applied and no xtb job is run for it."
            ),
        )

    @staticmethod
    def requested_solvent(args_parsed: dict[str, Any]) -> str | None:
        """
        Return the requested ALPB solvent, or None for the gas phase.

        Parameters
        ----------
        args_parsed: dict[str, Any]
            Arguments parsed as defined in extend_parser
        """
        solvent = args_parsed.get("solvent")
        if solvent is None or solvent.strip().lower() in GAS_PHASE:
            return None
        return str(solvent)

    def set_xtb_path(self, calc_data: CalculationData, prog: str | None) -> None:
        """
        Set the xtb executable from -e/--exe, or from PATH if not given.

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        prog: str | None
            Path to the xtb executable given with -e/--exe
        """
        calc_data.set_program_path(prog)
        if calc_data.prog_path:
            print(f"Using executable {calc_data.prog_path}")
        else:
            raise FileNotFoundError(
                f"Could not find a valid executable from standard program names: {self.PROGRAM_NAMES}"
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
    ) -> tuple[float, list[float], Path]:
        """
        Run the xtb program once with a dedicated namespace, redirecting its
        STDOUT and STDERR to a per-run logfile, and read its energy and gradient.

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
        float
            The computed energy (Eh)
        list[float]
            The gradient (X,Y,Z) for each atom (Eh/Bohr), if computed, otherwise empty
        Path
            Path to the per-run logfile that was written
        """
        # For the xtb package, not all additional command-line arguments (e.g., specifying a solvation model)
        # can be passed before the structure input, so pass the additional args at the end.
        args = (
            [
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
            ]
            + method_args
            + extra_args
        )
        nue = mult_to_nue(calc_data.mult)
        if nue:
            args += ["-u", str(nue)]
        if calc_data.dograd:
            args += ["--grad"]
        if not calc_data.prog_path:
            raise RuntimeError("Path to program is None.")
        run_log = calc_data.tmp_dir / f"{namespace}.out"
        run_command(calc_data.prog_path, run_log, args)
        energy, gradient = self.read_xtbout(
            output_file=run_log,
            grad_file=f"{namespace}.gradient",
            natoms=calc_data.natoms,
            dograd=calc_data.dograd,
        )
        return energy, gradient, run_log

    def add_alpb_correction(
        self,
        calc_data: CalculationData,
        method: str,
        energy: float,
        gradient: list[float],
        solvent: str | None,
        sub_outputs: list[tuple[str, Path]],
        extra_args: list[str],
    ) -> tuple[float, list[float]]:
        """
        Add the GFN2-ALPB correction to a gas-phase energy and gradient, if a solvent is requested,
        and write all xtb logs and the energy summary to the output file printed to STDOUT.
        The xtb executable must already be set with `set_xtb_path`.

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        method: str
            Name of the gas-phase method, used in the summary
        energy: float
            Gas-phase energy of the method (Eh)
        gradient: list[float]
            Gas-phase gradient of the method (Eh/Bohr), empty if not computed
        solvent: str | None
            ALPB solvent, or None for the gas phase
        sub_outputs: list[tuple[str, Path]]
            Labelled logfiles of the gas-phase method to include in the output
        extra_args: list[str]
            Additional pass-through arguments from the ORCA input, forwarded to both xtb runs

        Returns
        -------
        float
            The combined energy (Eh): E(method) + [E(GFN2,ALPB) - E(GFN2,gas)],
            or just E(method) if no solvent was requested
        list[float]
            Combined flattened gradient vector (Eh/Bohr), if computed, otherwise empty
        """
        if solvent is None:
            summary = [
                "\n========== gas phase (no solvation correction) ==========\n",
                f"E({method}, gas) = {energy:.10f} Eh\n",
            ]
            total_energy, total_gradient = energy, gradient
        else:
            e_alpb, g_alpb, out_alpb = self.run_xtb(
                calc_data=calc_data,
                namespace=f"{calc_data.basename}_alpb",
                method_args=["--gfn", "2", "--alpb", solvent],
                extra_args=extra_args,
            )
            e_gas, g_gas, out_gas = self.run_xtb(
                calc_data=calc_data,
                namespace=f"{calc_data.basename}_gas",
                method_args=["--gfn", "2"],
                extra_args=extra_args,
            )
            total_energy = energy + (e_alpb - e_gas)
            if calc_data.dograd:
                total_gradient = [
                    g + (a - b) for g, a, b in zip(gradient, g_alpb, g_gas, strict=True)
                ]
            else:
                total_gradient = []
            sub_outputs = sub_outputs + [
                (f"GFN2 ALPB({solvent})", out_alpb),
                ("GFN2 (gas)", out_gas),
            ]
            summary = [
                "\n========== Δsolv combination ==========\n",
                f"{f'E({method}, gas)':<22} = {energy:.10f} Eh\n",
                f"{f'E(GFN2 ALPB,{solvent})':<22} = {e_alpb:.10f} Eh\n",
                f"{'E(GFN2, gas)':<22} = {e_gas:.10f} Eh\n",
                f"{'E(total)':<22} = {total_energy:.10f} Eh\n",
            ]

        with open(calc_data.output_file, "w") as fout:
            for label, fpath in sub_outputs:
                fout.write(f"\n========== {label}: {fpath.name} ==========\n")
                with open(fpath) as fin:
                    fout.write(fin.read())
            fout.writelines(summary)

        return total_energy, total_gradient
