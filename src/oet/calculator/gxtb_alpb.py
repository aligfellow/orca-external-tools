#!/usr/bin/env python3
"""
Wrapper that adds an ALPB(solvent) implicit-solvation correction (computed at the GFN2-xTB level)
on top of a g-xTB gas-phase energy and gradient, compatible with ORCA's ExtTool interface.

E_total = E(g-xTB, gas) + [E(GFN2, ALPB) - E(GFN2, gas)]
grad    = grad(g-xTB, gas) + [grad(GFN2, ALPB) - grad(GFN2, gas)]

Uses the modified xtb binary (github.com/grimme-lab/g-xTB), which supports both --gxtb
and standard GFN2-xTB. Three xtb single-points are performed per ORCA gradient call.

Without -s/--solvent (or with -s none), no correction is applied and only the single
gas-phase g-xTB run is performed. The correction is shared with other wrappers, see oet.core.alpb.

Provides
--------
class: GxtbAlpbCalc(AlpbCorrection)
    Class for performing a g-xTB calculation with a GFN2-ALPB Δsolv correction
main: function
    Main function
"""

from typing import Any

from oet.core.alpb import AlpbCorrection
from oet.core.base_calc import CalculationData


class GxtbAlpbCalc(AlpbCorrection):
    def calc(
        self,
        calc_data: CalculationData,
        args_parsed: dict[str, Any],
        args_not_parsed: list[str],
    ) -> tuple[float, list[float]]:
        """
        Routine for calculating energy and optional gradient.
        Runs g-xTB in the gas phase and adds the GFN2-ALPB correction if a solvent is requested.
        Writes ORCA output

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        args_parsed: dict[str, Any]
            Arguments parsed as defined in extend_parser
        args_not_parsed: list[str]
            Arguments not parsed so far. Forwarded to all xtb runs.

        Returns
        -------
        float
            The combined energy (Eh): E(g-xTB) + [E(GFN2,ALPB) - E(GFN2,gas)],
            or just E(g-xTB) if no solvent was requested
        list[float]
            Combined flattened gradient vector (Eh/Bohr), if computed, otherwise empty
        """
        solvent = self.requested_solvent(args_parsed)
        self.set_xtb_path(calc_data, args_parsed.get("prog"))

        # g-xTB gas-phase
        energy, gradient, out_gxtb = self.run_xtb(
            calc_data=calc_data,
            namespace=f"{calc_data.basename}_gxtb",
            method_args=["--gxtb"],
            extra_args=args_not_parsed,
        )

        return self.add_alpb_correction(
            calc_data=calc_data,
            method="g-xTB",
            energy=energy,
            gradient=gradient,
            solvent=solvent,
            sub_outputs=[("g-xTB (gas)", out_gxtb)],
            extra_args=args_not_parsed,
        )


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
