#!/usr/bin/env python3
"""
Wrapper that adds an ALPB(solvent) implicit-solvation correction (computed at the GFN2-xTB level)
on top of a UMA gas-phase energy and gradient, compatible with ORCA's ExtTool interface.

E_total = E(UMA, gas) + [E(GFN2, ALPB) - E(GFN2, gas)]
grad    = grad(UMA, gas) + [grad(GFN2, ALPB) - grad(GFN2, gas)]

UMA takes the same options as `oet_uma`; the two GFN2-xTB single points run through any
xtb binary (github.com/grimme-lab/xtb) on every call. In server mode (`oet_server uma_alpb`)
the UMA model stays loaded between calls.

Without -s/--solvent (or with -s none), no correction is applied and this is plain UMA.
The correction is shared with other wrappers, see oet.core.alpb.

Provides
--------
class: UmaAlpbCalc(AlpbCorrection, UmaCalc)
    Class for performing a UMA calculation with a GFN2-ALPB Δsolv correction
main: function
    Main function
"""

from typing import Any

from oet.calculator.uma import UmaCalc
from oet.core.alpb import AlpbCorrection
from oet.core.base_calc import CalculationData


class UmaAlpbCalc(AlpbCorrection, UmaCalc):
    def calc(
        self,
        calc_data: CalculationData,
        args_parsed: dict[str, Any],
        args_not_parsed: list[str],
    ) -> tuple[float, list[float]]:
        """
        Routine for calculating energy and optional gradient.
        Runs UMA in the gas phase and adds the GFN2-ALPB correction if a solvent is requested.
        Writes ORCA output

        Parameters
        ----------
        calc_data: CalculationData
            Object with calculation data for the run
        args_parsed: dict[str, Any]
            Arguments parsed as defined in extend_parser
        args_not_parsed: list[str]
            Arguments not parsed so far. Forwarded to both xtb runs.

        Returns
        -------
        float
            The combined energy (Eh): E(UMA) + [E(GFN2,ALPB) - E(GFN2,gas)],
            or just E(UMA) if no solvent was requested
        list[float]
            Combined flattened gradient vector (Eh/Bohr), if computed, otherwise empty
        """
        solvent = self.requested_solvent(args_parsed)
        # Find xtb before loading the model, so that a missing binary fails fast
        if solvent is not None:
            self.set_xtb_path(calc_data, args_parsed.get("prog"))

        # UMA gas-phase
        energy, gradient = super().calc(calc_data, args_parsed, args_not_parsed)

        return self.add_alpb_correction(
            calc_data=calc_data,
            method="UMA",
            energy=energy,
            gradient=gradient,
            solvent=solvent,
            sub_outputs=[],
            extra_args=args_not_parsed,
        )


def main() -> None:
    """
    Main routine for execution
    """
    calculator = UmaAlpbCalc()
    inputfile, args, args_not_parsed = calculator.parse_args()
    calculator.run(inputfile=inputfile, args_parsed=args, args_not_parsed=args_not_parsed)


# Python entry point
if __name__ == "__main__":
    main()
