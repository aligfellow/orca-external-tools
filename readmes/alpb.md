# g-xTB + GFN2-ALPB Δsolv calculator — options and usage

This document covers the `oet_gxtb_alpb` wrapper, which runs a **g-xTB**
gas-phase calculation and adds an implicit-solvation correction evaluated at
the **GFN2-xTB / ALPB** level. For installation, generic ORCA wiring, and the
`extinp` / `engrad` interface see the main [README](../README.md).

## What it computes

Per ORCA gradient call the wrapper performs **three** xtb single points and
combines them so that the electronic structure comes from g-xTB while
the solvation contribution comes from GFN2-ALPB:

```
E_total = E(g-xTB, gas) + [ E(GFN2, ALPB) − E(GFN2, gas) ]
grad    = grad(g-xTB, gas) + [ grad(GFN2, ALPB) − grad(GFN2, gas) ]
```

The three runs are: `xtb --gxtb`, `xtb --gfn2 --alpb <solvent>`, and
`xtb --gfn2`. Because g-xTB has no native ALPB, the GFN2 difference is used as
a Δsolv estimate. Note the cost: every energy/gradient evaluation ORCA requests
triggers three xtb jobs.

## Requirements

1. **A g-xTB `xtb` binary.** 
  - The wrapper drives `xtb --gxtb`, which requires [github.com/grimme-lab/g-xTB](https://github.com/grimme-lab/g-xTB) that accepts `--gxtb` in addition to standard `--gfn2`.
    - It must be on `PATH` (the wrapper searches for `xtb`, then `otool_xtb`),
     or pass the full path with `-e/--exe`.

## Arguments

| Argument | Values / Description | Default |
|---|---|---|
| `-s, --solvent` | Solvent name passed to xtb's `--alpb` for the Δsolv correction | `thf` |
| `-e, --exe` | Full path to the (g-xTB-enabled) xtb executable | searched on `PATH` |

Any additional, unrecognised arguments are forwarded verbatim to **all three**
xtb runs (*e.g.* `--acc 0.5`). Charge and multiplicity come from the ORCA
`extinp` file and are translated to xtb's `-c` / `-u` automatically.

## Examples

ORCA optimisation in acetonitrile, explicit xtb path:

```text
! ExtOpt
%method
  ProgExt "/path/to/oet_gxtb_alpb"
  Ext_Params "-s acetonitrile -e /path/to/xtb"
end
* xyzfile 0 1 mol.xyz
```

ORCA transition-state search with a numerical Hessian (the reason for this wrapper):

```text
! ExtOpt OptTS xyzfile
%maxcore 600
%pal
  nprocs 24
end
%method
  ProgExt "/path/to/oet_gxtb_alpb"
  Ext_Params "-s acetonitrile"
end
%geom
  Calc_Hess true
  NumHess true
  TS_Mode {B 1 2} end
  TS_Active_Atoms {1 2} end
end
* xyzfile 0 1 ts.xyz
```

To check - can we run opt + freq in a single job ? Previously this did *not* run the frequency too.