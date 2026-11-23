# GFN2-ALPB Δsolv corrections — options and usage

This document covers the wrappers that add an implicit-solvation correction,
evaluated at the **GFN2-xTB / ALPB** level, to a gas-phase method that has no
ALPB of its own:

| Wrapper | Gas-phase method | Server |
|---|---|---|
| `oet_gxtb_alpb` | g-xTB, through the g-xTB-enabled `xtb` binary | — |
| `oet_uma_alpb` | UMA (fairchem), with the options of `oet_uma` | `oet_server uma_alpb` |

For installation, generic ORCA wiring, and the `extinp` / `engrad` interface
see the main [README](../README.md).

## What it computes

Per ORCA gradient call the wrapper computes the gas-phase method and two xtb
single points, and combines them so that the electronic structure comes from
the method while the solvation contribution comes from GFN2-ALPB:

```
E_total = E(method, gas) + [ E(GFN2, ALPB) − E(GFN2, gas) ]
grad    = grad(method, gas) + [ grad(GFN2, ALPB) − grad(GFN2, gas) ]
```

The two xtb runs are `xtb --gfn 2 --alpb <solvent>` and `xtb --gfn 2`.
Without `-s/--solvent` (or with `-s none`), no correction is applied and no xtb
job is run for it: `oet_gxtb_alpb` is then plain g-xTB and `oet_uma_alpb` plain UMA.

The correction lives in `oet.core.alpb.AlpbCorrection`. Another method gains it by
subclassing it before its own calculator and finishing `calc` with
`add_alpb_correction`, as `UmaAlpbCalc(AlpbCorrection, UmaCalc)` does.

## Requirements

1. **An `xtb` binary.**
   - Any xtb (github.com/grimme-lab/xtb) serves the GFN2 terms.
     `oet_gxtb_alpb` also drives `xtb --gxtb`, which requires
     [github.com/grimme-lab/g-xTB](https://github.com/grimme-lab/g-xTB).
   - It must be on `PATH` (the wrappers search for `xtb`, then `otool_xtb`),
     or pass the full path with `-e/--exe`.
2. **For `oet_uma_alpb`**, the UMA installation of `oet_uma` (`python install.py -e uma`)
   and access to the model files, see the main README.

## Arguments

Both wrappers:

| Argument | Values / Description | Default |
|---|---|---|
| `-s, --solvent` | Solvent name passed to xtb's `--alpb` for the Δsolv correction; `none` is the gas phase | none (gas phase) |
| `-e, --exe` | Full path to the xtb executable | searched on `PATH` |

`oet_uma_alpb` also takes the options of `oet_uma`:

| Argument | Values / Description | Default |
|---|---|---|
| `-m, --model` | UMA base model, e.g. `uma-s-1p1`, `uma-s-1p2` | `uma-s-1p1` |
| `-t, --task` | UMA task, e.g. `omol` | `omol` |
| `-d, --device` | `cpu` or `cuda` | `cpu` |
| `-c, --cachedir` | Cache directory of the model files (or `FAIRCHEM_CACHE_DIR`) | `oet` assets |
| `-o, --offline` | Force offline mode | off |

Any additional, unrecognised arguments are forwarded verbatim to **all** xtb
runs (*e.g.* `--acc 0.5`). Charge and multiplicity come from the ORCA
`extinp` file and are translated to xtb's `-c` / `-u` automatically.

## Examples

ORCA optimisation with g-xTB in acetonitrile, explicit xtb path:

```text
! ExtOpt
%method
  ProgExt "/path/to/oet_gxtb_alpb"
  Ext_Params "-s acetonitrile -e /path/to/xtb"
end
* xyzfile 0 1 mol.xyz
```

ORCA transition-state search with a numerical Hessian:

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

UMA in methanol, standalone (the model is loaded on every call):

```text
! ExtOpt Opt
%method
  ProgExt "/path/to/oet_uma_alpb"
  Ext_Params "-m uma-s-1p2 -s methanol -e /path/to/xtb"
end
* xyzfile 0 1 mol.xyz
```

UMA in methanol through a server, which keeps the model loaded between calls:

```bash
oet_server uma_alpb -n 4 -b 127.0.0.1:9000
```

```text
%method
  ProgExt "/path/to/oet_client"
  Ext_Params "-b 127.0.0.1:9000 -m uma-s-1p2 -s methanol -e /path/to/xtb"
end
```

To check - can we run opt + freq in a single job ? Previously this did *not* run the frequency too.
