# Table converter

`convert.py` converts calculated Dirac spectra or nuclear-scattering
distributions into tables read by CCGeant4. It needs Python 3.10+ and NumPy.
Ready-to-use tables are also available in the
[runtime release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107).

Run the commands below from the repository root. Input paths may point
anywhere on disk. Use a new or empty output directory.

## Dirac

Three inputs are needed:

- `--spectra`: the Dirac ZIP containing `grid.npz`, `bare_Z*.npz` and
  `Z*_N*.npz` under `full_sweep/results/`.
- `--global`: the GLOBAL transition CSV, including projectile and target Z,
  incident energy, bound-electron counts and transition cross sections.
- `--atima`: the ATIMA 1.41 range-spline archive `splines_gz_1.41.tar`.

```sh
python converters/convert.py dirac --spectra inputs/dirac.zip --global inputs/global.csv --atima inputs/splines_gz_1.41.tar --output runtime --projectiles 29 54 --targets 1 14
```

`--projectiles` and `--targets` select atomic numbers. Omit either flag to
use all available values for that selection. Omit both to convert all pairs.
The inputs must cover every selected projectile–target pair. Include every
element used in the application's materials, including hydrogen for
`G4_Galactic`.

The result can be connected directly through `G4CCPhysicsOptions`:

```text
runtime/
  delta/Z029.csv       # Electron spectra, one file per projectile
  delta/Z054.csv
  stopping/Z029.csv    # Stopping coefficients for the selected targets
  stopping/Z054.csv
  global.csv           # Selected GLOBAL rows
  atima_low.csv        # Complete supplied low-energy ATIMA table
```

The converter copies the complete low-energy table from
[`data/atima_low.csv.xz`](../data/atima_low.csv.xz). It does not recalculate or
filter it. All original spectrum energies and calculated charge states are
retained. The stopping grids are refined to reproduce the ATIMA equilibrium
stopping.

## Elastic scattering

`--source` points to a directory of calculated `.npz` distributions with
matching `.json` metadata files. The files may also be in its `packed/`
subdirectory.

```sh
python converters/convert.py elastic --source inputs/elastic --output runtime/elastic --projectiles 26 29 --targets 14
```

The same selection flags apply. The output contains one `.eltb` file per
selected pair and `tables.index.tsv`. Keep them together and pass the index
path as `cc.scattering.tables` in `TabulatedSingleScattering` mode.
The converter preserves the calculated arrays. It does not calculate a
nuclear-scattering model.

The [format description](FORMAT.md) gives the input fields, units and output
layout. See [model configuration](../README.md#connect-dirac-to-an-application),
[electron and nuclear scattering data](../docs/dirac-elastic-data.md), and
[GLOBAL and ATIMA data](../docs/global-atima-data.md) for their use and sources.
