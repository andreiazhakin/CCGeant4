# Ionization table converter

Transport users can run the bundled example directly. This directory is for
rebuilding the ionization tables from the project's frozen-shell spectra, the
GLOBAL cross-section CSV and the ATIMA 1.41 range splines. The converter needs
Python 3.11+ and NumPy. Set aside a fresh output directory for the binary packs
and another for the runtime CSVs.

From the repository root, with the three inputs already on disk:

```sh
python -m pip install -r converters/dirac/requirements.txt
python converters/dirac/convert.py --dirac Dirac_ALL_N_Z29-92_2026-11-01.zip --global-csv global.csv --atima-splines splines_gz_1.41.tar --output packs
python converters/dirac/validate.py --tables packs --raw-dirac Dirac_ALL_N_Z29-92_2026-11-01.zip --raw-global global.csv --raw-atima splines_gz_1.41.tar
python converters/dirac/export_dirac_csv.py --source packs --output runtime
python converters/dirac/validate_csv_export.py --source packs --csv-root runtime --report csv-check.json
```

For a single projectile, add `--z-min 54 --z-max 54` to conversion and
`--z 54` to the CSV exporter. The full production domain is Z = 29 through 92.
The export checker compares every CSV field with its source array. The converter records the inputs and numerical checks in its reports.

The export produces `delta/Znnn.csv` (electron spectra) and
`stopping/Znnn.csv` (charge-dependent soft coefficient and equilibrium
stopping, including left limits). Transport also needs the GLOBAL transition
CSV and the six-column low-energy ATIMA reference CSV. To add the reference
mass to an existing five-column low-energy project export:

```sh
python converters/dirac/convert_low_atima_reference.py --source low-five-columns.csv --atima-splines splines_gz_1.41.tar --output runtime/atima_low.csv
```

The [format specification](FORMAT_RU.md) describes the grids, units, source
columns and normalization. The [atomic calculation](../../docs/dirac-elastic-data.md)
and [GLOBAL and ATIMA records](../../docs/global-atima-data.md) give the input
identifiers and terms. A fresh calculation with stock upstream programs may
differ from the historical inputs. The [runtime manifest](../../data/README.md)
identifies the exact final tables used by this version.
