# Full transport tables

Download the five ZIP archives from the [CCGeant4 runtime release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107). The [download checksums](../README.md#download-checksums) identify these files. The archives contain **188 runtime files**, about 5.84 GB compressed and 11.24 GB unpacked. The installer also copies the Legacy stopping table supplied in this checkout. A small [transport example](../examples/transport/) can run without the full bank.

| Archive | Compressed bytes | Contents |
| --- | ---: | --- |
| `CCGeant4-runtime-20261107-ionization.zip` | 1,230,221,341 | 64 Dirac spectra, 64 stopping tables, GLOBAL transition cross sections and low-energy ATIMA stopping |
| `CCGeant4-runtime-20261107-elastic-Z26-40.zip` | 1,201,974,306 | Silicon elastic tables for projectile Z = 26–40 and the index for the full bank |
| `CCGeant4-runtime-20261107-elastic-Z41-55.zip` | 1,211,111,281 | Silicon elastic tables for Z = 41–55 |
| `CCGeant4-runtime-20261107-elastic-Z56-70.zip` | 1,216,191,021 | Silicon elastic tables for Z = 56–70 |
| `CCGeant4-runtime-20261107-elastic-Z71-82.zip` | 975,596,837 | Silicon elastic tables for Z = 71–82 |

The ionization archive covers the full Dirac domain, Z = 29–92. The GLOBAL and low-energy ATIMA CSVs retain their original, wider projectile grid. The 57 elastic binaries cover 10–5000 MeV/u on silicon. For tabulated single scattering use all four elastic archives together because the runtime index identifies every pair. For Dirac with multiple scattering or native single scattering the ionization archive is sufficient.

For charge-scaled ATIMA stopping, Legacy reads [`legacy_q2.csv`](legacy_q2.csv), included in the checkout and available as a [separate release download](https://github.com/andreiazhakin/CCGeant4/releases/download/runtime-20261107/legacy_q2.csv). The [Legacy configuration](../README.md#full-tables-and-application-setup) gives the three input paths and the screening setting.

The table contains 246,636 rows for projectile Z = 29–96 and target Z = 1–84 and 88–96 at 39 energies from 30 to 2000 MeV/u. The file is 19,353,020 bytes.

## Installation

Put the archives in one directory. From the repository root run either command below with Python 3. The destination `DATA` must be a new directory. Allow 11.24 GB of free space in addition to the archives for the full installation.

```sh
# Install the complete bank
python data/install_runtime.py --archives /path/to/archives --output /path/to/DATA

# Or install only the electron and charge-exchange data
python data/install_runtime.py --archives /path/to/archives --ionization-only --output /path/to/DATA
```

The complete installation contains 189 files, including `legacy_q2.csv`.
The ionization-only installation contains 131 files and supports both Dirac
and Legacy. `--verify-only` reads the ZIP members without extracting them.

The resulting layout matches [the Geant4 configuration example](../docs/physics-configuration.md):

```text
DATA/
  delta/Z029.csv ... Z092.csv
  stopping/Z029.csv ... Z092.csv
  global.csv
  atima_low.csv
  legacy_q2.csv
  elastic/tables.index.tsv
  elastic/*_Si28.eltb
```

These archives contain the numerical transport inputs. Sources, upstream terms and the treatment of calculated results are described for [GLOBAL and ATIMA](../docs/global-atima-data.md) and [Dirac spectra and nuclear scattering](../docs/dirac-elastic-data.md). The archives exclude upstream solvers, spline databases, potential coefficients, journal PDFs and validation events.
