# Transport tables

Download the prepared tables from the [runtime release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107).
Extract the required ZIP files into one directory, preserving their folders.
For Legacy, copy the supplied [`legacy_q2.csv`](legacy_q2.csv) there as well.
The C++ connection is shown in the [main README](../README.md#connect-dirac-to-an-application).

| Archive | Contents |
| --- | --- |
| `CCGeant4-runtime-20261107-ionization.zip` | `delta/`, `stopping/`, `global.csv`, `atima_low.csv` |
| `CCGeant4-runtime-20261107-elastic-Z26-40.zip` | Silicon elastic tables for projectile Z = 26–40 and `elastic/tables.index.tsv` |
| `CCGeant4-runtime-20261107-elastic-Z41-55.zip` | Silicon elastic tables for Z = 41–55 |
| `CCGeant4-runtime-20261107-elastic-Z56-70.zip` | Silicon elastic tables for Z = 56–70 |
| `CCGeant4-runtime-20261107-elastic-Z71-82.zip` | Silicon elastic tables for Z = 71–82 |

The ionization archive is sufficient for Dirac with multiple scattering or
native single scattering. Tabulated single scattering uses all four elastic
archives, whose common index refers to the full bank.

Dirac spectra cover projectile Z = 29–92 and energies of 30–2000 MeV in the
table's energy-per-mass coordinate. The GLOBAL and low-energy ATIMA files
contain a wider projectile grid. The elastic tables cover Z = 26–82 on
silicon at 10–5000 MeV/u.

`legacy_q2.csv` contains the charge-scaled stopping inputs for projectile
Z = 29–96 and target Z = 1–84 and 88–96 at 39 energies from 30 to 2000 MeV/u.
It is also available as a separate release download.

The five archives occupy 5.84 GB compressed and 11.24 GB unpacked.
[Download checksums](../README.md#download-checksums),
[GLOBAL and ATIMA sources](../docs/global-atima-data.md), and
[electron and nuclear scattering sources](../docs/dirac-elastic-data.md)
are supplied with the package.
