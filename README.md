# CCGeant4

CCGeant4 adds charge-state-resolved heavy-ion transport to Geant4 11.4.0.
Dirac uses tabulated electron spectra and stopping powers for each ionic
charge. Legacy uses charge-scaled ATIMA stopping and a Mott-corrected
electron spectrum. Both modes simulate electron capture and loss using
GLOBAL cross sections.

The package contains the Geant4 patch, its added and modified source files,
table converters, and a [transport example](examples/transport/README.md)
with reference results and a script for analyzing and comparing calculations.

## Apply the patch and build

Apply the patch to a clean Geant4 11.4.0 source tree:

```sh
git -C /path/to/geant4-v11.4.0 apply /path/to/CCGeant4/patches/ccgeant4-geant4-11.4.0.patch
```

`-C` selects the Geant4 source directory. The argument after `apply` is the
patch file. Then build and install Geant4 normally:

```sh
cmake -S /path/to/geant4-v11.4.0 -B build-geant4 -DCMAKE_BUILD_TYPE=Release -DGEANT4_INSTALL_DATA=ON -DCMAKE_INSTALL_PREFIX=/path/to/geant4-install
cmake --build build-geant4 --config Release --parallel
cmake --install build-geant4 --config Release
```

Build applications against this patched installation.

## Tables

The full tables are supplied in the [table release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107).
Extract `CCGeant4-runtime-20261107-ionization.zip` into a directory such as
`DATA`, preserving its folders. For Legacy, also copy the supplied
[`data/legacy_q2.csv`](data/legacy_q2.csv) into that directory.

```text
DATA/
  delta/          # Dirac electron spectra
  stopping/       # Dirac stopping tables
  global.csv      # Capture and loss cross sections
  atima_low.csv   # Low-energy equilibrium stopping
  legacy_q2.csv   # Legacy charge-scaled stopping
```

The example in the full release ZIP includes a smaller table set for Cu, Xe and U in silicon and
hydrogen, plus Cu–Si and Fe–Si elastic tables for its single-scattering
cases. These files are already unpacked in `examples/transport/data/`.
A Git clone contains the example sources. Copy `examples/transport/data/`
from the [full release ZIP](https://github.com/andreiazhakin/CCGeant4/releases),
or supply prepared tables with `--tables`.
For the full silicon elastic bank, extract the four elastic archives into
the same `DATA` directory. Archive contents and coverage are listed in
[data/README.md](data/README.md).

## Connect Dirac to an application

Pass the table paths through `G4CCPhysicsOptions`. For an existing modular
physics list, replace its electromagnetic constructor before initialization:

```cpp
#include "FTFP_BERT.hh"
#include "G4CCPhysicsOptions.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

void SetPhysics(G4RunManager& runManager)
{
    G4CCPhysicsOptions cc;
    cc.mode = G4CCMode::Dirac;
    cc.tables.delta    = "DATA/delta";
    cc.tables.stopping = "DATA/stopping";
    cc.tables.global   = "DATA/global.csv";
    cc.tables.atimaLow = "DATA/atima_low.csv";

    auto* physics = new FTFP_BERT;
    physics->ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge(cc));
    physics->SetDefaultCutValue(25*um);
    runManager.SetUserInitialization(physics);
}
```

`cc.mode` selects the model and defaults to Dirac. `delta` and `stopping`
point to directories of CSV tables, or to a single projectile CSV.
`global` and `atimaLow` point to files. Paths are relative to the application's
working directory, or can be absolute. Choose the production cut for your
detector and materials.

For an empty `G4VModularPhysicsList`, use
`RegisterPhysics(new G4EmStandardPhysics_option4_CurrentCharge(cc))` instead
of `ReplacePhysics`. The example does this directly in
[`main.cc`](examples/transport/main.cc), with separate branches for Dirac,
Legacy and standard emopt4.

Multiple scattering is the default. To select native ion single scattering,
set this before creating the physics constructor:

```cpp
cc.scattering.mode = G4CCScatteringMode::NativeSingleScattering;
```

For tabulated single scattering, use instead:

```cpp
cc.scattering.mode = G4CCScatteringMode::TabulatedSingleScattering;
cc.scattering.tables = "DATA/elastic/tables.index.tsv";
```

## Connect Legacy

Legacy uses the same constructor, with its charge-scaled stopping table:

```cpp
G4CCPhysicsOptions cc;
cc.mode = G4CCMode::Legacy;
cc.tables.global = "DATA/global.csv";
cc.tables.legacyStopping = "DATA/legacy_q2.csv";
cc.tables.atimaLow = "DATA/atima_low.csv";

physics->ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge(cc));
```

All three paths point to CSV files. The scattering settings above apply
to both modes. Set the mode and paths before initialization.
The lower Legacy transition remains fixed at 30–40 MeV/u.

## Build and run the example

```sh
cd examples/transport
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/geant4-install
cmake --build build --config Release --parallel
./build/cc_transport --macro macros/cu_400.mac --output results/dirac_cu.csv
python analyze.py results/dirac_cu.csv --reference reference/summary.csv --output results/summary.csv
```

Initialize the environment for the patched Geant4 installation before running.
With Visual Studio, the executable is `build/Release/cc_transport.exe`.
The default mode is Dirac. Use `--model legacy` or `--model emopt4` to change
it. The C++ application reads the supplied `data/` directly. Only analysis
requires Python 3.10+. Add `--plot` to save the distributions. Only that
option requires Matplotlib.
Beam settings, macros, scattering choices and the output columns are
described in the [example README](examples/transport/README.md).

## Converters and article data

The [table converter](converters/README.md) prepares Dirac or elastic
transport tables from calculated input data. It requires NumPy and can
select individual projectile and target elements.

The [article-data release](https://github.com/andreiazhakin/CCGeant4/releases/tag/article-data-20261109)
contains simulation events, run settings, experimental and reference data,
and the article's figures.

## Download checksums

| Table archive | SHA-256 |
| --- | --- |
| `CCGeant4-runtime-20261107-ionization.zip` | `0ed618a78ffe98f06fd6125a2692b1ba0b891432534753d86bda04fb0e56ba1d` |
| `CCGeant4-runtime-20261107-elastic-Z26-40.zip` | `38b5f2f859ba06e440ae4a80f6c622823c0d1f52a743d646dfdb4ddbde6b4b9b` |
| `CCGeant4-runtime-20261107-elastic-Z41-55.zip` | `e803dff4502bfa9ccfb5dced17bc881256f83d99cde08a0770ef270cf98aa358` |
| `CCGeant4-runtime-20261107-elastic-Z56-70.zip` | `eb5668bb0a25b4d800067773440145f82b769839d13c88689a39868fa6c9b6a3` |
| `CCGeant4-runtime-20261107-elastic-Z71-82.zip` | `3f919ab7d576508c76e25331dfe41dcc96ba57f8fd76916a319087affb3fa0c3` |

## License

See [LICENSING.md](LICENSING.md) for the project license, Geant4 terms and
numerical input credits.

This product includes software developed by Members of the Geant4 Collaboration (http://cern.ch/geant4).
