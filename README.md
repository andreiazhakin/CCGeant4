# CCGeant4

CCGeant4 adds charge-state-resolved heavy-ion transport to Geant4 11.4.0.
It combines GLOBAL capture and loss cross sections, ATIMA stopping powers,
and explicit delta-electron production. The Dirac mode uses tabulated
electron spectra with projectile-electron screening. Legacy uses
charge-scaled stopping and a Mott-corrected electron spectrum.
Multiple scattering, native ion single scattering and tabulated nuclear
elastic scattering can be selected independently.

## Contents

| Directory | Contents |
| --- | --- |
| `patches/` | Patch for Geant4 11.4.0 and its file list. |
| `source_changes/` | The same 64 added or modified source files. |
| `examples/transport/` | One transport application, runner, analysis script, small input table set and reference results. |
| `converters/` | Ionization and elastic table converters and optional numerical validators. |
| `data/` | Full-bank installer, archive index and the Legacy stopping table. |
| `docs/` | Physics configuration, data sources and attribution. |

The [article data](https://github.com/andreiazhakin/CCGeant4/releases/tag/article-data-20261109)
are distributed separately from this software package.

## Build Geant4

Use a clean [Geant4 11.4.0 source tree](https://gitlab.cern.ch/geant4/geant4/-/tree/v11.4.0),
a C++17 compiler and CMake. Apply the patch before configuring Geant4:

```sh
git -C /path/to/geant4-v11.4.0 apply /absolute/path/to/CCGeant4/patches/ccgeant4-geant4-11.4.0.patch
cmake -S /path/to/geant4-v11.4.0 -B build-g4 -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DBUILD_SHARED_LIBS=ON -DGEANT4_BUILD_MULTITHREADED=ON -DGEANT4_INSTALL_DATA=ON -DCMAKE_INSTALL_PREFIX=/path/to/patched-install
cmake --build build-g4 --config Release --parallel
cmake --install build-g4 --config Release
```

Initialize the environment for that Geant4 installation before running an
application. Rebuild Geant4 and linked applications when updating the patch.
An optional source check is available as
`python tools/verify_patch.py --base /path/to/geant4-v11.4.0.zip`.

## Build and run the example

The example uses ordinary Geant4 detector, physics-list and action classes.
It records energy loss, energy deposited in the slab, outgoing charge and
angle. The [example README](examples/transport/README.md) explains the files
and how to change the beam, material, thickness and scattering model.

```sh
cmake -S examples/transport -B build-example -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/patched-install
cmake --build build-example --config Release --parallel
python examples/transport/run_example.py --executable build-example/cc_transport --output results --events 10000
python examples/transport/analyze.py results/dirac_cu.csv results/dirac_xe.csv results/dirac_u.csv --reference examples/transport/reference/summary.csv --output results/summary.csv
```

With Visual Studio, use `build-example/Release/cc_transport.exe`.
The runner uses the supplied Cu, Xe and U tables and does not require the
full bank. Select `--model legacy` or `--model emopt4` for the other modes.
Each configuration runs in a separate process.

Running and numerical comparison require only Python 3.10+ with its standard
library. Add `--plot` to the analysis command for distributions with the
reference histograms; only plotting requires Matplotlib.

For single scattering, use `--case ss_cu ss_native fe_ss`. The tabulated
cases require the separate elastic bank:

```sh
python examples/transport/run_example.py --executable build-example/cc_transport --case ss_cu ss_native fe_ss --elastic-tables /path/to/DATA/elastic/tables.index.tsv --output results-ss --events 10000
python examples/transport/analyze.py results-ss/dirac_ss_cu.csv results-ss/dirac_ss_native.csv results-ss/dirac_fe_ss.csv --reference examples/transport/reference/summary.csv --output results-ss/summary.csv --plot
```

Reference results contain 100,000 events per case. Comparisons use statistical
uncertainties rather than requiring identical random sequences. The Cu
single-scattering energy widths use an iterative ±3σ core; angles always
use the full exiting sample. The analysis reads this selection from the
reference summary.

## Full tables and application setup

The [runtime release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107)
contains five table archives, about 5.84 GB compressed. See
[data/README.md](data/README.md) for installation and coverage.
The large tables are release downloads, not part of the source ZIP.

The [configuration guide](docs/physics-configuration.md) shows how to pass
Dirac table paths through `G4CCPhysicsOptions` in another application.
Legacy uses these three environment variables before initialization:

| Variable | Input |
| --- | --- |
| `G4CC_GLOBAL_TABLE` | `DATA/global.csv` |
| `G4CC_ATIMA_Q_DEDX_TABLE` | `DATA/legacy_q2.csv` |
| `G4CC_ATIMA_EQUILIBRIUM_DEDX_TABLE` | `DATA/atima_low.csv` |

The example sets these paths itself. Legacy's binding-energy override defaults
to 1000 GeV. To set it explicitly, use `G4CC_BINDING_OVERRIDE_GEV=1000` with
`G4CC_BINDING_OVERRIDE_MEV` unset. The supplied low-energy table includes the
reference target masses needed for isotope-dependent stopping.

## Converters

- [Dirac/GLOBAL/ATIMA](converters/dirac/README.md): convert calculated spectra
  and stopping inputs into runtime tables; optional validators compare the
  resulting moments and exported values.
- [Elastic scattering](converters/elastic/README.md): convert calculated
  angular distributions or generate a simplified Coulomb table for an
  isotope pair.

The converters require NumPy. They are not needed to use the supplied tables.

## Download checksums

| Runtime archive | SHA-256 |
| --- | --- |
| `CCGeant4-runtime-20261107-ionization.zip` | `0ed618a78ffe98f06fd6125a2692b1ba0b891432534753d86bda04fb0e56ba1d` |
| `CCGeant4-runtime-20261107-elastic-Z26-40.zip` | `38b5f2f859ba06e440ae4a80f6c622823c0d1f52a743d646dfdb4ddbde6b4b9b` |
| `CCGeant4-runtime-20261107-elastic-Z41-55.zip` | `e803dff4502bfa9ccfb5dced17bc881256f83d99cde08a0770ef270cf98aa358` |
| `CCGeant4-runtime-20261107-elastic-Z56-70.zip` | `eb5668bb0a25b4d800067773440145f82b769839d13c88689a39868fa6c9b6a3` |
| `CCGeant4-runtime-20261107-elastic-Z71-82.zip` | `3f919ab7d576508c76e25331dfe41dcc96ba57f8fd76916a319087affb3fa0c3` |

## License

See [LICENSING.md](LICENSING.md) for the MIT project code, Geant4 terms and
numerical input credits. Cite the corresponding physics methods when using
their data.

This product includes software developed by Members of the Geant4 Collaboration (http://cern.ch/geant4).
