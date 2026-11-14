# Ion transport through a slab

This example records primary-ion energy loss, energy deposited in a slab,
outgoing charge and angle. It uses electromagnetic physics only.

## Build

Use Geant4 11.4.0 built with the CCGeant4 patch. From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/patched-install
cmake --build build --config Release --parallel
```

Initialize the environment for that Geant4 installation before running.
The executable is `build/cc_transport` on Linux and
`build/Release/cc_transport.exe` with Visual Studio.

## Run and analyze

Run from this directory so that the application finds `data/` and `macros/`.
These commands simulate Cu, Xe and U at 400 MeV/u in 100 µm of silicon,
with a 25 µm production cut and 10,000 events per case:

```sh
./build/cc_transport --macro macros/cu_400.mac --output results/dirac_cu.csv
./build/cc_transport --macro macros/xe_400.mac --output results/dirac_xe.csv
./build/cc_transport --macro macros/u_400.mac --output results/dirac_u.csv
python analyze.py results/dirac_cu.csv results/dirac_xe.csv results/dirac_u.csv --reference reference/summary.csv --output results/summary.csv
```

Only `analyze.py` requires Python 3.10+. It summarizes the events and compares
them with the supplied reference results. Add `--plot` to save energy and
angle distributions with reference curves. Only plotting requires Matplotlib
(`python -m pip install matplotlib`).

The default mode is Dirac. Add `--model legacy` or `--model emopt4` for
another mode. For reference comparisons, name the output `legacy_cu.csv`,
`emopt4_cu.csv`, and so on to match the selected model and case.
On Windows, use `./build/Release/cc_transport.exe` in these commands.

## Tables and scattering

In the full release ZIP, the tables in `data/` are already unpacked and ready
to read. The example
includes Cu, Xe and U data for silicon and hydrogen (the vacuum material),
the complete low-energy ATIMA table, and Cu–Si and Fe–Si elastic tables.
A Git clone contains the example sources. Copy `data/` from the
[full release ZIP](https://github.com/andreiazhakin/CCGeant4/releases) into this
example directory.
Use `--tables /path/to/tables` to supply another prepared set.
The full tables are in the [table release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107).

Dirac reads `delta/`, `stopping/`, `global.csv` and `atima_low.csv`.
Legacy reads `global.csv`, `legacy_q2.csv` and `atima_low.csv`.
Standard emopt4 needs no CCGeant4 tables. In [main.cc](main.cc), an empty
`G4VModularPhysicsList` registers the selected physics constructor.
The Dirac and Legacy branches set the table paths and scattering choice
through `G4CCPhysicsOptions`. `cc.mode` selects `G4CCMode::Dirac`
or `G4CCMode::Legacy`. The
[main README](../../README.md#connect-dirac-to-an-application) shows how to
connect the same physics to another application.

The scattering choices are `msc` (option-4 multiple scattering), `ss`
(native single scattering) and `tabulated_ss` (tabulated single scattering).
The emopt4 control uses `msc`. Cu–Si and Fe–Si elastic tables are included
for the tabulated-scattering cases:

```sh
./build/cc_transport --scattering tabulated_ss --macro macros/cu_200.mac --output results-ss/dirac_ss_cu.csv
./build/cc_transport --scattering ss --macro macros/cu_200.mac --output results-ss/dirac_ss_native.csv
python analyze.py results-ss/dirac_ss_cu.csv results-ss/dirac_ss_native.csv --reference reference/summary.csv --output results-ss/summary.csv --plot
```

These cases use Cu-63 at 200 MeV/u in 100 µm of silicon. For Fe-56, use
`macros/fe_200.mac` with `--scattering tabulated_ss` and output
`results-ss/dirac_fe_ss.csv`. Both included elastic tables
cover 10–5000 MeV/u. The full bank covers Z = 26–82 on silicon. Its index
must be at `elastic/tables.index.tsv` under the directory passed to `--tables`.
Other pairs and energies use native single scattering.

## Use your own macro

```sh
./build/cc_transport --model dirac --scattering msc --tables /path/to/tables --output results/cu.csv --macro macros/smoke.mac
```

Select `legacy` or `emopt4` with `--model`, and `ss` or `tabulated_ss`
with `--scattering`. Set threads, random seeds and production cut in the
macro before `/run/initialize`:

```text
/run/numberOfThreads 1
/random/setSeeds 6110101 6110102
/run/setCut 25 um
```

The macro sets material and thickness through
`/example/detector/material` and `/example/detector/thickness`, and the beam
through the standard `/gun/` commands. `/gun/energy` is the total kinetic
energy per ion: Cu-63 at 200 MeV/u uses 12600 MeV. `/run/beamOn` sets the
event count. Use one `/run/beamOn` per invocation because each run opens
its output CSV anew.

The supplied `cu_400.mac`, `xe_400.mac`, `u_400.mac`, `cu_200.mac` and
`fe_200.mac` each request 10,000 events. `smoke.mac` runs 100 Cu events
at 400 MeV/u for a quick check.

## Where to edit the application

Headers are in `include/`, implementations in `src/`, and macros in `macros/`.

| File | Purpose |
| --- | --- |
| `main.cc` | Physics mode, scattering, table paths and run manager. |
| `DetectorConstruction.cc` | Slab material and geometry. |
| `PrimaryGeneratorAction.cc` | Particle gun. |
| `ActionInitialization.cc` | Create the actions for each worker. |
| `RunAction.cc` | Open and close the run output. |
| `EventAction.cc` | Collect and save each event. |
| `SteppingAction.cc` | Score energy deposition, ion exit and escaping energy. |

## Results and reference comparison

The event CSV stores energy loss and energy deposition separately, including
secondary-particle deposits in the slab. It also records charge, angle,
whether the primary exited, and the event-energy balance.

`analyze.py` writes means, variances and their statistical uncertainties,
68% and 95% angular quantiles, and discrete charge fractions. The reference
set contains 100,000 events per case. Comparison matches event-file names
and includes uncertainties from both samples. Deviations above four combined
standard errors are marked `REVIEW`. Angular comparisons also include the
fraction above 5 mrad.

The Cu single-scattering references use an iterative ±3σ energy core.
The analysis reads this selection from the reference summary. Without a
reference, `--core-sigma 3` selects it explicitly. Other reference cases
use the complete energy distributions.
