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

The two scripts require Python 3.10+. The first runs Cu, Xe and U at
400 MeV/u in 100 µm of silicon with a 25 µm production cut. The second
summarizes the events and compares them with the supplied reference results.

```sh
python run_example.py --executable build/cc_transport --output results --events 1000
python analyze.py results/dirac_cu.csv results/dirac_xe.csv results/dirac_u.csv --reference reference/summary.csv --output results/summary.csv
```

Add `--plot` to the second command to save energy and angle distributions
with reference curves. Only plotting requires Matplotlib
(`python -m pip install matplotlib`).

The runner uses Dirac by default. Select `--model legacy` or
`--model emopt4` for another mode, `--case cu` for one projectile,
`--threads 8` for eight threads, and `--seed 12345` to set the random seed.
On Windows, pass the Release executable above; `--geant4-bin` can specify
the directory containing the matching Geant4 DLLs.

## Tables and scattering

The runner unpacks the tables included in `data/` into `results/tables/`.
This small set contains Cu, Xe and U data for silicon and hydrogen
(the vacuum material). Use `--tables /path/to/DATA` to supply other tables.
The full tables are in the [runtime release](https://github.com/andreiazhakin/CCGeant4/releases/tag/runtime-20261107).

Dirac reads `delta/`, `stopping/`, `global.csv` and `atima_low.csv`.
Legacy reads `global.csv`, `legacy_q2.csv` and `atima_low.csv`.
Standard emopt4 needs no CCGeant4 tables. The paths are passed to the physics
constructor in [PhysicsList.cc](src/PhysicsList.cc); the
[main README](../../README.md#connect-dirac-to-an-application) shows how to
connect the same physics to another application.

The scattering choices are `msc` (option-4 multiple scattering), `ss`
(native single scattering) and `tabulated_ss` (tabulated single scattering).
The emopt4 control uses `msc`. Tabulated scattering needs the separately
supplied elastic tables:

```sh
python run_example.py --executable build/cc_transport --case ss_cu ss_native --elastic-tables /path/to/DATA/elastic/tables.index.tsv --output results-ss --events 10000
python analyze.py results-ss/dirac_ss_cu.csv results-ss/dirac_ss_native.csv --reference reference/summary.csv --output results-ss/summary.csv --plot
```

These cases use Cu-63 at 200 MeV/u in 100 µm of silicon. Add `fe_ss` to
`--case` for Fe-56 under the same conditions. The supplied elastic tables
cover Z = 26–82 on silicon at 10–5000 MeV/u; other pairs and energies use
native single scattering.

## Use your own macro

```sh
./build/cc_transport --model dirac --scattering msc --tables /path/to/DATA --cut-um 25 --threads 1 --seed 12345 --output results/cu.csv --macro macros/smoke.mac
```

Select `legacy` or `emopt4` with `--model`, and `ss` or `tabulated_ss`
with `--scattering`. The macro sets material and thickness through
`/example/detector/material` and `/example/detector/thickness`, and the beam
through the standard `/gun/` commands. `/gun/energy` is the total kinetic
energy per ion: Cu-63 at 200 MeV/u uses 12600 MeV. `/run/beamOn` sets the
event count. Use one `/run/beamOn` per invocation because each run opens
its output CSV anew.

## Where to edit the application

Headers are in `include/`, implementations in `src/`, and macros in `macros/`.

| File | Purpose |
| --- | --- |
| `main.cc` | Command-line options and run manager. |
| `DetectorConstruction.cc` | Slab material and geometry. |
| `PhysicsList.cc` | Physics mode, scattering and table paths. |
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
and includes uncertainties from both samples; deviations above four combined
standard errors are marked `REVIEW`. Angular comparisons also include the
fraction above 5 mrad.

The Cu single-scattering references use an iterative ±3σ energy core.
The analysis reads this selection from the reference summary; without a
reference, `--core-sigma 3` selects it explicitly. Other reference cases
use the complete energy distributions.
