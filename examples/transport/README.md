# Ion transport through a slab

This example sends an ion beam through a thin slab and records the primary-ion
energy loss, energy deposited in the slab, outgoing charge and angle. It uses
electromagnetic physics only. Nuclear reactions are not registered.

The application and its result checks are separate. Running the C++ executable
produces events. The Python analysis summarizes the events and can compare
them with a reference summary. Plotting is optional.

## Where to edit the example

| File | Purpose |
| --- | --- |
| `main.cc` | Read the command-line options and create the run manager. |
| `DetectorConstruction.cc` | Define the slab material, thickness and geometry. |
| `PhysicsList.cc` | Select Dirac, Legacy or standard emopt4 and connect the tables. |
| `PrimaryGeneratorAction.cc` | Define the particle gun. |
| `ActionInitialization.cc` | Create the run, event and stepping actions. |
| `RunAction.cc` | Open and close the output for a run. |
| `EventAction.cc` | Collect and save one event. |
| `SteppingAction.cc` | Score energy deposition, primary-ion exit and energy escaping the world. |

Headers are in `include/`, implementations in `src/`. Beam and target settings
are in ordinary Geant4 macro files under `macros/`.

## Build

Use a Geant4 11.4.0 installation built with the CCGeant4 patch.
From this directory on Linux:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/patched-install
cmake --build build --parallel
```

With Visual Studio, use the same CMake commands and select the Release build:

```powershell
cmake -S . -B build -DCMAKE_PREFIX_PATH=C:/path/to/patched-install
cmake --build build --config Release --parallel
```

The executable is `build/cc_transport` on Linux and
`build/Release/cc_transport.exe` with Visual Studio. Initialize the Geant4
environment before running it.

## Run and analyze

The runner unpacks the small supplied Cu, Xe and U table set, then runs
the three ions at 400 MeV/u in 100 µm of silicon with a 25 µm production cut:

```sh
python run_example.py --executable build/cc_transport --output results --events 1000
python analyze.py results/dirac_cu.csv results/dirac_xe.csv results/dirac_u.csv --reference reference/summary.csv --output results/summary.csv
```

Use `--case cu` to run only copper, or `--model legacy` or `--model emopt4`
to select another model. `--tables /path/to/DATA` uses an existing table
directory instead of unpacking the supplied files. On Windows, pass the
Release executable above; `--geant4-bin` can specify the directory containing
the matching Geant4 DLLs.

The small bundled table set covers Cu, Xe and U in silicon and hydrogen
(the vacuum material). For other materials or projectiles, supply the
corresponding tables with `--tables`.
The low-energy ATIMA file is supplied in full because the Legacy loader
requires its complete projectile--target grid.

For a tabulated/native scattering comparison, supply the elastic tables:

```sh
python run_example.py --executable build/cc_transport --case ss_cu ss_native --elastic-tables /path/to/DATA/elastic/tables.index.tsv --output results-ss --events 10000
```

These two cases use Cu-63 at 200 MeV/u in 100 µm of silicon. The additional
`fe_ss` case uses Fe-56 under the same conditions.

## Tables and physics

`--tables` points to a directory with the following contents:

```text
DATA/
  delta/
  stopping/
  global.csv
  atima_low.csv
  legacy_q2.csv
  elastic/tables.index.tsv
  elastic/*.eltb
```

Dirac uses `delta/`, `stopping/`, `global.csv` and `atima_low.csv`.
Legacy uses the three CSV files. The example sets the Legacy table-path
environment variables before creating the physics constructor. Its default
binding-energy override remains 1000 GeV.
Standard emopt4 does not require CCGeant4 tables.

The scattering choices are `msc` (option-4 multiple scattering), `ss`
(native ion single scattering), and `tabulated_ss` (tabulated single
scattering). The emopt4 control supports `msc`.

Only `tabulated_ss` needs the elastic files. By default it reads the index
under `DATA/elastic/`; `--elastic-tables` can point to another index or
directory. The supplied tables cover projectile Z = 26–82 on silicon at
10–5000 MeV/u. Other pairs and energies use native single scattering.
An ion can therefore use tabulated scattering even when its ionization is
outside the Dirac coverage, as for iron. Conversely, uranium uses native
single scattering with the supplied elastic tables.
The elastic tables are downloaded separately; they are not included in
this example's small data directory.

## Run a macro

```sh
./build/cc_transport --model dirac --scattering msc --tables /path/to/DATA --cut-um 25 --threads 1 --seed 12345 --output results/cu.csv --macro macros/smoke.mac
```

Replace `--model dirac` with `legacy` or `emopt4` to select another model.
For an elastic-table comparison, run the same macro with `--scattering ss`
and `--scattering tabulated_ss` into separate output CSV files.

The macro sets the material and thickness with
`/example/detector/material` and `/example/detector/thickness`, and the beam
with the standard `/gun/` commands. `/gun/energy` is the **total kinetic
energy per ion**, not MeV per nucleon. For example, Cu-63 at 200 MeV/u uses
12600 MeV. `/run/beamOn` sets the event count.
Use one `/run/beamOn` per invocation: each run opens a new output CSV.

## Results

The event CSV contains the primary-ion energy loss and the slab energy
deposit separately. Charge, energy and angle are recorded at the primary's
last step in the slab. `exited` distinguishes ions that left the slab from
ions that stopped in it; exits through any face are included. Energy deposited
by secondary particles is included in the slab deposit. Escaping energy is
retained for the event-energy balance.

`analyze.py` writes `summary.csv` with the means, variances, their statistical
uncertainties and the 68% and 95% angular quantiles. It also writes
`charges.csv` with discrete final charge fractions. Samples containing
noninteger effective charges, as in native Geant4 ionization, are omitted
from this discrete-state table. The analysis reports malformed
event data or an energy-balance residual above 10⁻⁶ MeV.

To compare with a saved summary, add `--reference reference/summary.csv`.
Matching uses the event-file names. The comparison checks energy means and
variances, the fractions below the reference 68% and 95% angular quantiles,
the fraction above 5 mrad, and the charge fractions in the adjacent
`charges.csv`, when supplied.
Differences above four combined standard errors are marked `REVIEW`.
Uncertainties from both samples are included; this flag does not stop a run.

The Cu single-scattering references use an iterative ±3σ energy core.
For each energy observable the analysis repeatedly selects values within
three standard deviations of the mean until the selection stops changing.
The reported mean and width use that core. Angles and charge fractions
always use all exiting events. The reference summary records the selection;
without a reference, `--core-sigma 3` selects it explicitly. The other
reference cases use the complete energy distributions.

Add `--plot` to save energy and angle histograms, including reference curves
from `reference/histograms.csv` and the angular tail. Only this option requires
Matplotlib (`python -m pip install matplotlib`); running and numerical
analysis use the Python standard library.

Each physics configuration runs in its own application process. The same
seed makes a run repeatable, but comparisons across builds or platforms use
statistical uncertainties rather than requiring identical events.

For scattering, compare the outgoing angular distribution with the
reference for the selected scattering model. Tabulated and native scattering
are different physical descriptions and need not give the same distribution.
Use the complete exiting sample for this comparison; do not remove the
largest-angle events.
