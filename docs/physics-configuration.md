# Configure transport

Replace the electromagnetic constructor of a modular physics list before
initialization. In this example `DATA` denotes the extracted production
directory. Use absolute paths or paths relative to the application's working
directory.

```cpp
#include "FTFP_BERT.hh"
#include "G4CCPhysicsOptions.hh"
#include "G4EmStandardPhysics_option4_CurrentCharge.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

void ConfigurePhysics(G4RunManager& runManager)
{
    G4CCPhysicsOptions cc;
    cc.tables.delta = "DATA/delta";
    cc.tables.stopping = "DATA/stopping";
    cc.tables.global = "DATA/global.csv";
    cc.tables.atimaLow = "DATA/atima_low.csv";
    cc.transition = {30*MeV, 40*MeV, 1500*MeV, 2000*MeV};
    cc.scattering.mode = G4CCScatteringMode::Opt4MultipleScattering;

    auto* physics = new FTFP_BERT;
    physics->ReplacePhysics(new G4EmStandardPhysics_option4_CurrentCharge(cc, 0));
    runManager.SetUserInitialization(physics);
}
```

The four transition coordinates are expressed in `931.494 MeV × (γ − 1)`.
They control how transport blends between the tabulated and fallback
models. For example the lower transition can be set to 100–120 MeV if the
application needs the explicit model only at higher energies. Frozen-shell
ionization tables cover projectile Z = 29–92 from 30 to 2000 MeV in this
coordinate and up to 27 bound electrons. Large electron counts use squared
ionic-charge continuation from the last directly calculated state. Keep the
electron production threshold at or above 0.99 keV for applicable
material–cuts couples. The threshold affects the partition between explicit
electron production and continuous fluctuations.

For tabulated elastic scattering, set the following before creating the
physics constructor:

```cpp
cc.scattering.mode = G4CCScatteringMode::TabulatedSingleScattering;
cc.scattering.tables = "DATA/elastic/tables.index.tsv";
```

The supplied silicon tables span projectile Z = 26–82 and 10–5000 MeV/u in
kinetic energy per nucleon. Other pairs or energies use the native screened
ion Coulomb model. For native single scattering throughout, select
`G4CCScatteringMode::NativeSingleScattering`. In single scattering,
collisions transfer recoil energy explicitly. A recoil above its production
threshold is tracked; below threshold its energy is deposited locally as
nonionizing energy. The hadronic part of the selected physics list controls
nuclear reactions.

Run different configurations in separate application processes. The tables
are shared; each ion track has its own charge state. The included
[silicon example](../examples/transport/README.md) shows the scoring setup for
both primary energy loss and total energy deposited in the target.
