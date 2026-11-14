# Geant4 source notices

This product includes software developed by Members of the Geant4 Collaboration (http://cern.ch/geant4).

Copyright (c) Copyright Holders of the Geant4 Collaboration, 1994–2006.
Geant4 Software License, version 1.0, 28 June 2006. When redistributing the
patch or the readable source changes, include the full licence and preserve
the source notices and acknowledgement. The MIT grant for independently
authored CCGeant4 work is described in the [licensing statement](../LICENSING.md).
[Bundled Geant4 licence](../LICENSES/Geant4.txt) ·
[Licence and copyright holders](https://geant4.web.cern.ch/download/license.html)

## Code from Geant4 11.4.0

The patch and `source_changes/` contain the same 11 modified upstream files
and 53 added paths. The [file manifest](../patches/BASE.json) lists
these paths. Added paths can also contain
adaptations of upstream code. The examples below identify the origins of
specific parts of the implementation.

| CCGeant4 file or component | Upstream contribution |
| --- | --- |
| `G4CCNuclearStopping.{hh,cc}` | Adapted from `G4NuclearStopping` by Vladimir Ivanchenko. Both files retain the Geant4 headers. The CC class changes the handling of the local nuclear energy deposit. |
| `G4CCLindhardSorensenIonModel.cc` | Adapts the Geant4 ion model by Alexander Bagulya and Vladimir Ivanchenko. This includes projectile parameter setup and parts of `SampleSecondaries`, including energy and direction sampling and the primary momentum update. |
| `G4CCBetheBlochModel.cc` and `G4CCIonDeltaElectronProcess.cc` | Use portions of the projectile form-factor calculation in `G4BetheBlochModel`, credited to Vladimir Ivanchenko, based on code by Laszlo Urban. The CC Bethe–Bloch branch also adapts the restricted stopping calculation. |
| `G4CCTabulatedIonScatteringModel.cc` | Adapts the target-selection and recoil treatment of `G4IonCoulombScatteringModel`. It delegates the fallback cross section and angular law to Geant4. The upstream scattering model and `G4IonCoulombCrossSection` credit Cristina Consolandi. |
| Direct modifications to Geant4 | The original headers credit Laszlo Urban for `G4UrbanMscModel`, Vladimir Ivanchenko for `G4EmCalculator` and `G4PhysListFactory`, Vladimir Ivanchenko based on Laszlo Urban's code for `G4VEnergyLossProcess`, and R. Hatcher for `G4RegisterPhysLists.icc`. The modified `G4VProcess.hh` credits G. Cosmo and H. Kurashige; `G4HadronicProcess.{hh,cc}` credit H. P. Wellisch and J. L. Chuma. The other three modified paths are upstream CMake files. |

[Upstream electromagnetic models](https://gitlab.cern.ch/geant4/geant4/-/tree/v11.4.0/source/processes/electromagnetic/standard) ·
[Upstream electromagnetic utilities](https://gitlab.cern.ch/geant4/geant4/-/tree/v11.4.0/source/processes/electromagnetic/utils) ·
[Upstream physics lists](https://gitlab.cern.ch/geant4/geant4/-/tree/v11.4.0/source/physics_lists/lists)

## Geant4 components used at build time

`G4CCMottCorrection.cc` reads `fMottCoef` from the installed
`G4MottData.hh`. The upstream header credits Vladimir Ivanchenko for the
data class and M. Tacconi for recalculating the fit coefficients. It cites
M. J. Boschini et al. and the fitting procedure of T. Lijian et al.
The coefficients are supplied by the Geant4 installation.
[Upstream data header and references](https://gitlab.cern.ch/geant4/geant4/-/blob/v11.4.0/source/processes/electromagnetic/standard/include/G4MottData.hh)

The CC physics-list constructors extend `FTFP_BERT` (Gunter Folger) and call
`G4EmStandardPhysics_option4::ConstructProcess()` (V. Ivanchenko). The
Geant4 installation supplies these constructors, their implementations and
the standard models selected by them.
[Reference list](https://gitlab.cern.ch/geant4/geant4/-/blob/v11.4.0/source/physics_lists/lists/src/FTFP_BERT.cc) ·
[Option 4 constructor](https://gitlab.cern.ch/geant4/geant4/-/blob/v11.4.0/source/physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option4.cc)
