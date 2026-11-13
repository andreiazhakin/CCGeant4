# Sources of the electron and nuclear scattering data

The project calculated two distinct products. The Dirac spectra describe electrons scattered by the field of a fast ion and enter the charge-dependent ionization model. The silicon elastic tables describe scattering of two nuclei and are used only when the application selects tabulated single scattering.

## Dirac electron spectra

For a projectile with nuclear charge `Z` and `N` bound electrons, the calculation used spherical, frozen Dirac–Hartree–Fock–Slater charge densities. The projectile was a point nucleus. The scattering potential contained the nuclear and electronic Hartree fields. Dynamic shell polarization and rearrangement were omitted. The partial-wave solver produced the bare differential spectrum and the screened-to-bare ratio. The grid has 39 projectile energy coordinates from 30 to 2000 MeV. Direct calculations cover `N = 0…10` for `Z = 29…50` and `N = 0…15` for `Z = 51…92`. At larger `N`, transport uses squared ionic-charge scaling from the last calculated state.

| Offline input | Source and terms | Use in the calculation |
| --- | --- | --- |
| RADIAL with DHFS, 2018 code, 2019 deposit v1 | F. Salvat and J. M. Fernández-Varea. [Distribution](https://doi.org/10.17632/h6nbws5h6s.1), [paper](https://doi.org/10.1016/j.cpc.2019.02.011). CC BY-NC 3.0. | The project adapted the density builder to specify ionic occupancies, set numerical controls and export self-consistent densities. |
| ELSEPA 2020, deposit v1 | F. Salvat, A. Jablonski and C. J. Powell. [Distribution](https://doi.org/10.17632/w4hm5vymym.1), [method](https://doi.org/10.1016/j.cpc.2004.09.006). CC BY-NC 3.0. | A project driver used modified partial-wave controls and the frozen-ion potentials to calculate the angular differential cross sections and convert them to energy spectra. |

The deposits above provide the upstream solvers. This repository contains the transport converter, which takes calculated spectra as input. It does not contain the drivers for the offline atomic calculations.

The input to that converter is the project-generated `Dirac_ALL_N_Z29-92_2026-11-01.zip` (723,024,033 bytes). It contains 914 charge–projectile combinations. The transport `delta` export consists of 64 CSV files, one per projectile, totalling 1,587,633,267 bytes. Their payload is projectile and energy coordinates and calculated differential cross sections. The separate 64 `stopping` CSV files (202,919,584 bytes) combine these results with GLOBAL and ATIMA inputs. Their provenance therefore also includes [GLOBAL and ATIMA](global-atima-data.md).

**Distribution scope.** The numerical Dirac spectra are project calculation results. The two solver deposits license their published materials under [CC BY-NC 3.0](https://creativecommons.org/licenses/by-nc/3.0/legalcode.en). CCGeant4's MIT grant does not cover those upstream distributions, their bundled databases or copied Fortran code.

## Ion–silicon elastic scattering

The production export covers 57 projectile elements (`Z = 26…82`) on silicon, with 207 energy nodes from 10 to 5000 MeV/u. The model combines nuclear and Coulomb scattering within each branch, interpolates their cross-section ratios across the transition bands, then applies atomic screening. The converter writes the sampled quantities used by transport, including the ratio to Rutherford scattering, energy-dependent total cross section, cell integrals and cumulative probabilities.

| Nominal energy (MeV/u) | Input used in the production branch | Terms and attribution |
| --- | --- | --- |
| 10–50 | REGINA SPP2 real potential with an absorptive component equal to 0.80 times the real potential, solved with FRESCOX. REGINA also supplies the spherical nuclear density input used in the Glauber calculation. | L. C. Chamon, B. V. Carlson and L. R. Gasques, [REGINA deposit v1, GPLv3](https://doi.org/10.17632/vfkkjb8dv7.1), [paper](https://doi.org/10.1016/j.cpc.2021.108061). |
| 50–400 | Furumoto MGOP v1.02, CEG07b; the silicon parameter set supplies the optical potential in the reciprocal orientation. FRESCOX solves the scattering problem. | The [authors' download page](https://www2.yukawa.kyoto-u.ac.jp/~takenori.furumoto/code/mgop/index.html) requests citation of [the global potential](https://doi.org/10.1103/PhysRevC.85.044607) and both underlying interaction papers, [2009](https://doi.org/10.1103/PhysRevC.80.044614) and [2008](https://doi.org/10.1103/PhysRevC.78.044610). |
| 400–1300 | A scalar forward nucleon–nucleon amplitude was constructed from the SAID SM16 observables SGT, ALFA and the forward DSG slope for the optical-limit Glauber calculation. | R. L. Workman, W. J. Briscoe and I. I. Strakovsky, [SM16](https://doi.org/10.1103/PhysRevC.94.065203), [GWU SAID service](https://gwdac.phys.gwu.edu/). |
| 1300–5000 | The same Glauber construction uses the published Uzhinsky forward nucleon–nucleon parametrization. | [Forward-amplitude model](https://doi.org/10.1103/PhysRevC.94.064003), [Glauber and Coulomb framework](https://doi.org/10.1016/j.nuclphysa.2017.02.002). |

The offline optical calculation used Ian Thompson's FRESCOX, revision `eca88f3d3483c4b817051734c698845abcecee5d`. Its `copyright.txt` identifies GPL version 2 or later and gives LLNL's notices. Geant4 11.4 supplies the reference nuclear masses and the atomic-screening convention. [FRESCOX source](https://github.com/LLNL/Frescox/tree/eca88f3d3483c4b817051734c698845abcecee5d) · [Geant4 Physics Reference Manual](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/elastic_scattering/ion_coulomb.html)

**Distribution scope.** The transport files contain calculated scattering arrays, reference isotope identifiers and nuclear masses. The upstream solvers, REGINA density grids, MGOP coefficients and SAID service responses are obtained separately under their distributors' terms.

For the exact input files, packing command and Geant4 configuration see the [silicon converter](../converters/elastic/README.md). Its input is the calculated NPZ dataset. A fresh calculation from nuclear potentials requires the external offline generator and source inputs. The full numerical banks are separate from this source checkout.
