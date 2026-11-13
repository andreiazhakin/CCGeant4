# Charge-changing cross sections and equilibrium stopping data

CCGeant4 consumes numerical transition cross sections and stopping powers. The
GLOBAL and ATIMA programs and the upstream ATIMA spline archive are obtained
from their respective distributors. This document identifies the inputs used
to prepare the production tables and their attribution.

## GLOBAL

The production input contains 26,143,416 transition cross sections in cm² for
39 nominal energy points between 30 and 2000 MeV/u. The CCGeant4 Dirac
conversion selects 24,605,568 rows for projectile charges 29 through 92.
The table records the nominal energy, the energy supplied to GLOBAL, the
number of bound projectile electrons before and after a transition, and the
run identifier `GLOBAL_DLL_CROSS_20260914`. Preserve both energy columns when
preparing an input for the converter.

| Production input | Bytes |
| --- | ---: |
| `global_transition_sigmas_Z29_96_Zt1_84_88_96_E30_2000_grid39.csv` | 3,143,713,752 |

GSI credits Walter E. Meyerhof for GLOBAL, Bertram Blank for the Fortran
implementation and Oleg Tarasov for the C++ implementation. Cite
C. Scheidenberger, Th. Stöhlker, W. E. Meyerhof, H. Geissel, P. H. Mokler
and B. Blank, *Charge states of relativistic heavy ions in matter*,
Nucl. Instrum. Methods B **142** (1998) 441–462.
[Paper](https://doi.org/10.1016/S0168-583X(98)00244-4) ·
[GSI GLOBAL page and versions](https://web-docs.gsi.de/~weick/charge_states/)

GSI distributes the May 2017 Fortran implementation. FRIB distributes the
LISE++ implementation under its user agreement. The agreement requires public
reports and papers using its results to cite
O. B. Tarasov and D. Bazin, Nucl. Instrum. Methods B **266** (2008) 4657–4664,
and the LISE++ website.
[GSI Fortran archive](https://web-docs.gsi.de/~weick/charge_states/global2017.tar.gz) ·
[FRIB downloads](https://lise.frib.msu.edu/download/) ·
[LISE++ agreement](https://lise.frib.msu.edu/doc/License.pdf) ·
[Tarasov and Bazin](https://doi.org/10.1016/j.nimb.2008.05.110) ·
[LISE++](https://lise.frib.msu.edu/)

## ATIMA 1.41

GSI credits Hans Geissel, Christoph Scheidenberger, Peter Malzacher,
Jörg Kunzendorf and Helmut Weick, and asks users to cite the authors and the
ATIMA website. Obtain version **1.41** from the link for the full spline set
on that page. Inside the distribution, retain the tar of individually gzipped
splines in the form expected by the CCGeant4 converter. The input contains
14,280 compressed data members.
[ATIMA](https://web-docs.gsi.de/~weick/atima/) ·
[Authors' installation notes](https://web-docs.gsi.de/~weick/atima/atima_readme.txt)

| Production file | Role |
| --- | --- |
| `splines_gz_1.41.tar` (1,571,000,320 bytes) | GSI range and related spline coefficients |
| `atima141_equilibrium_dedx_Z29_96_Zt1_84_88_96_E2_40.csv` | Project low-energy stopping values with both one-sided limits |
| `atima_equilibrium_per_atom_reference.csv` | Same values with the target molar mass from the ATIMA spline header |

The production converter selects 5,952 elemental splines in the Dirac
domain. It derives the equilibrium mass stopping power from the reciprocal
derivative of the range with respect to energy per nucleon, retaining both
limits at a jump. The per-atom output uses the target molar mass from the
corresponding header. Thus the output contains evaluated physical quantities
and reference masses. It omits the upstream range coefficients. The output
mesh includes source knot energies to represent the jumps in stopping at the
correct physical energies.

The low-energy production export contains 674,856 rows for 6,324 pairs at
2–40 MeV/u. Its molar-mass column comes from the ATIMA spline headers.
For every pair, the output grid contains the interior range-spline knot
energies and 62 common non-knot energies.

## Distribution scope

The distributed tables contain evaluated cross sections and stopping powers,
with the attribution above. CCGeant4's MIT grant covers its independently
authored code. It grants no rights to either upstream program, coefficient
archive or third-party database. Obtain those inputs and their terms from
the respective distributors.
