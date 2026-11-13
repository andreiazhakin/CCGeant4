# Reference results

These results use Geant4 11.4.0 with the CCGeant4 sources distributed in this
package and the supplied tables. Each case contains 100,000 events, with
seed 6111401 and eight worker threads.

| File | Contents |
| --- | --- |
| `cases.csv` | Beam, target, cut, physics mode and energy-core selection. |
| `summary.csv` | Energy means and variances with statistical errors, angular quantiles and fraction above 5 mrad. |
| `charges.csv` | Discrete exit-charge counts and fractions. |
| `histograms.csv` | Energy and angle densities for plotting the reference curves. |

The Cu cases with native or tabulated single scattering use an iterative
±3σ core for energy loss and deposition separately. All other energy
summaries use complete distributions. Angular summaries and discrete charge
fractions always use all exiting events. Native Geant4 effective charges
are not treated as discrete charge states.

Run the matching case with `run_example.py`, then pass
`--reference reference/summary.csv` to `analyze.py`. The analysis reads the
energy-core selection from the reference. The default runner seed differs
from the reference seed, and comparisons include uncertainties in both samples.
Add `--plot` to overlay the reference histograms.

Histogram columns are `sample,observable,bin_low,bin_high,density`.
Energy axes use MeV and angles use radians. Energy densities are normalized
to the selected energy sample; angle densities use all exiting events.
