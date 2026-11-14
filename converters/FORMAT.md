# Converter input and output formats

`Zp` and `Zt` are projectile and target atomic numbers. `N` is the number of
bound projectile electrons, so the ionic charge is `q = Zp − N`. The Dirac
spectra depend on the projectile and its charge state. GLOBAL transitions
and stopping coefficients also depend on the target element.

## Dirac inputs

### Spectrum ZIP

The ZIP contains `full_sweep/results/`, optionally under a top-level
directory. The results directory contains:

| File | Required contents |
| --- | --- |
| `grid.npz` | `E_u_MeV_u`: incident energies; `offsets`: boundaries of the concatenated electron grids; `T_eV`: electron kinetic energies in eV |
| `bare_Z{Zp}.npz` | `Z`; `D0_cm2_per_eV`: bare-ion differential cross section per target electron; `valid_energy` |
| `Z{Zp}_N{N}.npz` | `Z`, `N`; `delta_Z`: screening correction on the same grid; `valid_energy` |

Energy slice `i` is `offsets[i]:offsets[i+1]`. The converter reconstructs
`D_N = D0 * (1 − delta_Z/Zp)^2`, retaining negative values of `delta_Z`.
It multiplies cross sections by `1e6` to obtain cm²/MeV and electron energies
by `1e-6` to obtain MeV.

The supplied spectra have 39 incident energies from 30 to 2000 MeV/u.
Direct calculations cover `N = 0…10` for `Zp ≤ 50` and `N = 0…15` above
that. For normalization at larger `N ≤ 27`, the last calculated spectrum
is scaled by `((Zp − N)/(Zp − Nmax))^2`, as in transport.

### GLOBAL CSV

```text
Zp,Ap,Zt,At,E_MeV_u,Qe_in,Qe_out,dQe,kind,sigma_cm2,GLOBAL_E_MeV_u,source_run
```

`Qe_in` and `Qe_out` are bound-electron counts, not ionic charges.
`dQe = Qe_out − Qe_in`; capture increases this count and loss decreases it.
`sigma_cm2` is the cross section per target atom. Each energy node contains
106 transitions with changes of one or two electrons between `N = 0…27`.
`Ap` and `At` are the reference mass parameters.
Rows for each projectile, target and energy must form a contiguous group.

`E_MeV_u` is the nominal table coordinate. `GLOBAL_E_MeV_u` records the
energy supplied to GLOBAL; the nominal 30 MeV/u point was calculated at
30.1 MeV/u. Selected rows, including both coordinates, are copied to the
output `global.csv`.

### ATIMA range splines

`splines_gz_1.41.tar` contains individually gzipped files such as
`Cu63-Si.dat.gz`. Their headers specify projectile Z and mass, target
composition and molar mass. The converter uses the first B-spline block:
the range `R(E)` in mg/cm², with `E` in MeV/u. Other integrated quantities
in the files are not used.

The range derivative gives the mass stopping power
`S_mass = 1000 m_ref / R′(E)` in MeV cm²/g. Multiplication by
`A_target / N_A` gives the stopping cross section per atom, where
`m_ref` and `A_target` are the projectile mass in u and target molar mass
in g/mol recorded by ATIMA.

## Dirac output

The converter writes CSV files directly, without an intermediate binary
format.

`delta/Znnn.csv` has the header:

```text
Zp,Nmax,E_MeV_u,T_MeV,D0_cm2_MeV,D1_cm2_MeV,...,D15_cm2_MeV
```

Rows are ordered by incident energy, then electron energy. Columns above
`Nmax` are zero; transport uses the charge scaling described above instead.
Stored spectra already contain the screening correction. No additional
Mott or charge-squared factor should be applied to them.

`stopping/Znnn.csv` has the header:

```text
Zp,Zt,E_MeV_u,B_MeV_cm2,B_left_MeV_cm2,ATIMA_MeV_cm2,ATIMA_left_MeV_cm2
```

The GLOBAL equilibrium fractions `p_N` determine the soft coefficient:

```text
h_N = Zt * integral(T * D_N(T), T0, Tmax)
B = (ATIMA − sum_N p_N h_N) / sum_N p_N (Zp − N)^2
full stopping per atom = B (Zp − N)^2 + h_N
```

Here `T0 = 0.001 MeV`. Between incident-energy nodes, the two neighbouring
spectra are mapped to the new kinematic interval and mixed with weights
linear in log energy. Between electron-energy nodes, each spectrum follows
a power law. The stopping grid is refined to a relative equilibrium
normalization tolerance of `1e-4`.

Each target has its own stopping grid. The `_left` columns retain the
values immediately below a discontinuity in the derivative of an ATIMA
range spline. Interpolation uses `B` at the lower node and `B_left` at the
upper node, linearly in log incident energy.

The complete supplied `atima_low.csv` has the header:

```text
Zp,Zt,At_ref_g_mol,E_MeV_u,dEdx_MeV_cm2_g,dEdx_left_MeV_cm2_g
```

It contains equilibrium stopping over 2–40 MeV/u. `At_ref_g_mol` is the
target molar mass from ATIMA. This file is copied unchanged for every
Dirac conversion, irrespective of the selected pairs.

## Elastic input and output

Each input `.npz` has a matching `.json`. The NPZ supplies scalar
`projectile_Z`, `projectile_A`, `target_Z`, `target_A`, `target_mass_MeV`
and the arrays below. The JSON supplies `projectile.Z`, `projectile.A`
and `projectile.mass_MeV`.

| Array | Length | Meaning |
| --- | --- | --- |
| `energies_MeVu` | `nE` | Incident energy in MeV/u |
| `offsets` | `nE + 1` | Boundaries of the concatenated angular grids |
| `sigma_elastic_mb` | `nE` | Integrated elastic cross section in mb |
| `As` | `nE` | Atomic-screening parameter |
| `p_cm_MeV` | `nE` | Centre-of-mass momentum in MeV/c |
| `rutherford_C_mb` | `nE` | Rutherford normalization in mb |
| `x` | `nP` | Angular coordinate `sin²(theta_CM/2)` |
| `R` | `nP` | Cross-section ratio to Rutherford |
| `CDF`, `survival` | `nP` | Cumulative probability and its complement |
| `cell_mb` | `nP` | Integrated cross section in each angular cell; final entry of each energy slice is zero |

The transport measure is `dσ/dx = C R(x)/(As + x)^2`, with `R` linear
between angular nodes. Each energy slice spans `x = 0…1`.

An `.eltb` starts with a 256-byte little-endian header. Its fields are
`<8sIIiiiiQQdd`: magic `G4ELT001`, version 1, header size 256,
projectile Z and A, target Z and A, `nE`, `nP`, target mass and projectile
mass in MeV/c². The remaining header bytes are zero. The arrays follow in
the order listed above, with `CDF` before `survival`. `offsets` uses int64;
all other arrays use float64.

`tables.index.tsv` begins with `# G4ELASTIC_TABLE_INDEX 1 COUNT` and lists
the filename, projectile Z and A, and target Z and A, separated by tabs.
Only selected pairs are listed. No energy or angular grid is changed.
