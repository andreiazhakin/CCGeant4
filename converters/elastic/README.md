# Convert the calculated silicon elastic tables

`convert_tables.py` converts CCGeant4's calculated silicon scattering distributions from NumPy archives to the binary format used by `G4CCElasticTableData`. Each projectile has one table with 207 energy nodes spanning 10–5000 MeV/u. This exporter handles the 57 reference projectiles, Fe through Pb, with a silicon-28 target. Other physical models and target materials require their own input tables.

## Rebuild from the packed data

Install Python 3 and NumPy. The input is the project data archive `Si28_Fe-Pb_dense_R_CDF.zip` (4,587,609,717 bytes). This archive is separate from the source repository. With an existing copy, run from the repository root:

```sh
python -m zipfile -e Si28_Fe-Pb_dense_R_CDF.zip silicon-packed
python converters/elastic/convert_tables.py --source silicon-packed --output silicon-runtime
```

Use an empty output directory and leave Python optimizations disabled so the exporter's numerical assertions run. The input root must contain `packed/*_Si28.npz` and matching JSON files. The outputs include `tables.index.tsv` and 57 `.eltb` files totalling 6,251,348,960 bytes. Keep the index and binaries together. The transport reader uses the filename, projectile Z and A, and target Z and A in each index row. Applications select these data through:

```cpp
cc.scattering.mode = G4CCScatteringMode::TabulatedSingleScattering;
cc.scattering.tables = "/path/to/silicon-runtime/tables.index.tsv";
```

The distributed tables can be used directly. The reader also accepts the additional columns in the original silicon index.

## Checks

The converter checks cell integrals and the normalization of the angular distributions. It rereads the written arrays and compares them with the calculated inputs.

For the physics and upstream sources see [electron and nuclear scattering data](../../docs/dirac-elastic-data.md).

## Generate a sparse Coulomb table for one isotope pair

For pairs without a dedicated elastic table, `generate_grazing_table.py` estimates large-angle suppression by smoothly reducing the screened Coulomb cross section as the distance of closest approach nears the sum of the nuclear radii. This approximation omits Coulomb–nuclear interference and can underestimate the large-angle tail.

For Xe-132 on Si-28, run from the repository root with Python 3.10+ and NumPy:

```sh
python converters/elastic/generate_grazing_table.py --projectile 54 132 --target 14 28 --projectile-rms-fm 4.7859 --target-rms-fm 3.1224 --output xe-si-smooth
```

Set `cc.scattering.tables` to `xe-si-smooth/tables.index.tsv` in `TabulatedSingleScattering` mode.
