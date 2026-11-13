"""Write a sparse G4ELT001 table for geometrically suppressed Coulomb scattering.

R(x) is a probability of elastic survival at the classical distance of closest
approach. The runtime interpolates R linearly in x = sin²(theta_CM/2).
All angles are CM angles, energies are kinetic MeV/nucleon and lengths are fm.
The atomic screening length applies to incident ions with Z >= 2.
"""

import argparse
import json
import math
import tempfile
from pathlib import Path

import numpy as np

from convert_tables import HEADER, HEADER_BYTES, KEYS, MAGIC, integral_t


ALPHA = 7.2973525693e-3
HBARC = 197.32698046
BOHR = 52917.7210903
AMU = 931.49410242
TF = 0.5 * (3 * math.pi / 4) ** (2 / 3)


def kinematics(energy, projectile, target):
    m, M = projectile["mass_MeV"], target["mass_MeV"]
    T = energy * projectile["A"]
    s = (m + M) ** 2 + 2 * M * T
    pcm = M * math.sqrt(T * (T + 2 * m) / s)
    reduced = m * M / math.sqrt(s)
    beta = pcm / math.hypot(pcm, reduced)
    zp, zt = projectile["Z"], target["Z"]
    eta = ALPHA * zp * zt / beta
    a = HBARC * eta / pcm
    atomic_radius = TF * BOHR / (zp ** .23 + zt ** .23)
    screening = (HBARC / (2 * pcm * atomic_radius)) ** 2 * (1.13 + 3.76 * eta ** 2)
    return pcm, a, screening, 10 * math.pi * a ** 2


def energy_grid(anchors, max_ratio):
    anchors = np.asarray(anchors, dtype=float)
    if len(anchors) < 2 or np.any(~np.isfinite(anchors)) or np.any(anchors <= 0):
        raise ValueError("at least two positive, finite energies are required")
    if not np.all(np.diff(anchors) > 0):
        raise ValueError("energies must be strictly increasing")
    if not math.isfinite(max_ratio) or max_ratio <= 1:
        raise ValueError("max-energy-ratio must exceed one")
    counts = [max(1, math.ceil((math.log(b) - math.log(a)) / math.log(max_ratio)))
              for a, b in zip(anchors, anchors[1:])]
    if 1 + sum(counts) > 10000:
        raise ValueError("energy settings would produce more than 10,000 nodes")
    pieces = [np.geomspace(a, b, n + 1)[:-1] for a, b, n in zip(anchors, anchors[1:], counts)]
    return np.concatenate([*pieces, anchors[-1:]])


def angular_profile(a, radius, profile, width, n):
    if profile == "sharp":
        if radius <= 2 * a:
            x = np.array([0., .5, .99, .9999, 1.])
            return x, np.ones_like(x)
        edge = (a / (radius - a)) ** 2
        if not 0 < edge < 1:
            raise ValueError("grazing edge cannot be represented at this precision")
        abscissae = np.unique([0., .5 * edge, .99 * edge, .9999 * edge,
                               edge, np.nextafter(edge, 1.), 1.])
        return abscissae, (abscissae <= edge).astype(float)

    radii = radius + width * np.linspace(-12., 12., n + 1)
    local_x = (a / (radii[radii > 2 * a] - a)) ** 2
    base_x = np.geomspace(1e-14, 1., n)
    x = np.unique(np.clip(np.r_[0., local_x, base_x], 0., 1.))
    z = (a * (1 + 1 / np.sqrt(x[1:])) - radius) / width
    t = np.exp(-np.abs(z))
    r = np.r_[1., np.where(z >= 0, 1 / (1 + t), t / (1 + t))]
    return x, r


def build_node(energy, projectile, target, profile, width, angular_nodes):
    pcm, a, screening, C = kinematics(energy, projectile, target)
    radius = projectile["geometric_radius_fm"] + target["geometric_radius_fm"]
    x, r = angular_profile(a, radius, profile, width, angular_nodes)
    dx = np.diff(x)
    ax = screening + x[:-1]
    cells = C * (r[:-1] * dx / (ax * (ax + dx)) + np.diff(r) / dx * integral_t(dx / ax))
    if not np.all(np.isfinite(cells)) or np.any(cells < 0) or x[0] != 0 or x[-1] != 1:
        raise ValueError(f"invalid angular cells at {energy:g} MeV/u")
    if np.any((cells == 0) & ((r[1:] > 0) | (r[:-1] > 0))):
        raise ValueError(f"a nonzero angular cell underflowed at {energy:g} MeV/u")
    sigma = math.fsum(cells)
    if not sigma > 0 or not math.isfinite(sigma):
        raise ValueError(f"invalid integrated cross section at {energy:g} MeV/u")
    cdf = np.r_[0., np.minimum(np.cumsum(cells) / sigma, 1.)]
    tail = np.r_[np.minimum(np.cumsum(cells[::-1])[::-1] / sigma, 1.), 0.]
    cdf[-1] = tail[0] = 1.
    if (np.any(np.diff(cdf) < 0) or np.any(np.diff(tail) > 0)
            or np.max(abs(cdf + tail - 1)) > 1e-12):
        raise ValueError(f"inconsistent cumulative probabilities at {energy:g} MeV/u")
    return dict(energies_MeVu=energy, sigma_elastic_mb=sigma, As=screening,
                p_cm_MeV=pcm, rutherford_C_mb=C, x=x, R=r, CDF=cdf,
                survival=tail, cell_mb=np.r_[cells, 0.])


def isotope(pair, mass, rms, radius_coefficient):
    Z, A = pair
    if not 1 <= Z <= A <= 1000:
        raise ValueError("expected an isotope with 1 <= Z <= A <= 1000")
    estimated = mass is None
    mass = A * AMU if estimated else mass
    if not math.isfinite(mass) or mass <= 0:
        raise ValueError("nuclear masses must be finite and positive")
    if rms is not None and (not math.isfinite(rms) or rms <= 0):
        raise ValueError("charge RMS radii must be finite and positive")
    radius = math.sqrt(5 / 3) * rms if rms is not None else radius_coefficient * A ** (1 / 3)
    return dict(Z=Z, A=A, mass_MeV=mass, mass_source="A * amu" if estimated else "input",
                charge_rms_fm=rms, geometric_radius_fm=radius,
                radius_source="sqrt(5/3) * input RMS" if rms is not None else "r0 * A^(1/3)")


def write_pair(root, projectile, target, nodes, settings):
    ids = (projectile["Z"], projectile["A"], target["Z"], target["A"])
    name = "Z%dA%d_Z%dA%d" % ids
    arrays = {}
    for key in KEYS:
        if key == "offsets":
            arrays[key] = np.cumsum([0] + [len(node["x"]) for node in nodes], dtype="<i8")
        elif np.isscalar(nodes[0][key]):
            arrays[key] = np.array([node[key] for node in nodes], dtype="<f8")
        else:
            arrays[key] = np.concatenate([node[key] for node in nodes]).astype("<f8", copy=False)
    packed = root / (name + ".npz")
    np.savez(packed, **arrays, projectile_Z=ids[0], projectile_A=ids[1],
             target_Z=ids[2], target_A=ids[3], target_mass_MeV=target["mass_MeV"])
    meta = root / (name + ".json")
    report = dict(model="Screened Rutherford with a classical closest-approach survival factor",
                  projectile=projectile, target=target, **settings,
                  constants=dict(alpha=ALPHA, hbarc_MeV_fm=HBARC, bohr_fm=BOHR, amu_MeV=AMU),
                  energy_nodes=len(nodes), angular_points=len(arrays["x"]))
    meta.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="ascii", newline="\n")
    table = root / (name + ".eltb")
    header = HEADER.pack(MAGIC, 1, HEADER_BYTES, *ids, len(nodes), len(arrays["x"]),
                         target["mass_MeV"], projectile["mass_MeV"])
    with table.open("xb") as out:
        out.write(header.ljust(HEADER_BYTES, b"\0"))
        for key in KEYS:
            arrays[key].tofile(out)
    with table.open("rb") as check:
        if check.read(HEADER_BYTES) != header.ljust(HEADER_BYTES, b"\0"):
            raise ValueError("table header verification failed")
        for key in KEYS:
            if check.read(arrays[key].nbytes) != arrays[key].tobytes():
                raise ValueError(f"table verification failed for {key}")
        if check.read(1):
            raise ValueError("unexpected trailing bytes in table")
    row = "\t".join(map(str, [table.name, *ids]))
    (root / "tables.index.tsv").write_text("# G4ELASTIC_TABLE_INDEX 1 1\n" + row + "\n",
                                            encoding="ascii", newline="\n")
    return table.name


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True, help="new directory for one pair; must not exist")
    for side in ("projectile", "target"):
        parser.add_argument("--" + side, nargs=2, type=int, required=True, metavar=("Z", "A"))
        parser.add_argument(f"--{side}-mass-mev", type=float, help="nuclear rest mass in MeV (default: A * amu)")
        parser.add_argument(f"--{side}-rms-fm", type=float, help="charge RMS radius in fm; use sqrt(5/3) times this for contact")
    parser.add_argument("--radius-coefficient-fm", type=float, default=1.2,
                        help="geometric r0 in r0*A^(1/3) when no RMS radius is given (default: 1.2 fm)")
    parser.add_argument("--profile", choices=("smooth", "sharp"), default="smooth")
    parser.add_argument("--width-fm", type=float, default=0.65, help="smooth logistic length in fm (default: 0.65)")
    parser.add_argument("--energies", type=float, nargs="+", default=[100, 300, 1000, 4000], metavar="MEV/U",
                        help="exact energy anchors and bounds (default: 100 300 1000 4000)")
    parser.add_argument("--max-energy-ratio", type=float, default=1.10, help="maximum ratio of adjacent energies (default: 1.10)")
    parser.add_argument("--angular-nodes", type=int, default=32,
                        help="grid resolution for each of two smooth angular grids (default: 32)")
    args = parser.parse_args()
    try:
        if not math.isfinite(args.width_fm) or args.width_fm <= 0:
            raise ValueError("width must be positive and finite")
        if not math.isfinite(args.radius_coefficient_fm) or args.radius_coefficient_fm <= 0:
            raise ValueError("radius coefficient must be positive and finite")
        if not 8 <= args.angular_nodes <= 4096:
            raise ValueError("angular-nodes must be between 8 and 4096")
        root = args.output.absolute()
        if root.exists():
            raise ValueError("output directory already exists; choose a new directory")
        if args.projectile[0] < 2:
            raise ValueError("this screening model requires a projectile with Z >= 2")
        projectile = isotope(args.projectile, args.projectile_mass_mev, args.projectile_rms_fm, args.radius_coefficient_fm)
        target = isotope(args.target, args.target_mass_mev, args.target_rms_fm, args.radius_coefficient_fm)
        energies = energy_grid(args.energies, args.max_energy_ratio)
        nodes = [build_node(e, projectile, target, args.profile, args.width_fm, args.angular_nodes) for e in energies]
        settings = dict(profile=args.profile, width_fm=args.width_fm if args.profile == "smooth" else None,
                        radius_coefficient_fm=args.radius_coefficient_fm, energy_anchors_MeVu=args.energies,
                        max_energy_ratio=args.max_energy_ratio, angular_resolution=args.angular_nodes if args.profile == "smooth" else None)
        root.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=root.parent) as stage:
            product = Path(stage) / "pair"
            product.mkdir()
            name = write_pair(product, projectile, target, nodes, settings)
            product.rename(root)
        print(f"{root / name}  ({len(nodes)} energies, {sum(len(n['x']) for n in nodes)} angular nodes)")
    except (OSError, ValueError, OverflowError) as exc:
        parser.exit(2, f"{parser.prog}: {exc}\n")


if __name__ == "__main__":
    main()
