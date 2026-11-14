#!/usr/bin/env python3
"""Convert Dirac/GLOBAL/ATIMA or elastic-scattering data to CCGeant4 tables."""
import argparse
from dataclasses import dataclass
from functools import cached_property
import gzip
import io
import json
import lzma
import os
from pathlib import Path
import re
import shutil
import struct
import tarfile
import zipfile

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np


AVOGADRO = 6.02214076e23
N_STATES = 28
SYMBOLS = (
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni '
    'Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe '
    'Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au '
    'Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf '
    'Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og E119 E120'
).split()
SYMBOL_TO_Z = {s: i+1 for i, s in enumerate(SYMBOLS)}
TRANSITIONS = np.array([(n, m) for n in range(N_STATES)
                        for m in range(max(0, n-2), min(N_STATES, n+3))
                        if m != n], dtype=np.int16)


T0 = 0.001
ME_MEV = 0.51099895
MU_MEV = 931.494


def tmax(energy):
    return 2*ME_MEV*((1+np.asarray(energy)/MU_MEV)**2-1)


def bracket(grid, query):
    query = np.atleast_1d(np.asarray(query, dtype=float))
    if query.ndim != 1 or np.any(~np.isfinite(query)) or np.any(query < grid[0]) or np.any(query > grid[-1]):
        raise ValueError('Energy outside table coverage')
    hi = np.clip(np.searchsorted(grid, query, side='right'), 1, len(grid)-1)
    lo = hi-1
    w = np.log(query/grid[lo])/np.log(grid[hi]/grid[lo])
    return query, lo, hi, w


def powerlaw_cells(t, density, order):
    """Integral T**order D(T)dT in each log-log linear cell."""
    log_t = np.log(t[1:]/t[:-1])
    v = np.log(density[..., 1:]/density[..., :-1])+(order+1)*log_t
    factor = np.ones_like(v)
    np.divide(np.expm1(v), v, out=factor, where=np.abs(v)>1e-14)
    return density[..., :-1]*t[:-1]**(order+1)*log_t*factor


def integrate_spectrum(t, density):
    result = []
    for order in range(3):
        cells = powerlaw_cells(t, density, order)
        tail = np.zeros_like(density)
        tail[..., :-1] = np.cumsum(cells[..., ::-1], axis=-1)[..., ::-1]
        result.append(tail)
    return result


def delta_moments(data, energy, zp):
    """Return [query,N=0..27,moment=0..2] in cm2*MeV**moment."""
    e, lo, hi, w = bracket(data['energy_MeV_u'], energy)
    idx = data['offsets'][:-1]
    native = np.stack([data[k][:, idx].T for k in
                       ['tail0_cm2', 'tail1_MeV_cm2', 'tail2_MeV2_cm2']], axis=-1)
    lengths = tmax(data['energy_MeV_u'])-T0
    means_u = np.empty_like(native)
    means_u[..., 0] = native[..., 0]
    means_u[..., 1] = (native[..., 1]-T0*native[..., 0])/lengths[:, None]
    means_u[..., 2] = (native[..., 2]-2*T0*native[..., 1]+T0*T0*native[..., 0])/lengths[:, None]**2
    mixed = (1-w[:, None, None])*means_u[lo]+w[:, None, None]*means_u[hi]
    L = (tmax(e)-T0)[:, None]
    out = np.empty((len(e), 28, 3))
    ncalc = native.shape[1]
    out[:, :ncalc, 0] = mixed[..., 0]
    out[:, :ncalc, 1] = T0*mixed[..., 0]+L*mixed[..., 1]
    out[:, :ncalc, 2] = T0*T0*mixed[..., 0]+2*T0*L*mixed[..., 1]+L*L*mixed[..., 2]
    for n in range(ncalc, 28):
        out[:, n] = out[:, ncalc-1]*((zp-n)/(zp-(ncalc-1)))**2
    return out


def _numbers(line: str) -> np.ndarray:
    return np.fromstring(line.replace('D', 'E').replace('d', 'e'), sep=' ', dtype=np.float64)


def _bspline(knots: np.ndarray, coeff: np.ndarray, degree: int, energy):
    x = np.asarray(energy, dtype=np.float64)
    if not np.all(np.isfinite(x)) or np.any(x < knots[degree]) or np.any(x > knots[len(coeff)]):
        raise ValueError('ATIMA energy outside spline domain; extrapolation disabled')
    indices = np.clip(np.searchsorted(knots, x, side='right')-1, degree, len(coeff)-1)
    values = np.stack([coeff[indices-degree+j] for j in range(degree+1)])
    for level in range(1, degree+1):
        for j in range(degree, level-1, -1):
            low = knots[indices-degree+j]
            den = knots[indices+1+j-level]-low
            fraction = np.divide(x-low, den, out=np.zeros_like(x), where=den != 0)
            values[j] = (1.0-fraction)*values[j-1]+fraction*values[j]
    return values[degree]


@dataclass
class AtimaSpline:
    mass_ref_u: float
    target_mass_g_mol: float
    degree: int
    knots: np.ndarray
    derivative: np.ndarray

    @cached_property
    def discontinuities_MeV_u(self):
        knots, counts = np.unique(self.knots, return_counts=True)
        low, high = self.knots[self.degree], self.knots[len(self.derivative)+1]
        return knots[(counts >= self.degree) & (knots > low) & (knots < high)]

    def atomic_stopping(self, energy, side='right'):
        energy = np.asarray(energy, dtype=float)
        if side == 'left':
            energy = np.where(np.isin(energy, self.discontinuities_MeV_u),
                              np.nextafter(energy, -np.inf), energy)
        slope = _bspline(self.knots[1:-1], self.derivative, self.degree-1, energy)
        if np.any(~np.isfinite(slope)) or np.any(slope <= 0):
            raise ValueError('ATIMA range derivative must be positive')
        return (1000.0*self.mass_ref_u/slope)*self.target_mass_g_mol/AVOGADRO


def read_atima(path, projectiles, targets):
    splines = {}
    pattern = r'(?:(E119|E120|[A-Z][a-z]?)(\d+)|(\d{3})X(\d+))-([A-Z][a-z]?)\.dat\.gz'
    with tarfile.open(path, mode='r|*') as archive:
        for member in archive:
            match = re.fullmatch(pattern, Path(member.name).name)
            if not member.isfile() or match is None:
                continue
            symbol, _, numeric_z, _, target_symbol = match.groups()
            zp = int(numeric_z) if numeric_z else SYMBOL_TO_Z.get(symbol)
            zt = SYMBOL_TO_Z.get(target_symbol)
            if zp not in projectiles or zt not in targets:
                continue
            with archive.extractfile(member) as stream:
                lines = gzip.decompress(stream.read()).decode('ascii').splitlines()
            projectile, components, target = map(_numbers, lines[:3])
            if (len(projectile) != 2 or components.tolist() != [1]
                    or len(target) != 4 or projectile[0] != zp
                    or target[0] != zt or target[1] != 1):
                raise ValueError(f'Invalid ATIMA element header: {member.name}')
            ncoeff, order = map(int, lines[4].split())
            degree = order-1
            knots, coefficients = map(_numbers, lines[5:7])
            if (degree < 1 or ncoeff < order or len(knots) != ncoeff+order
                    or len(coefficients) != ncoeff or np.any(np.diff(knots) < 0)):
                raise ValueError(f'Invalid ATIMA spline: {member.name}')
            if not all(np.isfinite(m) and m > 0 for m in (projectile[1], target[3])):
                raise ValueError(f'Invalid ATIMA masses: {member.name}')
            den = knots[degree+1:degree+ncoeff]-knots[1:ncoeff]
            derivative = np.divide(degree*np.diff(coefficients), den,
                                   out=np.zeros(ncoeff-1), where=den != 0)
            if (zp, zt) in splines:
                raise ValueError(f'More than one ATIMA range curve for {(zp, zt)}')
            splines[zp, zt] = AtimaSpline(float(projectile[1]), float(target[3]),
                                         degree, knots, derivative)
    missing = {(zp, zt) for zp in projectiles for zt in targets} - splines.keys()
    if missing:
        raise ValueError(f'Missing ATIMA pairs: {sorted(missing)}')
    return splines


def read_global(path, projectiles, targets, output):
    header = ['Zp', 'Ap', 'Zt', 'At', 'E_MeV_u', 'Qe_in', 'Qe_out',
              'dQe', 'kind', 'sigma_cm2', 'GLOBAL_E_MeV_u', 'source_run']
    transition_index = {tuple(pair): i for i, pair in enumerate(TRANSITIONS)}
    curves = {}
    current = None
    values = np.full(len(TRANSITIONS), np.nan)

    def finish_group():
        if current is not None:
            if not np.all(np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f'Incomplete or invalid GLOBAL transitions: {current}')
            zp, zt, energy = current
            curves.setdefault((zp, zt), []).append((energy, values.copy()))

    with path.open('rb') as source, output.open('wb') as destination:
        raw_header = next(source)
        if raw_header.decode('utf-8-sig').strip().split(',') != header:
            raise ValueError('Unsupported GLOBAL CSV header')
        destination.write(raw_header)
        for raw in source:
            fields = raw.rstrip(b'\r\n').split(b',')
            zp, zt = int(fields[0]), int(fields[2])
            if zp not in projectiles or (targets is not None and zt not in targets):
                continue
            energy = float(fields[4])
            if not np.isfinite(energy) or energy <= 0:
                raise ValueError('GLOBAL energy must be positive')
            key = zp, zt, energy
            if key != current:
                finish_group()
                values.fill(np.nan)
                current = key
            n, m, dn = map(int, fields[5:8])
            index = transition_index.get((n, m))
            if (index is None or dn != m-n
                    or fields[8] != (b'capture' if dn > 0 else b'loss')):
                raise ValueError(f'Invalid GLOBAL transition: {key}, {n}, {m}')
            if not np.isnan(values[index]):
                raise ValueError(f'Duplicate GLOBAL transition: {key}, {n}, {m}')
            values[index] = float(fields[9])
            destination.write(raw)
        finish_group()
    if targets is None:
        targets = sorted({zt for zp, zt in curves})
    missing = {(zp, zt) for zp in projectiles for zt in targets} - curves.keys()
    if not curves or missing:
        raise ValueError(f'Missing GLOBAL pairs: {sorted(missing)}')
    for key, rows in curves.items():
        rows.sort(key=lambda row: row[0])
        energies = np.array([row[0] for row in rows])
        if len(energies) < 2 or np.any(np.diff(energies) <= 0):
            raise ValueError(f'Invalid GLOBAL energy grid: {key}')
        curves[key] = energies, np.array([row[1] for row in rows])
    return curves, targets


def interpolate_global_rates(lower, upper, weight):
    lower, upper, weight = np.broadcast_arrays(lower, upper, weight)
    positive = (lower > 0) & (upper > 0)
    result = lower + weight*(upper-lower)
    result[positive] = np.exp(np.log(lower[positive]) + weight[positive]
                              * (np.log(upper[positive])-np.log(lower[positive])))
    return result


def equilibrium_from_rates(rates):
    probabilities = np.empty((len(rates), N_STATES))
    diagonal = np.arange(N_STATES)
    for start in range(0, len(rates), 2048):
        block = rates[start:start+2048]
        generator = np.zeros((len(block), N_STATES, N_STATES))
        generator[:, TRANSITIONS[:, 0], TRANSITIONS[:, 1]] = block
        generator[:, diagonal, diagonal] = -generator.sum(axis=2)
        scale = np.max(np.abs(generator), axis=(1, 2))
        if np.any(scale <= 0):
            raise ValueError('GLOBAL rates do not define an equilibrium')
        generator /= scale[:, None, None]
        matrix = generator.swapaxes(-1, -2).copy()
        matrix[:, -1, :] = 1
        rhs = np.zeros((len(block), N_STATES, 1))
        rhs[:, -1, 0] = 1
        p = np.linalg.solve(matrix, rhs)[..., 0]
        if np.any(~np.isfinite(p)) or np.min(p) < -1e-10:
            raise ValueError('Invalid GLOBAL equilibrium probabilities')
        p = np.maximum(p, 0)
        p /= p.sum(axis=1)[:, None]
        probabilities[start:start+len(block)] = p
    return probabilities


class DiracArchive:
    def __init__(self, path):
        self.archive = zipfile.ZipFile(path)
        grids = [n for n in self.archive.namelist()
                 if n.endswith('full_sweep/results/grid.npz')]
        if len(grids) != 1:
            raise ValueError('Expected one full_sweep/results/grid.npz in spectra ZIP')
        self.prefix = grids[0][:-len('grid.npz')]
        self.projectiles = sorted(int(m.group(1)) for name in self.archive.namelist()
                                  if (m := re.fullmatch(re.escape(self.prefix)
                                                       + r'bare_Z(\d+)\.npz', name)))
        self.grid = self.arrays('grid.npz')
        energies, offsets, transfer = (self.grid[k] for k in ('E_u_MeV_u', 'offsets', 'T_eV'))
        if (energies.ndim != 1 or len(energies) < 2 or np.any(np.diff(energies) <= 0)
                or not np.all(np.isfinite(energies)) or offsets.dtype.kind not in 'iu'
                or offsets.shape != (len(energies)+1,) or offsets[0] != 0
                or offsets[-1] != len(transfer) or np.any(np.diff(offsets) < 2)):
            raise ValueError('Invalid Dirac energy grid or offsets')
        for i, energy in enumerate(energies):
            t = transfer[offsets[i]:offsets[i+1]]
            if (t[0] != 1000 or np.any(np.diff(t) <= 0) or not np.all(np.isfinite(t))
                    or abs(t[-1]/1e6/tmax(energy)-1) > 1e-12):
                raise ValueError('Invalid Dirac electron grid or maximum energy')

    def arrays(self, name):
        with np.load(io.BytesIO(self.archive.read(self.prefix+name)), allow_pickle=False) as data:
            return dict(data)

    def projectile(self, zp):
        nmax = 10 if zp <= 50 else 15
        bare = self.arrays(f'bare_Z{zp}.npz')
        d0 = bare['D0_cm2_per_eV']
        if (int(bare['Z']) != zp or not np.all(bare['valid_energy'])
                or d0.shape != self.grid['T_eV'].shape):
            raise ValueError(f'Invalid bare Dirac spectrum: Z={zp}')
        density = []
        for n in range(nmax+1):
            data = self.arrays(f'Z{zp}_N{n}.npz')
            delta = data['delta_Z'].astype(float)
            if (int(data['Z']) != zp or int(data['N']) != n
                    or not np.all(data['valid_energy']) or delta.shape != d0.shape
                    or np.any(delta >= zp) or (n == 0 and np.any(delta != 0))):
                raise ValueError(f'Invalid Dirac spectrum: Z={zp}, N={n}')
            density.append(d0*(1-delta/zp)**2*1e6)
        d = np.array(density)
        if np.any(~np.isfinite(d)) or np.any(d <= 0):
            raise ValueError(f'Nonpositive Dirac spectrum: Z={zp}')
        arrays = {'energy_MeV_u': self.grid['E_u_MeV_u'], 'offsets': self.grid['offsets'],
                  'T_MeV': self.grid['T_eV']/1e6, 'D_cm2_MeV': d}
        tails = [np.empty_like(d) for _ in range(3)]
        for first, last in zip(arrays['offsets'][:-1], arrays['offsets'][1:]):
            sl = slice(first, last)
            for tail, values in zip(tails, integrate_spectrum(arrays['T_MeV'][sl], d[:, sl])):
                tail[:, sl] = values
        arrays.update(zip(('tail0_cm2', 'tail1_MeV_cm2', 'tail2_MeV2_cm2'), tails))
        return arrays


def normalization(zp, zt, curve, spline, delta, energy):
    global_energy, sigma = curve
    energy, lo, hi, w = bracket(global_energy, energy)
    p = equilibrium_from_rates(interpolate_global_rates(sigma[lo], sigma[hi], w[:, None]))
    stopping = spline.atomic_stopping(energy)
    stopping_left = spline.atomic_stopping(energy, side='left')
    hard = zt*delta_moments(delta, energy, zp)[..., 1]
    q2 = (zp-np.arange(N_STATES, dtype=float))**2
    q2eq = np.sum(p*q2, axis=-1)
    heq = np.sum(p*hard, axis=-1)
    soft = (stopping-heq)/q2eq
    soft_left = (stopping_left-heq)/q2eq
    if (np.any(~np.isfinite(soft)) or np.any(soft < 0)
            or np.any(~np.isfinite(soft_left)) or np.any(soft_left < 0)):
        raise ValueError(f'Negative/nonfinite soft stopping: Z={zp}, Zt={zt}')
    return soft, soft_left, stopping, stopping_left, q2eq


def prepare_stopping(zp, zt, curve, spline, delta):
    energy = delta['energy_MeV_u']
    jumps = spline.discontinuities_MeV_u
    jumps = jumps[(jumps > energy[0]) & (jumps < energy[-1])]
    grid = np.unique(np.r_[curve[0], energy, jumps])
    fractions = np.array([.25, .5, .75])
    for _ in range(13):
        b, b_left, stopping, stopping_left, _ = normalization(zp, zt, curve, spline, delta, grid)
        probes = np.exp(np.log(grid[:-1, None])
                        + np.log(grid[1:, None]/grid[:-1, None])*fractions)
        actual, _, actual_stopping, _, q2eq = normalization(zp, zt, curve, spline, delta, probes.ravel())
        interpolated = ((1-fractions)*b[:-1, None] + fractions*b_left[1:, None]).ravel()
        error = np.abs((interpolated-actual)*q2eq/actual_stopping).reshape(-1, 3).max(axis=1)
        if error.max() <= 0.9e-4:
            return np.column_stack((np.full(len(grid), zp), np.full(len(grid), zt),
                                    grid, b, b_left, stopping, stopping_left))
        grid = np.unique(np.r_[grid, np.sqrt(grid[:-1]*grid[1:])[error > 0.9e-4]])
        if len(grid) > 10000:
            break
    raise ValueError(f'Stopping grid refinement failed: Z={zp}, Zt={zt}')


def write_delta(path, zp, delta):
    nmax = delta['D_cm2_MeV'].shape[0]-1
    with path.open('w', encoding='ascii', newline='\n') as output:
        output.write('Zp,Nmax,E_MeV_u,T_MeV,' + ','.join(f'D{n}_cm2_MeV' for n in range(16)) + '\n')
        for i, energy in enumerate(delta['energy_MeV_u']):
            first, last = delta['offsets'][i:i+2]
            rows = np.zeros((last-first, 20))
            rows[:, :3] = zp, nmax, energy
            rows[:, 3] = delta['T_MeV'][first:last]
            rows[:, 4:5+nmax] = delta['D_cm2_MeV'][:, first:last].T
            np.savetxt(output, rows, delimiter=',', fmt='%.17g')


def convert_dirac(args):
    spectra = DiracArchive(args.spectra)
    try:
        projectiles = sorted(set(args.projectiles or spectra.projectiles))
        missing = set(projectiles) - set(spectra.projectiles)
        if not projectiles or missing or not set(projectiles).issubset(range(29, 93)):
            raise ValueError(f'Requested projectile spectra are unavailable or outside Z=29..92: {projectiles}')
        targets = sorted(set(args.targets)) if args.targets else None
        print('Reading GLOBAL and ATIMA inputs...', flush=True)
        global_data, targets = read_global(args.global_csv, projectiles, targets, args.output/'global.csv')
        splines = read_atima(args.atima, projectiles, targets)
        (args.output/'delta').mkdir(exist_ok=True)
        (args.output/'stopping').mkdir(exist_ok=True)
        for zp in projectiles:
            delta = spectra.projectile(zp)
            write_delta(args.output/'delta'/f'Z{zp:03d}.csv', zp, delta)
            with (args.output/'stopping'/f'Z{zp:03d}.csv').open('w', encoding='ascii', newline='\n') as output:
                output.write('Zp,Zt,E_MeV_u,B_MeV_cm2,B_left_MeV_cm2,ATIMA_MeV_cm2,ATIMA_left_MeV_cm2\n')
                for zt in targets:
                    rows = prepare_stopping(zp, zt, global_data[zp, zt], splines[zp, zt], delta)
                    np.savetxt(output, rows, delimiter=',', fmt='%.17g')
            print(f'Z={zp}: delta spectrum and stopping for {len(targets)} targets', flush=True)
        low_atima = Path(__file__).resolve().parent.parent/'data/atima_low.csv.xz'
        with lzma.open(low_atima, 'rb') as source, (args.output/'atima_low.csv').open('wb') as output:
            shutil.copyfileobj(source, output)
    finally:
        spectra.archive.close()


ELASTIC_HEADER = struct.Struct('<8sIIiiiiQQdd')
ELASTIC_KEYS = ('energies_MeVu', 'offsets', 'sigma_elastic_mb', 'As',
                'p_cm_MeV', 'rutherford_C_mb', 'x', 'R', 'CDF',
                'survival', 'cell_mb')


def validate_elastic(arrays):
    energy, offsets, x = (arrays[k] for k in ('energies_MeVu', 'offsets', 'x'))
    for key, values in arrays.items():
        size = (len(energy) + 1 if key == 'offsets' else
                len(energy) if key in ELASTIC_KEYS[:6] else len(x))
        if values.shape != (size,) or not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError(f'Invalid elastic array: {key}')
    if (len(energy) < 2 or np.any(energy <= 0) or np.any(np.diff(energy) <= 0)
            or offsets[0] != 0 or offsets[-1] != len(x)
            or np.any(np.diff(offsets) < 2)):
        raise ValueError('Invalid elastic energy grid or offsets')
    if any(np.any(arrays[key] <= 0) for key in ELASTIC_KEYS[2:6]):
        raise ValueError('Elastic cross sections, screening and momentum must be positive')
    for lo, hi in zip(offsets[:-1], offsets[1:]):
        xx, cdf, tail, cell = (arrays[k][lo:hi] for k in ('x', 'CDF', 'survival', 'cell_mb'))
        if (xx[0] != 0 or xx[-1] != 1 or np.any(np.diff(xx) <= 0)
                or cdf[0] != 0 or cdf[-1] != 1 or np.any(np.diff(cdf) < 0)
                or tail[0] != 1 or tail[-1] != 0 or np.any(np.diff(tail) > 0)
                or cell[-1] != 0):
            raise ValueError('Invalid elastic angular grid or cumulative distribution')


def write_elastic(source, target):
    metadata = json.loads(source.with_suffix('.json').read_text(encoding='utf-8-sig'))
    with np.load(source, allow_pickle=False) as data:
        arrays = {key: np.asarray(data[key], dtype='<i8' if key == 'offsets' else '<f8')
                  for key in ELASTIC_KEYS}
        ids = tuple(int(data[key]) for key in ('projectile_Z', 'projectile_A',
                                               'target_Z', 'target_A'))
        masses = (float(data['target_mass_MeV']), float(metadata['projectile']['mass_MeV']))
    if ids[:2] != (metadata['projectile']['Z'], metadata['projectile']['A']):
        raise ValueError(f'Projectile differs between NPZ and JSON: {source.name}')
    if not np.isfinite(masses).all() or min(masses) <= 0:
        raise ValueError(f'Invalid nuclear mass: {source.name}')
    validate_elastic(arrays)
    header = ELASTIC_HEADER.pack(b'G4ELT001', 1, 256, *ids,
                                 len(arrays['energies_MeVu']), len(arrays['x']), *masses)
    temporary = target.with_suffix('.partial')
    with temporary.open('wb') as stream:
        stream.write(header)
        stream.write(bytes(256 - len(header)))
        for values in arrays.values():
            values.tofile(stream)
    temporary.replace(target)
    return ids


def convert_elastic(args):
    source = args.source / 'packed' if (args.source / 'packed').is_dir() else args.source
    selected = []
    for path in sorted(source.glob('*.npz')):
        with np.load(path, allow_pickle=False) as data:
            zp, zt = int(data['projectile_Z']), int(data['target_Z'])
        if ((args.projectiles is None or zp in args.projectiles)
                and (args.targets is None or zt in args.targets)):
            selected.append((path, zp, zt))
    if not selected:
        raise ValueError('No elastic tables match the requested elements')
    if len({(zp, zt) for _, zp, zt in selected}) != len(selected):
        raise ValueError('Only one elastic table per projectile-Z/target-Z pair is supported')
    missing_p = set(args.projectiles or ()) - {zp for _, zp, _ in selected}
    missing_t = set(args.targets or ()) - {zt for _, _, zt in selected}
    if missing_p or missing_t:
        raise ValueError(f'Elastic inputs missing for projectile Z={sorted(missing_p)}, '
                         f'target Z={sorted(missing_t)}')
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for path, _, _ in selected:
        target = args.output / (path.stem + '.eltb')
        ids = write_elastic(path, target)
        rows.append('\t'.join(map(str, (target.name, *ids))))
        print(target.name, flush=True)
    index = args.output / 'tables.index.tsv'
    index.write_text(f'# G4ELASTIC_TABLE_INDEX 1 {len(rows)}\n' + '\n'.join(rows) + '\n',
                     encoding='ascii', newline='\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='mode', required=True)
    dirac = commands.add_parser('dirac', help='Dirac spectra, GLOBAL rates and ATIMA range splines')
    dirac.add_argument('--spectra', type=Path, required=True, help='Dirac NPZ archive (ZIP)')
    dirac.add_argument('--global', dest='global_csv', type=Path, required=True, help='GLOBAL transition CSV')
    dirac.add_argument('--atima', type=Path, required=True, help='ATIMA range-spline TAR archive')
    dirac.set_defaults(convert=convert_dirac)
    elastic = commands.add_parser('elastic', help='Elastic NPZ angular tables with JSON mass metadata')
    elastic.add_argument('--source', type=Path, required=True, help='Directory containing NPZ and JSON files')
    elastic.set_defaults(convert=convert_elastic)
    for command in (dirac, elastic):
        command.add_argument('--output', type=Path, required=True, help='Empty output directory')
        command.add_argument('--projectiles', type=int, nargs='+', help='Projectile atomic numbers (default: all)')
        command.add_argument('--targets', type=int, nargs='+', help='Target atomic numbers (default: all)')
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Choose an empty output directory to avoid mixing different table selections')
    args.output.mkdir(parents=True, exist_ok=True)
    args.convert(args)


if __name__ == '__main__':
    main()
