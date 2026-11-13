"""Readers for original GLOBAL and ATIMA 1.41 inputs.

ATIMA range splines use the de Boor algorithm. Energies retain the source
coordinate (MeV/u), without isotope or finite-projectile-mass remapping.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
import gzip
import math
from pathlib import Path
import re
import tarfile
from typing import Callable, Iterable

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


def _emit(callback, **message):
    if callback is not None:
        callback(message)


def interpolate_global_rates(lower, upper, log_energy_weight):
    """Interpolate finite nonnegative GLOBAL rates, matching the runtime law.

    Positive endpoint pairs retain log-log interpolation. If either endpoint
    is zero, sigma is linear in log(E); two explicit zeros remain zero. Inputs
    broadcast together. No logarithm is evaluated on the zero-rate branch.
    """
    lo, hi, w = np.broadcast_arrays(np.asarray(lower, dtype=np.float64),
                                    np.asarray(upper, dtype=np.float64),
                                    np.asarray(log_energy_weight, dtype=np.float64))
    if (not np.all(np.isfinite(lo)) or not np.all(np.isfinite(hi))
            or np.any(lo < 0) or np.any(hi < 0)):
        raise ValueError('GLOBAL interpolation rates must be finite and nonnegative')
    if not np.all(np.isfinite(w)) or np.any(w < 0) or np.any(w > 1):
        raise ValueError('GLOBAL interpolation weight must be finite and within [0,1]')
    positive = (lo > 0) & (hi > 0)
    result = np.empty(lo.shape, dtype=np.float64)
    result[positive] = np.exp(np.log(lo[positive])
                              + w[positive]*(np.log(hi[positive])-np.log(lo[positive])))
    linear = ~positive
    result[linear] = lo[linear]+w[linear]*(hi[linear]-lo[linear])
    return result


def equilibrium_from_rates(sigma: np.ndarray, transitions: np.ndarray = TRANSITIONS,
                           *, return_diagnostics: bool = False,
                           batch_size: int = 2048):
    """Stationary probabilities for rates [..., transition], N=0..27.

    Rates may be microscopic cross sections of one element or already summed
    n_t*sigma_t rates of a mixture.  A common positive scale cancels.  Tiny
    negative probabilities from floating-point solves are clipped only after a
    strict tolerance check; the clipped probability mass is reported.
    """
    rates = np.asarray(sigma, dtype=np.float64)
    pairs = np.asarray(transitions)
    if rates.ndim < 1 or pairs.ndim != 2 or pairs.shape[1] != 2 or rates.shape[-1] != len(pairs):
        raise ValueError('Equilibrium: incompatible rate/transition shapes')
    if not np.all(np.isfinite(pairs)) or not np.all(pairs == np.rint(pairs)):
        raise ValueError('Equilibrium: transitions must contain integer states')
    pairs = pairs.astype(np.intp)
    if np.any(pairs < 0) or np.any(pairs >= N_STATES) or np.any(pairs[:, 0] == pairs[:, 1]):
        raise ValueError('Equilibrium: invalid N=0..27 transition')
    if len(set(map(tuple, pairs))) != len(pairs):
        raise ValueError('Equilibrium: duplicate transitions')
    if not np.all(np.isfinite(rates)) or np.any(rates < 0):
        raise ValueError('Equilibrium: rates must be finite and nonnegative')
    if batch_size < 1:
        raise ValueError('Equilibrium: batch_size must be positive')
    flat = rates.reshape(-1, rates.shape[-1])
    out = np.empty((len(flat), N_STATES), dtype=np.float64)
    diagnostics = dict(nodes=len(flat), min_raw_probability=1.0,
                       max_clipped_negative_mass=0.0,
                       max_stationary_residual=0.0, max_sum_residual=0.0)
    diagonal = np.arange(N_STATES)
    for start in range(0, len(flat), batch_size):
        block = flat[start:start+batch_size]
        generator = np.zeros((len(block), N_STATES, N_STATES), dtype=np.float64)
        generator[:, pairs[:, 0], pairs[:, 1]] = block
        generator[:, diagonal, diagonal] = -generator.sum(axis=2)
        scale = np.max(np.abs(generator), axis=(1, 2))
        if np.any(scale <= 0):
            raise ValueError('Equilibrium: zero total transition generator')
        generator /= scale[:, None, None]
        matrix = generator.swapaxes(-1, -2).copy()
        matrix[:, -1, :] = 1.0
        rhs = np.zeros((len(block), N_STATES, 1)); rhs[:, -1, 0] = 1.0
        try:
            p = np.linalg.solve(matrix, rhs)[..., 0]
        except np.linalg.LinAlgError as exc:
            raise ValueError('Equilibrium generator has no unique stationary solution') from exc
        raw_min = float(np.min(p))
        sum_error = float(np.max(np.abs(p.sum(axis=1)-1.0)))
        if not np.all(np.isfinite(p)) or raw_min < -1e-10 or sum_error > 1e-10:
            raise ValueError(f'Invalid stationary solution: min={raw_min}, sum error={sum_error}')
        clipped = float(np.max(-np.minimum(p, 0.0).sum(axis=1)))
        p = np.maximum(p, 0.0); p /= p.sum(axis=1)[:, None]
        residual = float(np.max(np.abs(np.einsum('bi,bij->bj', p, generator))))
        if residual > 1e-10:
            raise ValueError(f'Stationary generator residual exceeds tolerance: {residual}')
        diagnostics['min_raw_probability'] = min(diagnostics['min_raw_probability'], raw_min)
        diagnostics['max_clipped_negative_mass'] = max(diagnostics['max_clipped_negative_mass'], clipped)
        diagnostics['max_stationary_residual'] = max(diagnostics['max_stationary_residual'], residual)
        diagnostics['max_sum_residual'] = max(diagnostics['max_sum_residual'], sum_error)
        out[start:start+len(block)] = p
    result = out.reshape(rates.shape[:-1]+(N_STATES,))
    return (result, diagnostics) if return_diagnostics else result


@dataclass
class GlobalData:
    projectile_z: np.ndarray
    target_z: np.ndarray
    target_mass_g_mol: np.ndarray
    energies_MeV_u: np.ndarray
    actual_energies_MeV_u: np.ndarray
    transitions: np.ndarray
    sigma_cm2: np.ndarray
    source_path: str
    source_bytes: int
    source_rows: int
    selected_rows: int
    source_runs: tuple[str, ...]

    def equilibrium(self, **kwargs):
        return equilibrium_from_rates(self.sigma_cm2, self.transitions, **kwargs)

    def interpolate_rates(self, energy_MeV_u):
        """GLOBAL interpolation on source labelled energy; no extrapolation.

        Positive endpoints use log-log; explicit zeros use linear sigma in logE.

        Scalar energy returns [projectile,target,transition].  A one-dimensional
        energy array returns [projectile,target,energy,transition].
        """
        e = np.asarray(energy_MeV_u, dtype=np.float64)
        if e.ndim > 1 or not np.all(np.isfinite(e)):
            raise ValueError('GLOBAL interpolation expects finite scalar or 1-D energy')
        if np.any(e < self.energies_MeV_u[0]) or np.any(e > self.energies_MeV_u[-1]):
            raise ValueError('GLOBAL energy outside source coverage; extrapolation disabled')
        index = np.clip(np.searchsorted(self.energies_MeV_u, e, side='right')-1,
                        0, len(self.energies_MeV_u)-2)
        e0, e1 = self.energies_MeV_u[index], self.energies_MeV_u[index+1]
        w = np.log(e/e0)/np.log(e1/e0)
        lo, hi = self.sigma_cm2[:, :, index, :], self.sigma_cm2[:, :, index+1, :]
        return interpolate_global_rates(lo, hi, np.expand_dims(w, -1))


def read_global(path: str | Path, zmin: int = 29, zmax: int = 92,
                *, progress: Callable | None = None) -> GlobalData:
    """Read a complete GLOBAL CSV without assuming a fixed number of targets.

    The selected Z interval must be rectangular in target and energy coverage.
    Each group must contain the 106 distinct +/-1,+/-2 transitions within
    N=0..27.  Nominal and actual GLOBAL energies are kept separately.
    """
    path = Path(path)
    if not 1 <= zmin <= zmax:
        raise ValueError('Invalid selected projectile Z range')
    required = ['Zp', 'Ap', 'Zt', 'At', 'E_MeV_u', 'Qe_in', 'Qe_out',
                'dQe', 'kind', 'sigma_cm2', 'GLOBAL_E_MeV_u', 'source_run']
    signature = {tuple(map(int, tr)): i for i, tr in enumerate(TRANSITIONS)}
    curves = {}
    masses = {}
    source_runs = set()
    last = None
    group_values = np.full(len(TRANSITIONS), np.nan)
    group_actual = None
    rows = selected_rows = 0

    def finish_group():
        if last is None:
            return
        if not np.all(np.isfinite(group_values)):
            raise ValueError(f'Incomplete GLOBAL transition group {last}')
        zp, zt, energy = last
        curves.setdefault((zp, zt), []).append((energy, group_actual, group_values.copy()))

    with path.open('rb', buffering=8*1024*1024) as stream:
        raw = next(stream)
        header = raw.decode('utf-8-sig').strip().split(',')
        if header != required:
            raise ValueError(f'Unexpected GLOBAL CSV schema: {header}')
        for rows, raw in enumerate(stream, 1):
            fields = raw.rstrip(b'\r\n').split(b',')
            if len(fields) != len(required):
                raise ValueError(f'GLOBAL row {rows+1}: wrong field count')
            zp = int(fields[0])
            if zp < zmin or zp > zmax:
                continue
            zt = int(fields[2]); at = float(fields[3]); energy = float(fields[4])
            actual = float(fields[10])
            if zt < 1 or not math.isfinite(at) or at <= 0 or not math.isfinite(energy) or energy <= 0 or not math.isfinite(actual) or actual <= 0:
                raise ValueError(f'GLOBAL row {rows+1}: invalid target/energy')
            if zt in masses and masses[zt] != at:
                raise ValueError(f'Inconsistent GLOBAL target mass for Zt={zt}')
            masses[zt] = at
            key = (zp, zt, energy)
            if key != last:
                finish_group()
                group_values.fill(np.nan)
                last = key; group_actual = actual
            elif group_actual != actual:
                raise ValueError(f'Inconsistent actual GLOBAL energy in group {key}')
            n, m, dn = map(int, fields[5:8]); kind = fields[8]
            transition = (n, m)
            if transition not in signature or dn != m-n or kind != (b'capture' if dn > 0 else b'loss'):
                raise ValueError(f'GLOBAL row {rows+1}: invalid transition')
            idx = signature[transition]
            if not math.isnan(group_values[idx]):
                raise ValueError(f'Duplicate GLOBAL transition in {key}: {transition}')
            value = float(fields[9])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'GLOBAL row {rows+1}: negative/nonfinite cross section')
            group_values[idx] = value
            source_runs.add(fields[11].decode('utf-8'))
            selected_rows += 1
            if selected_rows % 3_000_000 == 0:
                _emit(progress, phase='read_global', selected_rows=selected_rows, source_rows=rows)
        finish_group()
    if not curves:
        raise ValueError('No selected GLOBAL data')
    projectile_z = np.arange(zmin, zmax+1, dtype=np.int16)
    target_z = np.array(sorted(masses), dtype=np.int16)
    first = sorted(next(iter(curves.values())), key=lambda x: x[0])
    energies = np.array([x[0] for x in first], dtype=np.float64)
    actual = np.array([x[1] for x in first], dtype=np.float64)
    if len(energies) < 2 or np.any(np.diff(energies) <= 0):
        raise ValueError('Invalid/duplicate GLOBAL energy grid')
    sigma = np.empty((len(projectile_z), len(target_z), len(energies), len(TRANSITIONS)), dtype=np.float64)
    for pi, zp in enumerate(projectile_z):
        for ti, zt in enumerate(target_z):
            key = int(zp), int(zt)
            if key not in curves:
                raise ValueError(f'Missing GLOBAL projectile/target pair {key}')
            curve = sorted(curves.pop(key), key=lambda x: x[0])
            if len(curve) != len(energies) or not np.array_equal([r[0] for r in curve], energies) or not np.array_equal([r[1] for r in curve], actual):
                raise ValueError(f'Inconsistent/incomplete GLOBAL energy grid for {key}')
            sigma[pi, ti] = np.array([r[2] for r in curve])
    if curves or selected_rows != sigma.size:
        raise ValueError('Unexpected GLOBAL coverage/row count')
    _emit(progress, phase='read_global_complete', shape=list(sigma.shape), source_rows=rows)
    return GlobalData(projectile_z, target_z, np.array([masses[int(z)] for z in target_z]),
                      energies, actual, TRANSITIONS.copy(), sigma,
                      str(path.resolve()), path.stat().st_size,
                      rows, selected_rows, tuple(sorted(source_runs)))


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
    projectile_z: int
    target_z: int
    reference_A: int
    mass_ref_u: float
    target_mass_g_mol: float
    density_g_cm3: float
    degree: int
    knots: np.ndarray
    coefficients: np.ndarray
    derivative_coefficients: np.ndarray
    source_member: str
    header_lines: tuple[str, ...]

    @property
    def energy_domain_MeV_u(self):
        return float(self.knots[self.degree]), float(self.knots[len(self.coefficients)])

    @cached_property
    def discontinuities_MeV_u(self):
        """Exact internal source knots at which R' can be discontinuous.

        A degree-p spline has continuity C^(p-m) at a knot of multiplicity m.
        Therefore m >= p permits a discontinuity in the range derivative and
        hence in stopping.  Both limits must be retained by any interpolating
        converter; increasing a continuous interpolation grid cannot reproduce
        the jump.  Domain endpoints are deliberately excluded.
        """
        knots, multiplicities = np.unique(self.knots, return_counts=True)
        low, high = self.energy_domain_MeV_u
        result = knots[(multiplicities >= self.degree) & (knots > low) & (knots < high)]
        result.setflags(write=False)
        return result

    def _side_energy(self, energy_MeV_u, side):
        energy = np.asarray(energy_MeV_u, dtype=np.float64)
        if side not in ('left', 'right'):
            raise ValueError('ATIMA spline side must be left or right')
        if side == 'left':
            energy = np.where(np.isin(energy, self.discontinuities_MeV_u),
                              np.nextafter(energy, -np.inf), energy)
        return energy

    def evaluate_range(self, energy_MeV_u):
        """Reference-projectile CSDA range, mg/cm2, in native ATIMA coordinate."""
        return _bspline(self.knots, self.coefficients, self.degree, energy_MeV_u)

    def range_derivative(self, energy_MeV_u, *, side='right'):
        return _bspline(self.knots[1:-1], self.derivative_coefficients, self.degree-1,
                        self._side_energy(energy_MeV_u, side))

    def mass_stopping(self, energy_MeV_u, *, side='right'):
        """Stopping derived from range, MeV cm2/g, without changing energy axis."""
        derivative = self.range_derivative(energy_MeV_u, side=side)
        if not np.all(np.isfinite(derivative)) or np.any(derivative <= 0):
            raise ValueError(f'Nonpositive ATIMA range derivative: {self.source_member}')
        return 1000.0*self.mass_ref_u/derivative

    def atomic_stopping(self, energy_MeV_u, *, side='right'):
        """Stopping cross section per target atom, MeV cm2."""
        return self.mass_stopping(energy_MeV_u, side=side)*self.target_mass_g_mol/AVOGADRO


class AtimaData(dict):
    """Mapping (Zp,Zt)->AtimaSpline plus archive provenance attributes."""
    source_path: str
    source_bytes: int
    source_members: int


def read_atima(path: str | Path, selected_Z: Iterable[int], selected_targets: Iterable[int],
               *, progress: Callable | None = None) -> AtimaData:
    path = Path(path)
    selected_Z = set(map(int, selected_Z)); selected_targets = set(map(int, selected_targets))
    if not selected_Z or not selected_targets:
        raise ValueError('ATIMA selection must be nonempty')
    if not selected_Z.issubset(range(1, 121)) or not selected_targets.issubset(range(1, 100)):
        raise ValueError('ATIMA selection exceeds source element coverage')
    result = AtimaData()
    members = 0
    with tarfile.open(path, mode='r|') as archive:
        for member in archive:
            members += 1
            if not member.isfile():
                continue
            match = re.fullmatch(r'(?:(E119|E120|[A-Z][a-z]?)(\d+)|(\d{3})X(\d+))-([A-Z][a-z]?)\.dat\.gz', Path(member.name).name)
            if match is None:
                continue
            symbol, atext, numeric_z, numeric_a, target_symbol = match.groups()
            zp = int(numeric_z) if numeric_z is not None else SYMBOL_TO_Z.get(symbol)
            if numeric_a is not None:
                atext = numeric_a
            zt = SYMBOL_TO_Z.get(target_symbol)
            if zp not in selected_Z or zt not in selected_targets:
                continue
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f'Cannot read ATIMA archive member {member.name}')
            raw = stream.read()
            lines = gzip.decompress(raw).decode('ascii').splitlines()
            if len(lines) < 7:
                raise ValueError(f'Truncated ATIMA spline {member.name}')
            projectile = _numbers(lines[0]); components = _numbers(lines[1])
            target = _numbers(lines[2]); material = _numbers(lines[3])
            if len(projectile) != 2 or len(components) != 1 or components[0] != 1 or len(target) != 4 or len(material) < 1:
                raise ValueError(f'Unsupported ATIMA element header {member.name}')
            if projectile[0] != zp or target[0] != zt or target[1] != 1.0:
                raise ValueError(f'ATIMA identity/composition mismatch: {member.name}')
            mass, target_mass, density = float(projectile[1]), float(target[3]), float(material[0])
            if not all(np.isfinite(x) and x > 0 for x in (mass, target_mass, density)):
                raise ValueError(f'Invalid ATIMA masses/density: {member.name}')
            ncoeff, order = map(int, lines[4].split())
            degree = order-1
            knots = _numbers(lines[5]); coefficients = _numbers(lines[6])
            if degree < 1 or ncoeff < order or len(knots) != ncoeff+order or len(coefficients) != ncoeff:
                raise ValueError(f'Invalid ATIMA spline dimensions: {member.name}')
            if not np.all(np.isfinite(knots)) or not np.all(np.isfinite(coefficients)) or np.any(np.diff(knots) < 0) or knots[degree] <= 0 or knots[ncoeff] <= knots[degree]:
                raise ValueError(f'Invalid ATIMA spline grid: {member.name}')
            den = knots[degree+1:degree+ncoeff]-knots[1:ncoeff]
            derivative = np.divide(degree*np.diff(coefficients), den,
                                   out=np.zeros(ncoeff-1), where=den != 0)
            key = int(zp), int(zt)
            if key in result:
                raise ValueError(f'Ambiguous duplicate ATIMA element pair {key}')
            result[key] = AtimaSpline(int(zp), int(zt), int(atext), mass, target_mass, density,
                                     degree, knots, coefficients, derivative, member.name,
                                     tuple(lines[:4]))
            if len(result) % 1000 == 0:
                _emit(progress, phase='read_atima', selected_members=len(result))
    expected = {(zp, zt) for zp in selected_Z for zt in selected_targets}
    if set(result) != expected:
        missing = sorted(expected-set(result))
        raise ValueError(f'Missing ATIMA element pairs: {missing[:20]} ({len(missing)} total)')
    result.source_path = str(path.resolve())
    result.source_bytes = path.stat().st_size
    result.source_members = members
    _emit(progress, phase='read_atima_complete', selected_members=len(result), archive_members=members)
    return result
