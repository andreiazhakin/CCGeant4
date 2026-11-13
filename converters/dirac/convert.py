#!/usr/bin/env python3
"""Convert original Dirac/GLOBAL/ATIMA tables to per-atom CC data, schema v1."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import time
import zipfile
import numpy as np

from cc_binary import read_pack, write_pack
from cc_inputs import (read_atima, read_global, equilibrium_from_rates,
                       interpolate_global_rates)
from cc_physics import bracket, delta_moments, integrate_spectrum, tmax


def emit(message):
    print(json.dumps(message, ensure_ascii=False, allow_nan=False), flush=True)


def save_json(path, value):
    temporary = path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    temporary.replace(path)


class DiracArchive:
    def __init__(self, path):
        self.path = Path(path)
        self.archive = zipfile.ZipFile(path)
        manifests = [n for n in self.archive.namelist() if n.endswith('/archive_manifest.json')]
        if len(manifests) != 1:
            raise ValueError('Expected one Dirac archive manifest')
        self.prefix = manifests[0][:-len('archive_manifest.json')]
        self.grid = self.arrays('full_sweep/results/grid.npz')
        e, off, t = (self.grid[k] for k in ['E_u_MeV_u','offsets','T_eV'])
        if (e.ndim != 1 or len(e) < 2 or not np.all(np.isfinite(e))
                or np.any(np.diff(e) <= 0) or t.ndim != 1 or off.dtype.kind not in 'iu'
                or e[0] != 30 or e[-1] != 2000 or off.shape != (len(e)+1,)
                or off[0] != 0 or off[-1] != len(t) or np.any(np.diff(off) < 2)):
            raise ValueError('Unsupported Dirac energy/transfer grid')
        for i, energy in enumerate(e):
            ts = t[off[i]:off[i+1]]
            if (ts[0] != 1000 or np.any(np.diff(ts) <= 0)
                    or not np.all(np.isfinite(ts)) or abs(ts[-1]/1e6/tmax(energy)-1)>1e-12):
                raise ValueError('Invalid Dirac transfer grid/Tmax')

    def arrays(self, relative):
        raw = self.archive.read(self.prefix+relative)
        with np.load(io.BytesIO(raw), allow_pickle=False) as file:
            return {k: file[k] for k in file.files}

    def projectile(self, zp):
        nmax = 10 if zp <= 50 else 15
        bare = self.arrays(f'full_sweep/results/bare_Z{zp}.npz')
        d0 = bare['D0_cm2_per_eV']
        schema = 'delta-with-canonical-D0-v1'
        nenergy = len(self.grid['E_u_MeV_u'])
        if (str(bare['table_schema']) != schema or int(bare['Z']) != zp
                or bare['valid_energy'].shape != (nenergy,) or bare['valid_energy'].dtype != np.dtype('bool')
                or not np.all(bare['valid_energy']) or d0.shape != self.grid['T_eV'].shape
                or np.any(~np.isfinite(d0)) or np.any(d0 <= 0)):
            raise ValueError(f'Invalid bare spectrum Z={zp}')
        for key in ['numerical_error_abs','numerical_error_rel']:
            if bare[key].shape != d0.shape or np.any(~np.isfinite(bare[key])) or np.any(bare[key] < 0):
                raise ValueError(f'Invalid bare numerical diagnostics Z={zp}')
        density = []
        metadata = {'Zp': zp, 'N_calculated': [0,nmax], 'N_technical_extension': [nmax+1,27],
                    'maximum_bare_numerical_relative_error': float(bare['numerical_error_rel'].max()),
                    'maximum_screened_ratio_numerical_relative_error': 0.,
                    'signed_delta_min': 0., 'signed_delta_max': 0.}
        for n in range(nmax+1):
            data = self.arrays(f'full_sweep/results/Z{zp}_N{n}.npz')
            delta = data['delta_Z'].astype(float)
            if (str(data['table_schema']) != schema or int(data['Z']) != zp or int(data['N']) != n
                    or data['valid_energy'].shape != (nenergy,) or data['valid_energy'].dtype != np.dtype('bool')
                    or not np.all(data['valid_energy']) or delta.shape != d0.shape
                    or np.any(~np.isfinite(delta)) or np.any(delta >= zp)
                    or (n == 0 and np.any(delta != 0))):
                raise ValueError(f'Invalid screened spectrum Z={zp}, N={n}')
            for key in ['numerical_error_abs','R_num_error_rel']:
                if data[key].shape != delta.shape or np.any(~np.isfinite(data[key])) or np.any(data[key] < 0):
                    raise ValueError(f'Invalid Dirac numerical diagnostics Z={zp} N={n}')
            density.append(d0*(1-delta/zp)**2*1e6)
            metadata['maximum_screened_ratio_numerical_relative_error'] = max(
                metadata['maximum_screened_ratio_numerical_relative_error'],float(data['R_num_error_rel'].max()))
            metadata['signed_delta_min'] = min(metadata['signed_delta_min'],float(delta.min()))
            metadata['signed_delta_max'] = max(metadata['signed_delta_max'],float(delta.max()))
        d = np.array(density)
        if not np.all(np.isfinite(d)) or np.any(d <= 0):
            raise ValueError('Nonpositive reconstructed spectrum')
        arrays = {'energy_MeV_u': self.grid['E_u_MeV_u'], 'offsets': self.grid['offsets'],
                  'T_MeV': self.grid['T_eV']/1e6, 'D_cm2_MeV': d}
        tails = [np.empty_like(d) for _ in range(3)]
        off = arrays['offsets']
        for i in range(len(off)-1):
            sl = slice(off[i],off[i+1])
            for k, tail in enumerate(integrate_spectrum(arrays['T_MeV'][sl],d[:,sl])):
                tails[k][:,sl] = tail
        arrays.update(zip(['tail0_cm2','tail1_MeV_cm2','tail2_MeV2_cm2'],tails))
        return arrays, metadata


def rates_at(g, zi, energy):
    _, lo, hi, w = bracket(g.energies_MeV_u, energy)
    a = g.sigma_cm2[zi][:,lo,:]
    b = g.sigma_cm2[zi][:,hi,:]
    return interpolate_global_rates(a, b, w[None,:,None])


def normalization_target(g, zi, ti, atima, delta, energy):
    zp = int(g.projectile_z[zi]); zt = int(g.target_z[ti])
    energy, lo, hi, w = bracket(g.energies_MeV_u, energy)
    a = g.sigma_cm2[zi,ti,lo,:]; b = g.sigma_cm2[zi,ti,hi,:]
    rates = interpolate_global_rates(a, b, w[:,None])
    p, diagnostic = equilibrium_from_rates(rates,g.transitions,return_diagnostics=True)
    spline = atima[zp,zt]
    stopping = spline.atomic_stopping(energy)
    stopping_left = spline.atomic_stopping(energy,side='left')
    hard = zt*delta_moments(delta,energy,zp)[...,1]
    q2 = (zp-np.arange(28,dtype=float))**2
    q2eq = np.sum(p*q2,axis=-1)
    heq = np.sum(p*hard,axis=-1)
    soft = (stopping-heq)/q2eq
    soft_left = (stopping_left-heq)/q2eq
    if (np.any(~np.isfinite(soft)) or np.any(soft < 0)
            or np.any(~np.isfinite(soft_left)) or np.any(soft_left < 0)):
        raise ValueError(f'Negative/nonfinite soft stopping: Z={zp} Zt={zt}')
    full = soft[...,None]*q2+hard
    residual = float(np.max(np.abs(np.sum(p*full,axis=-1)/stopping-1)))
    if residual > 1e-12 or np.any(full <= 0):
        raise ValueError('ATIMA mean stopping identity failed')
    return {'energy_MeV_u': np.asarray(energy,dtype='<f8'),
            'atima_MeV_cm2': stopping, 'atima_left_MeV_cm2': stopping_left,
            'B_MeV_cm2': soft, 'B_left_MeV_cm2': soft_left, 'q2_eq': q2eq,
            'H_eq_MeV_cm2': heq, 'p_eq': p, 'S_full_MeV_cm2': full}, diagnostic


def prepare_target(g, zi, ti, atima, delta, tolerance):
    zp,zt = int(g.projectile_z[zi]),int(g.target_z[ti])
    jumps = atima[zp,zt].discontinuities_MeV_u
    jumps = jumps[(jumps>30)&(jumps<2000)]
    grid = np.unique(np.r_[g.energies_MeV_u,delta['energy_MeV_u'],jumps])
    history = []
    for iteration in range(13):
        nodes, diagnostics = normalization_target(g,zi,ti,atima,delta,grid)
        fractions = np.array([.25,.5,.75])
        probes = np.exp(np.log(grid[:-1,None])
                        +np.log(grid[1:,None]/grid[:-1,None])*fractions)
        actual, _ = normalization_target(g,zi,ti,atima,delta,probes.ravel())
        interp_b = ((1-fractions[None,:])*nodes['B_MeV_cm2'][:-1,None]
                    +fractions[None,:]*nodes['B_left_MeV_cm2'][1:,None]).ravel()
        error = np.abs((interp_b-actual['B_MeV_cm2'])*actual['q2_eq']/actual['atima_MeV_cm2'])
        per_interval = np.max(error.reshape(len(grid)-1,3),axis=1)
        worst = float(per_interval.max())
        history.append({'iteration':iteration,'energy_nodes':len(grid),'max_relative_mean_error':worst})
        if worst <= 0.9*tolerance:
            return nodes, {'Zt':zt,'refinement':history,'stationary_solver':diagnostics,
                           'ATIMA_discontinuities':len(jumps),
                           'minimum_B_MeV_cm2':float(min(nodes['B_MeV_cm2'].min(),nodes['B_left_MeV_cm2'].min())),
                           'hard_fraction_min':float((nodes['H_eq_MeV_cm2']/nodes['atima_MeV_cm2']).min()),
                           'hard_fraction_max':float((nodes['H_eq_MeV_cm2']/nodes['atima_MeV_cm2']).max())}
        grid = np.unique(np.r_[grid,np.sqrt(grid[:-1]*grid[1:])[per_interval>0.9*tolerance]])
        if len(grid)>10000:
            break
    raise ValueError(f'Energy refinement failed for Z={zp} Zt={zt}: {history[-1]}')


def prepare_stopping(g, zi, atima, delta, tolerance):
    zp = int(g.projectile_z[zi])
    targets, diagnostics, offsets = [], [], [0]
    for ti, zt in enumerate(g.target_z):
        data, report = prepare_target(g,zi,ti,atima,delta,tolerance)
        targets.append(data); diagnostics.append(report)
        offsets.append(offsets[-1]+len(data['energy_MeV_u']))
    result = {'target_Z':g.target_z.astype('<i4'),
              'target_mass_g_mol':np.array([atima[zp,int(z)].target_mass_g_mol for z in g.target_z]),
              'offsets':np.array(offsets,dtype='<i8')}
    result.update({k:np.concatenate([t[k] for t in targets],axis=0) for k in targets[0]})
    summary = {'targets':diagnostics,
               'refinement':[{'energy_nodes':offsets[-1],
                              'max_relative_mean_error':max(d['refinement'][-1]['max_relative_mean_error'] for d in diagnostics)}],
               'nodes_per_target_min':min(np.diff(offsets).tolist()),
               'nodes_per_target_max':max(np.diff(offsets).tolist()),
               'minimum_B_MeV_cm2':min(d['minimum_B_MeV_cm2'] for d in diagnostics),
               'hard_fraction_min':min(d['hard_fraction_min'] for d in diagnostics),
               'hard_fraction_max':max(d['hard_fraction_max'] for d in diagnostics)}
    return result,summary


def ensure_pack(out, relative, kind, zp, arrays, resume):
    path = out/relative
    if path.exists() and resume:
        header, old = read_pack(path)
        from cc_binary import KINDS
        if header != {'kind':KINDS[kind],'Zp':zp} or old.keys() != arrays.keys():
            raise ValueError(f'Cannot resume incompatible output {path}')
        for key in arrays:
            if not np.array_equal(old[key],arrays[key]):
                raise ValueError(f'Cannot resume different output {path}: {key}')
        record = {'path':relative,'bytes':path.stat().st_size,'kind':kind,'Zp':zp,
                  'arrays':{key:{'dtype':a.dtype.kind+str(a.dtype.itemsize),'shape':list(a.shape)} for key,a in old.items()}}
    else:
        partial = path.with_suffix(path.suffix+'.partial')
        if resume and partial.exists():
            if (partial.is_symlink() or not partial.is_file()
                    or not partial.resolve().is_relative_to(out.resolve())):
                raise ValueError(f'Unsafe partial output path: {partial}')
            partial.unlink()
        record = write_pack(path,kind,zp,arrays)
        record['path'] = relative
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dirac',type=Path,required=True)
    parser.add_argument('--global-csv',type=Path,required=True)
    parser.add_argument('--atima-splines',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--z-min',type=int,default=29)
    parser.add_argument('--z-max',type=int,default=92)
    parser.add_argument('--normalization-tolerance',type=float,default=1e-4)
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args()
    if not 29 <= args.z_min <= args.z_max <= 92:
        parser.error('Dirac coverage requires 29 <= Zmin <= Zmax <= 92')
    if not 1e-9 <= args.normalization_tolerance <= 1e-2:
        parser.error('Tolerance must be in [1e-9,1e-2]')
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()) and not args.resume:
        parser.error('Output directory is not empty; choose a new path or explicit --resume')
    out.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    config = {'Z_min':args.z_min,'Z_max':args.z_max,'energy_MeV_u':[30.,2000.],
              'N':[0,27],'normalization_relative_tolerance':args.normalization_tolerance}
    emit({'phase':'read_inputs',**config})
    dirac = DiracArchive(args.dirac)
    g = read_global(args.global_csv,args.z_min,args.z_max,progress=emit)
    atima = read_atima(args.atima_splines,g.projectile_z,g.target_z,progress=emit)
    sources = {'dirac':{'path':args.dirac.name,'bytes':args.dirac.stat().st_size},
               'GLOBAL':{'path':args.global_csv.name,'bytes':g.source_bytes,
                         'source_rows':g.source_rows,'selected_rows':g.selected_rows,'source_runs':g.source_runs},
               'ATIMA':{'path':args.atima_splines.name,'bytes':atima.source_bytes,
                        'archive_members':atima.source_members,'selected_element_members':len(atima)}}
    previous_path = out/'manifest.json'
    if previous_path.exists():
        previous = json.loads(previous_path.read_text(encoding='utf-8'))
        previous_config = dict(previous['configuration'])
        if not any((out/'stopping').glob('*.ccbin')):
            previous_config['normalization_relative_tolerance'] = config['normalization_relative_tolerance']
        if previous_config != config:
            raise ValueError('Resume configuration differs')
    manifest = {'schema':'CC_Dirac_AtomicTables_v1','complete':False,
                'created_utc':datetime.now(timezone.utc).isoformat(), 'configuration':config,'sources':sources,
                'target_Z':g.target_z.tolist(),'files':[],'projectiles':[],
                'energy_convention':'Native numeric MeV/u of the inputs; no exact-isotope T/A remap. Dirac gamma=1+E/931.494.',
                'delta_Tmax':'2*0.51099895*((1+E_MeV_u/931.494)^2-1) MeV, infinite projectile mass',
                'delta_low_transfer_MeV':.001,
                'N_extension':'Nmax=10 for Z<=50, else15; N>Nmax: D_N=D_Nmax*((Z-N)/(Z-Nmax))^2; no truncation of GLOBAL states',
                'interpolation':{'GLOBAL':'log sigma vs log E for positive endpoints; linear sigma vs log E if either endpoint is zero; stationary p solved after rate interpolation',
                                 'delta_within_T':'positive piecewise power law, exact analytic cell moments',
                                 'delta_between_E':'linear mixture in log E after affine mapping u=(T-0.001)/(Tmax(E)-0.001), including Jacobian',
                                 'stopping':'Per-target ragged grids. Linear B vs log E from lower-node right value to upper-node left value; exact nodes use right value. S_N=B*q_N^2+Zt*M1_N; S_full contains right node values, not an independent interpolation law.'}}
    save_json(out/'manifest.json',manifest)
    atima_records = []
    for (zp,zt), s in sorted(atima.items()):
        atima_records.append({'Zp':zp,'Zt':zt,'member':s.source_member,
                              'projectile_reference_A':s.reference_A,'projectile_reference_mass_u':s.mass_ref_u,
                              'target_mass_g_mol':s.target_mass_g_mol,'density_g_cm3':s.density_g_cm3,
                              'spline_degree':s.degree,'header':s.header_lines})
    save_json(out/'atima_source_members.json',atima_records)
    for zi,zp0 in enumerate(g.projectile_z):
        zp = int(zp0)
        t0 = time.monotonic()
        delta, dm = dirac.projectile(zp)
        manifest['files'].append(ensure_pack(out,f'delta/Z{zp:03d}.ccbin','delta',zp,delta,args.resume))
        global_arrays = {'energy_MeV_u':g.energies_MeV_u,
                         'target_Z':g.target_z.astype('<i4'),'target_mass_g_mol':g.target_mass_g_mol,
                         'transitions_N':g.transitions.astype('<i4'),'sigma_cm2':g.sigma_cm2[zi],
                         'calc_energy_MeV_u':g.actual_energies_MeV_u}
        manifest['files'].append(ensure_pack(out,f'global/Z{zp:03d}.ccbin','global',zp,global_arrays,args.resume))
        stop, diagnostic = prepare_stopping(g,zi,atima,delta,args.normalization_tolerance)
        manifest['files'].append(ensure_pack(out,f'stopping/Z{zp:03d}.ccbin','stopping',zp,stop,args.resume))
        dm['stopping_energy_nodes'] = len(stop['energy_MeV_u'])
        dm['normalization'] = diagnostic
        manifest['projectiles'].append(dm)
        manifest['elapsed_seconds'] = time.monotonic()-started
        save_json(out/'manifest.json',manifest)
        emit({'phase':'projectile_complete','Zp':zp,'stopping_energy_nodes':dm['stopping_energy_nodes'],
              'max_offgrid_normalization_error':diagnostic['refinement'][-1]['max_relative_mean_error'],
              'projectile_seconds':round(time.monotonic()-t0,2),
              'total_seconds':round(time.monotonic()-started,2)})
    manifest['total_binary_bytes'] = sum(row['bytes'] for row in manifest['files'])
    manifest['complete'] = True
    manifest['elapsed_seconds'] = time.monotonic()-started
    save_json(out/'manifest.json',manifest)
    emit({'phase':'complete','binary_files':len(manifest['files']),
          'binary_bytes':manifest['total_binary_bytes'],'seconds':manifest['elapsed_seconds']})


if __name__=='__main__':
    main()
