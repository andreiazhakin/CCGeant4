"""Structural and numerical verification of CC Dirac table packs.

Parse the binary format independently and check moments by Gaussian quadrature.
"""
from __future__ import annotations
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
import csv
import io
import json
import math
from pathlib import Path
import struct
import tempfile
import time
import zipfile
import numpy as np

HEADER = struct.Struct('<8sIIIIQQ216s')
DESCRIPTOR = struct.Struct('<64s4sI6Q8s')
DTYPES = {'f8': np.dtype('<f8'), 'i8': np.dtype('<i8'), 'i4': np.dtype('<i4')}
T0 = .001
NA = 6.02214076e23


def check(condition, message):
    if not bool(condition):
        raise ValueError(message)


def parse_pack(path, expected_kind=None, expected_z=None):
    """Independent bounded reader of the published packed-array schema."""
    path = Path(path)
    size = path.stat().st_size
    with path.open('rb') as f:
        raw = f.read(HEADER.size)
        check(len(raw) == 256, f'{path}: truncated header')
        magic, version, kind, zp, count, header_bytes, file_bytes, reserved = HEADER.unpack(raw)
        check(magic == b'CCPACK01' and version == 1, f'{path}: unsupported schema')
        check(kind in (1, 2, 3) and 1 <= count <= 64, f'{path}: invalid kind/count')
        check(expected_kind is None or expected_kind == kind, f'{path}: incorrect kind')
        check(expected_z is None or expected_z == zp, f'{path}: incorrect projectile Z')
        check(file_bytes == size and header_bytes >= 256 + 128*count and header_bytes % 64 == 0,
              f'{path}: invalid header/file byte lengths')
        check(reserved == bytes(216), f'{path}: nonzero reserved header')
        arrays, spans = {}, []
        for _ in range(count):
            raw = f.read(128)
            check(len(raw) == 128, f'{path}: truncated descriptor')
            name_raw, dtype_raw, ndim, *rest = DESCRIPTOR.unpack(raw)
            shape_all = rest[:4]
            offset, nbytes, pad = rest[4:]
            name = name_raw.split(b'\0', 1)[0].decode('ascii')
            tag = dtype_raw.split(b'\0', 1)[0].decode('ascii')
            check(name and name not in arrays and tag in DTYPES and 1 <= ndim <= 4,
                  f'{path}: invalid descriptor')
            check(pad == bytes(8) and all(x == 0 for x in shape_all[ndim:]),
                  f'{path}: nonzero descriptor reserved fields')
            shape = tuple(shape_all[:ndim])
            check(all(x > 0 for x in shape), f'{path}: zero array dimension')
            check(nbytes == math.prod(shape)*DTYPES[tag].itemsize,
                  f'{path}: shape/byte-size mismatch for {name}')
            check(offset >= header_bytes and offset % 64 == 0 and offset+nbytes <= size,
                  f'{path}: invalid array span for {name}')
            spans.append((offset, offset+nbytes))
            arrays[name] = np.memmap(path, mode='r', dtype=DTYPES[tag], offset=offset, shape=shape)
        spans.sort()
        check(all(a[1] <= b[0] for a, b in zip(spans, spans[1:])), f'{path}: overlapping arrays')
    return arrays


def tmax(energy):
    return 2*.51099895*((1+np.asarray(energy)/931.494)**2-1)


def moments_at(delta, energies, zp):
    """Independent affine-support mixture of complete delta moments."""
    query = np.atleast_1d(np.asarray(energies, dtype=float))
    e = delta['energy_MeV_u']
    check(np.all((query >= e[0]) & (query <= e[-1])), 'delta query outside coverage')
    hi = np.clip(np.searchsorted(e, query, side='right'), 1, len(e)-1)
    lo = hi-1
    w = np.log(query/e[lo])/np.log(e[hi]/e[lo])
    idx = np.asarray(delta['offsets'][:-1])
    native = np.stack([delta[k][:, idx].T for k in
                       ['tail0_cm2', 'tail1_MeV_cm2', 'tail2_MeV2_cm2']], axis=-1)
    a, b = native[lo], native[hi]
    l0 = (tmax(e[lo])-T0)[:, None]
    l1 = (tmax(e[hi])-T0)[:, None]
    length = (tmax(query)-T0)[:, None]
    w = w[:, None]
    m0 = (1-w)*a[..., 0] + w*b[..., 0]
    u1 = (1-w)*(a[..., 1]-T0*a[..., 0])/l0 + w*(b[..., 1]-T0*b[..., 0])/l1
    u2 = ((1-w)*(a[..., 2]-2*T0*a[..., 1]+T0*T0*a[..., 0])/l0**2
          + w*(b[..., 2]-2*T0*b[..., 1]+T0*T0*b[..., 0])/l1**2)
    result = np.stack([m0, T0*m0+length*u1,
                       T0*T0*m0+2*T0*length*u1+length**2*u2], axis=-1)
    out = np.empty((len(query), 28, 3))
    out[:, :result.shape[1]] = result
    nmax = result.shape[1]-1
    for n in range(nmax+1, 28):
        out[:, n] = result[:, nmax]*((zp-n)/(zp-nmax))**2
    return out


def global_probabilities(global_data, query):
    """Build the rate generator directly and solve p K = 0."""
    query = np.atleast_1d(query)
    e = global_data['energy_MeV_u']
    hi = np.clip(np.searchsorted(e, query, side='right'), 1, len(e)-1)
    lo = hi-1
    w = np.log(query/e[lo])/np.log(e[hi]/e[lo])
    s = global_data['sigma_cm2']
    lower, upper = s[:, lo, :], s[:, hi, :]
    weight = np.broadcast_to(w[None, :, None], lower.shape)
    positive = (lower > 0) & (upper > 0)
    rates = lower+weight*(upper-lower)
    rates[positive] = np.exp((1-weight[positive])*np.log(lower[positive])
                             + weight[positive]*np.log(upper[positive]))
    transitions = global_data['transitions_N']
    n = 28
    matrix = np.zeros((*rates.shape[:-1], n, n))
    matrix[..., transitions[:, 0], transitions[:, 1]] = rates
    idx = np.arange(n)
    matrix[..., idx, idx] = -matrix.sum(axis=-1)
    scale = np.max(np.abs(matrix), axis=(-2, -1))
    check(np.all(scale > 0), 'zero GLOBAL transition generator')
    generator = matrix/scale[..., None, None]
    a = np.swapaxes(generator, -2, -1).copy()
    a[..., -1, :] = 1.
    rhs = np.zeros((*rates.shape[:-1], n)); rhs[..., -1] = 1.
    p = np.linalg.solve(a, rhs[..., None])[..., 0]
    check(p.min() > -1e-12, 'negative stationary probability')
    p = np.maximum(p, 0); p /= p.sum(axis=-1, keepdims=True)
    residual = np.max(np.abs(np.einsum('...i,...ij->...j', p, generator)))
    check(residual < 1e-12, 'stationary distribution residual too large')
    return p


def validate_delta(delta, zp, source_zip=None):
    expected = {'energy_MeV_u', 'offsets', 'T_MeV', 'D_cm2_MeV',
                'tail0_cm2', 'tail1_MeV_cm2', 'tail2_MeV2_cm2'}
    check(set(delta) == expected, f'Z{zp}: delta fields differ from schema')
    e, offsets, transfer, density = (delta[k] for k in
                                   ['energy_MeV_u', 'offsets', 'T_MeV', 'D_cm2_MeV'])
    ncalc = 11 if zp <= 50 else 16
    check(len(e) == 39 and e[0] == 30 and e[-1] == 2000 and np.all(np.diff(e)>0), 'bad delta energy grid')
    check(len(offsets) == len(e)+1 and offsets[0] == 0 and offsets[-1] == len(transfer), 'bad offsets')
    check(density.shape == (ncalc, len(transfer)) and np.all(np.isfinite(density)) and np.all(density>0), 'bad D')
    gx, gw = np.polynomial.legendre.leggauss(8)
    fractions, weights = (gx+1)/2, gw/2
    max_integral_error = 0.
    for ei in range(len(e)):
        sl = slice(offsets[ei], offsets[ei+1])
        t, d = transfer[sl], density[:, sl]
        check(t[0] == T0 and np.all(np.diff(t)>0), 'bad transfer grid')
        check(abs(t[-1]/tmax(e[ei])-1)<1e-13, 'bad Tmax')
        logratio = np.log(t[1:]/t[:-1])
        logdensityratio = np.log(d[:, 1:]/d[:, :-1])
        for order, key in enumerate(['tail0_cm2', 'tail1_MeV_cm2', 'tail2_MeV2_cm2']):
            actual = delta[key][:, sl]
            check(actual.shape == d.shape and np.all(np.isfinite(actual)) and np.all(actual>=0), 'bad cumulative moment')
            check(np.all(actual[:, -1] == 0) and np.all(np.diff(actual,axis=1)<=0), 'nonmonotone cumulative moment')
            exponent = logdensityratio+(order+1)*logratio[None, :]
            numerical_cell = d[:, :-1]*t[:-1]**(order+1)*logratio*np.sum(
                weights*np.exp(exponent[..., None]*fractions),axis=-1)
            reference = np.concatenate([np.cumsum(numerical_cell[:, ::-1],axis=1)[:, ::-1],
                                        np.zeros((len(d),1))],axis=1)
            error = float(np.max(np.abs(actual-reference)/reference[:, :1]))
            max_integral_error = max(max_integral_error,error)
    check(max_integral_error < 2e-11, f'Z{zp}: quadrature and stored moments disagree {max_integral_error}')
    max_source_error = None
    if source_zip is not None:
        prefix = 'Dirac_tables_all_Z29-92/full_sweep/results/'
        with np.load(io.BytesIO(source_zip.read(prefix+'grid.npz')),allow_pickle=False) as grid:
            check(np.array_equal(e,grid['E_u_MeV_u']) and np.array_equal(offsets,grid['offsets']),
                  'converted Dirac grids differ from source')
            check(np.allclose(transfer,grid['T_eV']*1e-6,rtol=3e-15,atol=0),
                  'converted transfer energies differ from source')
        with np.load(io.BytesIO(source_zip.read(prefix+f'bare_Z{zp}.npz')),allow_pickle=False) as bare:
            d0 = bare['D0_cm2_per_eV']
            check(np.all(bare['valid_energy']), 'invalid source bare energy')
        max_source_error = 0.
        for n in range(ncalc):
            with np.load(io.BytesIO(source_zip.read(prefix+f'Z{zp}_N{n}.npz')),allow_pickle=False) as pair:
                check(np.all(pair['valid_energy']), 'invalid source pair energy mask')
                reference = d0*(1-pair['delta_Z'].astype(float)/zp)**2*1e6
            error = float(np.max(np.abs(density[n]/reference-1)))
            max_source_error = max(max_source_error,error)
        check(max_source_error < 3e-14, f'Z{zp}: D differs from source reconstruction')
    return {'max_quadrature_relative_error':max_integral_error,
            'max_source_spectrum_relative_error':max_source_error}


def validate_global(g, zp):
    e, zt, tr, rates = (g[k] for k in ['energy_MeV_u','target_Z','transitions_N','sigma_cm2'])
    check(e[0] == 30 and e[-1] == 2000 and len(e)==39 and np.all(np.diff(e)>0), 'bad GLOBAL energy grid')
    expected = {(i,j) for i in range(28) for j in range(28) if abs(i-j) in (1,2)}
    check(set(map(tuple,tr)) == expected and tr.shape == (106,2), 'bad GLOBAL transitions')
    check(np.all(np.diff(zt)>0) and np.all((zt>=1)&(zt<=99)), 'bad target grid')
    check(rates.shape == (len(zt),len(e),len(tr)) and np.all(np.isfinite(rates)) and np.all(rates>=0), 'bad GLOBAL rates')
    check(np.all(np.isfinite(g['target_mass_g_mol'])) and np.all(g['target_mass_g_mol']>0), 'bad target masses')
    calc = g['calc_energy_MeV_u']
    check(calc.shape == e.shape and abs(calc[0]-30.1)<1e-12 and np.allclose(calc[1:],e[1:],rtol=0,atol=1e-12),
          'actual GLOBAL calculation energies lost')


def validate_cut_integrals(data,zp):
    """Numerically integrate spectra above cuts, including mapped supports."""
    from cc_physics import delta_tail
    gx,gw=np.polynomial.legendre.leggauss(8)
    f,wg=(gx+1)/2,gw/2
    grid=data['energy_MeV_u']
    queries=[float(grid[0]),float(np.sqrt(grid[0]*grid[1])),
             float(np.sqrt(grid[24]*grid[25])),float(grid[-1])]
    nmax=data['D_cm2_MeV'].shape[0]-1
    worst=0.;cases=0
    for energy in queries:
        hi=min(max(int(np.searchsorted(grid,energy,side='right')),1),len(grid)-1)
        lo=hi-1; ew=float(np.log(energy/grid[lo])/np.log(grid[hi]/grid[lo]))
        for n in [0,nmax,27]:
            nn=min(n,nmax);scale=((zp-n)/(zp-nn))**2
            for cut_fraction in [.1,.5,.9]:
                cut=T0*(float(tmax(energy))/T0)**cut_fraction
                reference=np.zeros(3)
                for index,weight in [(lo,1-ew),(hi,ew)]:
                    sl=slice(data['offsets'][index],data['offsets'][index+1])
                    source_t=data['T_MeV'][sl];source_d=data['D_cm2_MeV'][nn,sl]
                    a=float((tmax(energy)-T0)/(tmax(grid[index])-T0));b=T0*(1-a)
                    lower=T0+(cut-T0)/a
                    inside=source_t>lower
                    first=np.exp(np.interp(np.log(lower),np.log(source_t),np.log(source_d)))
                    t=np.r_[lower,source_t[inside]];d=np.r_[first,source_d[inside]]
                    log_t=np.log(t[1:]/t[:-1]);log_d=np.log(d[1:]/d[:-1])
                    node_t=t[:-1,None]*np.exp(log_t[:,None]*f)
                    values=np.exp((log_d+log_t)[:,None]*f)
                    common=d[:-1]*t[:-1]*log_t
                    for order in range(3):
                        reference[order]+=weight*np.sum(common*np.sum(wg*values*(a*node_t+b)**order,axis=1))
                reference*=scale
                actual=delta_tail(data,energy,n,cut,zp)
                check(np.all(actual>0),'nonpositive partial delta moment')
                worst=max(worst,float(np.max(np.abs(actual/reference-1))));cases+=1
    check(worst<3e-11,f'partial delta moments disagree with independent quadrature: {worst}')
    return {'partial_cut_cases':cases,'max_partial_cut_quadrature_relative_error':worst}


def validate_stopping(s, d, g, zp, atima=None, tolerance=1e-4):
    zt, offsets = s['target_Z'],s['offsets']
    points=len(s['energy_MeV_u'])
    check(np.array_equal(zt,g['target_Z']), 'target grids disagree')
    check(len(offsets)==len(zt)+1 and offsets[0]==0 and offsets[-1]==points
          and np.all(np.diff(offsets)>1),'bad stopping target offsets')
    for key,array in s.items():
        check(np.all(np.isfinite(array)),f'{key} contains nonfinite values')
    for key in ['B_MeV_cm2','B_left_MeV_cm2','atima_MeV_cm2','atima_left_MeV_cm2','q2_eq','H_eq_MeV_cm2']:
        check(s[key].shape==(points,),f'bad scalar array shape: {key}')
    check(s['p_eq'].shape==s['S_full_MeV_cm2'].shape==(points,28),'bad stopping state arrays')
    q2=(zp-np.arange(28))**2
    full_error=eq_error=left_eq_error=p_error=0.
    source_error=left_source_error=offgrid_error=0. if atima is not None else None
    offgrid_location=None
    for ti,target in enumerate(zt):
        sl=slice(offsets[ti],offsets[ti+1])
        e,b,bl,at,atl,p,full=(s[k][sl] for k in ['energy_MeV_u','B_MeV_cm2','B_left_MeV_cm2',
                                               'atima_MeV_cm2','atima_left_MeV_cm2','p_eq','S_full_MeV_cm2'])
        check(e[0]==30 and e[-1]==2000 and np.all(np.diff(e)>0),'bad target stopping energy grid')
        check(np.all(b>=0) and np.all(bl>=0) and np.all(at>0) and np.all(atl>0) and np.all(full>0),'nonpositive stopping')
        check(np.min(p)>=0 and np.max(np.abs(p.sum(axis=-1)-1))<1e-12,'bad equilibrium probabilities')
        hard=moments_at(d,e,zp)[...,1]*target
        expected_full=b[:,None]*q2+hard
        full_error=max(full_error,float(np.max(np.abs(full/expected_full-1))))
        eq_error=max(eq_error,float(np.max(np.abs(np.sum(p*full,axis=-1)/at-1))))
        check(np.max(np.abs(s['q2_eq'][sl]/np.sum(p*q2,axis=-1)-1))<2e-12,'q2 mean mismatch')
        check(np.max(np.abs(s['H_eq_MeV_cm2'][sl]/np.sum(p*hard,axis=-1)-1))<2e-12,'hard mean mismatch')
        left_eq_error=max(left_eq_error,float(np.max(np.abs((bl*s['q2_eq'][sl]+s['H_eq_MeV_cm2'][sl])/atl-1))))
        one_target=dict(g)
        one_target['sigma_cm2']=g['sigma_cm2'][ti:ti+1]
        one_target['target_Z']=zt[ti:ti+1]
        one_target['target_mass_g_mol']=g['target_mass_g_mol'][ti:ti+1]
        for start in range(0,len(e),128):
            reference_p=global_probabilities(one_target,e[start:start+128])[0]
            p_error=max(p_error,float(np.max(np.abs(p[start:start+128]-reference_p))))
        if atima is not None:
            spline=atima[zp,int(target)]
            actual_at=spline.atomic_stopping(e)
            source_error=max(source_error,float(np.max(np.abs(at/actual_at-1))))
            actual_left=spline.atomic_stopping(np.nextafter(e,-np.inf))
            left_source_error=max(left_source_error,float(np.max(np.abs(atl/actual_left-1))))
            frac=(3-np.sqrt(5))/2
            probe=np.exp(np.log(e[:-1])+(np.log(e[1:])-np.log(e[:-1]))*frac)
            for start in range(0,len(probe),128):
                en=probe[start:start+128]
                pp=global_probabilities(one_target,en)[0]
                hh=moments_at(d,en,zp)[...,1]*target
                bb=(1-frac)*b[start:start+len(en)]+frac*bl[start+1:start+len(en)+1]
                ref=spline.atomic_stopping(en)
                predicted=np.sum(pp*(bb[:,None]*q2+hh),axis=-1)
                errors=np.abs(predicted/ref-1)
                worst=int(np.argmax(errors))
                if errors[worst]>offgrid_error:
                    offgrid_error=float(errors[worst])
                    offgrid_location={'Zt':int(target),'E_MeV_u':float(en[worst]),
                                      'lower_energy_MeV_u':float(e[start+worst]),
                                      'upper_energy_MeV_u':float(e[start+worst+1])}
    check(full_error<2e-12 and eq_error<2e-12,'node stopping/normalization mismatch')
    check(left_eq_error<2e-12,'left-limit mean stopping mismatch')
    check(p_error<2e-11,'stored equilibrium differs from GLOBAL generator')
    if atima is not None:
        check(source_error<2e-12,'ATIMA source does not match stopping reference')
        check(left_source_error<2e-12,'ATIMA source does not match left limits')
        check(offgrid_error<=tolerance,
              f'offgrid normalization error {offgrid_error} exceeds {tolerance}: {offgrid_location}')
    return {'energy_nodes_total':points,'target_energy_nodes_min':int(np.diff(offsets).min()),
            'target_energy_nodes_max':int(np.diff(offsets).max()),'max_full_node_relative_error':full_error,
            'max_equilibrium_node_relative_error':eq_error,'max_equilibrium_probability_abs_error':p_error,
            'max_left_equilibrium_node_relative_error':left_eq_error,
            'max_ATIMA_source_relative_error':source_error,'max_ATIMA_left_source_relative_error':left_source_error,
            'max_offgrid_equilibrium_relative_error':offgrid_error,'offgrid_max_location':offgrid_location}


def validate_stopping_reader(s,d,zp):
    """Exercise the public Python reader at endpoints and ATIMA jumps."""
    from cc_physics import stopping_at_energy
    selections=[];queries=[30.,2000.]
    for ti,target in enumerate(s['target_Z']):
        lo,hi=map(int,s['offsets'][ti:ti+2])
        b=s['B_MeV_cm2'][lo:hi];bl=s['B_left_MeV_cm2'][lo:hi]
        j=1+int(np.argmax(np.abs(b[1:-1]-bl[1:-1])))
        energy=float(s['energy_MeV_u'][lo+j])
        selections.append((ti,lo+j,energy))
        queries.extend([float(np.nextafter(energy,-np.inf)),energy,
                        float(np.nextafter(energy,np.inf))])
    queries=np.unique(queries)
    lookup={float(e):i for i,e in enumerate(queries)}
    actual=stopping_at_energy(d,s,queries,zp)
    check(actual.shape==(len(s['target_Z']),len(queries),28),'unexpected Python stopping reader shape')
    q2=(zp-np.arange(28))**2
    worst=0.
    for ti,index,energy in selections:
        lo,hi=map(int,s['offsets'][ti:ti+2])
        for node in [lo,hi-1,index]:
            eq=float(s['energy_MeV_u'][node])
            worst=max(worst,float(np.max(np.abs(actual[ti,lookup[eq]]/s['S_full_MeV_cm2'][node]-1))))
        hard=moments_at(d,[energy],zp)[0,:,1]*s['target_Z'][ti]
        left_ref=s['B_left_MeV_cm2'][index]*q2+hard
        left_query=float(np.nextafter(energy,-np.inf))
        worst=max(worst,float(np.max(np.abs(actual[ti,lookup[left_query]]/left_ref-1))))
        right_query=float(np.nextafter(energy,np.inf))
        worst=max(worst,float(np.max(np.abs(actual[ti,lookup[right_query]]/s['S_full_MeV_cm2'][index]-1))))
    check(worst<2e-11,f'Python stopping reader fails at endpoint/jump: {worst}')
    return {'Python_reader_boundary_queries':len(queries),'max_Python_reader_boundary_relative_error':worst}


def negative_tests(global_path,work_dir=None):
    """Malformed data must be rejected before they can enter transport."""
    raw=Path(global_path).read_bytes()
    rejected=[]
    def expect_rejection(name, callback):
        try:
            callback()
        except (ValueError, OSError):
            rejected.append(name)
        else:
            raise ValueError(f'Invalid input was accepted: {name}')
    base=Path(work_dir or Path.cwd()).resolve()
    tmp=Path(tempfile.mkdtemp(prefix='validation_fixture_',dir=base))
    path=tmp/'invalid.ccbin'
    try:
        path.write_bytes(raw[:-1])
        expect_rejection('truncated_binary',lambda:parse_pack(path))
        changed=bytearray(raw)
        struct.pack_into('<Q',changed,256+64+4+4+32,0)
        path.write_bytes(changed)
        expect_rejection('array_overlaps_header',lambda:parse_pack(path))
    finally:
        path.unlink(missing_ok=True)
        tmp.rmdir()
    g=parse_pack(global_path,2)
    malformed=dict(g);tr=np.array(g['transitions_N']);tr[0,0]=-1
    malformed['transitions_N']=tr
    expect_rejection('invalid_GLOBAL_N',lambda:validate_global(malformed,29))
    return rejected


def validate_global_source(path, folder, zs):
    """Compare every selected raw CSV transition with the packed absolute rate.

    A uniform scaling error is invisible to equilibrium probabilities, so this
    independent direct comparison is separate from the stationary checks.
    """
    records={}
    for zp in zs:
        g=parse_pack(folder/'global'/f'Z{zp:03d}.ccbin',2,zp)
        records[zp]=(g,{int(z):i for i,z in enumerate(g['target_Z'])},
                     {float(e):i for i,e in enumerate(g['energy_MeV_u'])},
                     {tuple(map(int,tr)):i for i,tr in enumerate(g['transitions_N'])},
                     np.zeros(g['sigma_cm2'].shape,dtype=bool))
    matched=0
    with Path(path).open('rb',buffering=8*1024*1024) as stream:
        header=next(stream)
        columns=header.decode('utf-8-sig').strip().split(',')
        keys={name:columns.index(name) for name in
              ['Zp','Zt','At','E_MeV_u','Qe_in','Qe_out','sigma_cm2','GLOBAL_E_MeV_u']}
        for line_no,line in enumerate(stream,1):
            fields=line.strip().split(b',')
            zp=int(fields[keys['Zp']])
            if zp not in records:continue
            g,targets,energies,transitions,seen=records[zp]
            zt=int(fields[keys['Zt']]);en=float(fields[keys['E_MeV_u']])
            tr=(int(fields[keys['Qe_in']]),int(fields[keys['Qe_out']]))
            check(zt in targets and en in energies and tr in transitions,'raw GLOBAL coverage differs')
            idx=targets[zt],energies[en],transitions[tr]
            check(not seen[idx],'duplicate raw GLOBAL transition')
            check(float(fields[keys['sigma_cm2']])==float(g['sigma_cm2'][idx]),
                  f'raw GLOBAL rate differs: Z={zp},Zt={zt},E={en},N={tr}')
            check(float(fields[keys['At']])==float(g['target_mass_g_mol'][idx[0]]),'raw GLOBAL target mass differs')
            check(float(fields[keys['GLOBAL_E_MeV_u']])==float(g['calc_energy_MeV_u'][idx[1]]),
                  'raw GLOBAL actual energy differs')
            seen[idx]=True;matched+=1
            if line_no%4_000_000==0:
                print(json.dumps({'phase':'validate_raw_GLOBAL','input_rows':line_no,'matched_rows':matched}),flush=True)
    for zp,(_,_,_,_,seen) in records.items():
        check(np.all(seen),f'raw GLOBAL missing transitions for Z={zp}')
    return {'matched_transition_rows':matched,'all_rates_exactly_equal':True}


def validate_reference_atima(path,folder,zs):
    """Regression against existing independently prepared mass-stopping CSV."""
    selected=set(zs);previous_z=None;count=0;maximum=0.;worst=None;pack=None
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            zp=int(row['Zp'])
            if zp not in selected:continue
            if zp!=previous_z:
                pack=parse_pack(folder/'stopping'/f'Z{zp:03d}.ccbin',3,zp)
                target_index={int(z):i for i,z in enumerate(pack['target_Z'])}
                previous_z=zp
            zt=int(row['Zt']);energy=float(row['E_MeV_u'])
            check(zt in target_index,'reference ATIMA target missing')
            ti=target_index[zt];lo,hi=map(int,pack['offsets'][ti:ti+2])
            grid=pack['energy_MeV_u'][lo:hi]
            ei=int(np.searchsorted(grid,energy))
            check(ei<len(grid) and grid[ei]==energy,'native reference ATIMA energy missing')
            converted=float(pack['atima_MeV_cm2'][lo+ei])*NA/float(pack['target_mass_g_mol'][ti])
            reference=float(row['ATIMA_dEdx_eq_MeV_cm2_g'])
            error=abs(converted/reference-1)
            if error>maximum:
                maximum=error;worst={'Zp':zp,'Zt':zt,'E_MeV_u':energy}
            count+=1
    expected=sum(len(parse_pack(folder/'global'/f'Z{zp:03d}.ccbin',2,zp)['target_Z'])*39 for zp in zs)
    check(count==expected,f'reference ATIMA count mismatch: {count} vs{expected}')
    check(maximum<2e-12,f'reference ATIMA stopping mismatch {maximum}: {worst}')
    return {'path':str(path),'compared_rows':count,
            'max_relative_difference':maximum,'max_location':worst}


def run(args):
    began=time.monotonic(); folder=Path(args.tables)
    manifest=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
    check(args.z is not None or manifest.get('complete') is True,'conversion manifest is incomplete')
    paths=sorted((folder/'delta').glob('Z*.ccbin'))
    if args.z is not None:
        requested=set(args.z)
        paths=[p for p in paths if int(p.stem[1:]) in requested]
        check({int(p.stem[1:]) for p in paths}==requested,'selected projectile files are not ready')
    check(paths,'no delta tables')
    zs=[int(p.stem[1:]) for p in paths]
    check(len(set(zs))==len(zs),'duplicate projectile Z')
    configuration=manifest['configuration']
    if args.z is None:
        check(zs==list(range(configuration['Z_min'],configuration['Z_max']+1)),
              'projectile set differs from configured coverage')
    else:
        check(all(configuration['Z_min']<=z<=configuration['Z_max'] for z in zs),
              'selected projectile is outside configured coverage')
    atima=None; archive=None
    if not args.skip_sources:
        from cc_inputs import read_atima
        first_global=parse_pack(folder/'global'/paths[0].name,2,zs[0])
        print('Reading raw ATIMA reference splines...',flush=True)
        atima=read_atima(args.raw_atima,zs,list(map(int,first_global['target_Z'])))
        archive=zipfile.ZipFile(args.raw_dirac)
    results=[]; total_bytes=0
    try:
        for zp,path in zip(zs,paths):
            d=parse_pack(path,1,zp)
            gp=folder/'global'/path.name;sp=folder/'stopping'/path.name
            g=parse_pack(gp,2,zp);s=parse_pack(sp,3,zp)
            record={'Zp':zp,**validate_delta(d,zp,archive)}
            record.update(validate_cut_integrals(d,zp))
            validate_global(g,zp)
            record.update(validate_stopping(s,d,g,zp,atima,args.tolerance))
            record.update(validate_stopping_reader(s,d,zp))
            results.append(record)
            for pth in [path,gp,sp]:
                total_bytes+=pth.stat().st_size
            record['binary_bytes']={'delta':path.stat().st_size,'global':gp.stat().st_size,
                                    'stopping':sp.stat().st_size,
                                    'total':sum(pth.stat().st_size for pth in [path,gp,sp])}
            print(json.dumps(record),flush=True)
    finally:
        if archive:archive.close()
    global_source=None
    if not args.skip_sources:
        global_source=validate_global_source(args.raw_global,folder,zs)
    reference_atima=None
    if args.reference_atima_csv is not None:
        reference_atima=validate_reference_atima(args.reference_atima_csv,folder,zs)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    rejected=negative_tests(folder/'global'/paths[0].name,args.output.parent)
    summary={'passed':True,'schema':'CCPACK01','projectile_Z':zs,'files_checked':3*len(zs),
             'full_configured_coverage_checked':args.z is None,
             'bytes_checked':total_bytes,'raw_sources_checked':not args.skip_sources,
             'offgrid_requested_tolerance':args.tolerance,'results':results,
             'negative_cases_rejected':rejected,
             'raw_GLOBAL_validation':global_source,
             'reference_ATIMA_csv_validation':reference_atima,
             'elapsed_seconds':time.monotonic()-began}
    Path(args.output).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='results'}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables',type=Path,required=True)
    parser.add_argument('--raw-atima',type=Path)
    parser.add_argument('--raw-dirac',type=Path)
    parser.add_argument('--raw-global',type=Path)
    parser.add_argument('--reference-atima-csv',type=Path,help='Optional legacy ATIMA mass-stopping regression CSV')
    parser.add_argument('--output',type=Path,help='Default: TABLES/validation.json')
    parser.add_argument('--tolerance',type=float,default=1e-4)
    parser.add_argument('--skip-sources',action='store_true')
    parser.add_argument('--z',type=int,nargs='+',help='Validate only selected completed projectiles, also during conversion')
    args=parser.parse_args()
    if not args.skip_sources and any(value is None for value in [args.raw_atima,args.raw_dirac,args.raw_global]):
        parser.error('--raw-atima, --raw-dirac and --raw-global are required unless --skip-sources is set')
    if args.output is None:
        args.output=args.tables/'validation.json'
    run(args)
