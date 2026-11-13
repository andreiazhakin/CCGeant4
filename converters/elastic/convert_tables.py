"""Convert packed silicon elastic tables to the G4ELT001 binary format."""
from __future__ import annotations
import argparse, json, math, struct, time
from pathlib import Path
import numpy as np

MAGIC = b'G4ELT001'
HEADER = struct.Struct('<8sIIiiiiQQdd')
HEADER_BYTES = 256
KEYS = ['energies_MeVu', 'offsets', 'sigma_elastic_mb', 'As', 'p_cm_MeV',
        'rutherford_C_mb', 'x', 'R', 'CDF', 'survival', 'cell_mb']

def integral_t(t):
    return np.where(abs(t) < 1e-3,
        t*t*(.5+t*(-2/3+t*(.75+t*(-.8+t*5/6)))),
        np.log1p(t)-t/(1+t))

def validate(a):
    e, off = a['energies_MeVu'], a['offsets']
    for k, v in a.items():
        size = len(e)+1 if k == 'offsets' else len(e) if k in KEYS[:6] else len(a['x'])
        assert v.shape == (size,), k
        assert np.isfinite(v).all() and np.all(v >= 0), k
    assert len(e) == 207 and e[0] == 10 and e[-1] == 5000
    assert np.all(np.diff(e) > 0) and off[0] == 0
    assert np.all(np.diff(off) >= 2) and off[-1] == len(a['x'])
    assert np.all(a['p_cm_MeV'] > 0)
    normalization = measure = complement = 0.
    for i, (lo, hi) in enumerate(zip(off[:-1], off[1:])):
        x,r,cdf,surv,cell = (a[k][lo:hi] for k in ['x','R','CDF','survival','cell_mb'])
        assert x[0] == 0 and x[-1] == 1 and np.all(np.diff(x) > 0)
        assert cdf[0] == 0 and cdf[-1] == 1 and np.all(np.diff(cdf) >= 0)
        assert surv[0] == 1 and surv[-1] == 0 and np.all(np.diff(surv) <= 0)
        assert cell[-1] == 0 and np.all(cell[:-1] > 0)
        sigma, A, C = (a[k][i] for k in ['sigma_elastic_mb','As','rutherford_C_mb'])
        assert sigma > 0 and A > 0 and C > 0
        normalization = max(normalization, abs(math.fsum(cell)/sigma-1))
        complement = max(complement, float(np.max(abs(cdf+surv-1))))
        d = np.diff(x); ax = A+x[:-1]
        got = C*(r[:-1]*d/(ax*(ax+d))+np.diff(r)/d*integral_t(d/ax))
        measure = max(measure, float(np.max(abs(got-cell[:-1])/cell[:-1])))
    assert normalization < 1e-8 and complement < 1e-8 and measure < 1e-8
    return dict(max_sigma_relative_error=normalization,
                max_cell_measure_relative_error=measure,
                max_CDF_survival_complement_error=complement)

def convert_one(source, target):
    meta_path = source.with_suffix('.json')
    meta = json.loads(meta_path.read_text(encoding='utf-8-sig'))
    with np.load(source, allow_pickle=False) as z:
        a = {k:np.asarray(z[k], dtype='<i8' if k == 'offsets' else '<f8') for k in KEYS}
        ids = [int(z[k]) for k in ['projectile_Z','projectile_A','target_Z','target_A']]
        target_mass = float(z['target_mass_MeV'])
    projectile_mass = float(meta['projectile']['mass_MeV'])
    assert ids[:2] == [meta['projectile']['Z'],meta['projectile']['A']]
    assert ids[2:] == [14,28]
    assert all(math.isfinite(m) and m > 0 for m in (target_mass, projectile_mass))
    checks = validate(a)
    n, npoint = len(a['energies_MeVu']), len(a['x'])
    header = HEADER.pack(MAGIC,1,HEADER_BYTES,*ids,n,npoint,target_mass,
                         projectile_mass)
    tmp = target.with_suffix('.partial')
    with tmp.open('wb') as f:
        f.write(header); f.write(bytes(HEADER_BYTES-len(header)))
        for k in KEYS: a[k].tofile(f)
    tmp.replace(target)
    cursor = HEADER_BYTES
    for k, original in a.items():
        copied = np.memmap(target,mode='r',offset=cursor,dtype=original.dtype,shape=original.shape)
        assert np.array_equal(copied.view('u8'), original.view('u8')), (source,k)
        cursor += original.nbytes
        del copied
    assert cursor == target.stat().st_size
    return dict(file=target.name, projectile=meta['projectile'],
                target=dict(Z=14,A=28,mass_MeV=target_mass),
                source_npz=str(source.resolve()),
                bytes=cursor,energy_nodes=n,points=npoint,
                all_arrays_bitwise_equal=True,checks=checks)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True); start=time.time(); rows=[]
    files=sorted((args.source/'packed').glob('*_Si28.npz'))
    assert len(files)==57
    for i, source in enumerate(files):
        row=convert_one(source,args.output/(source.stem+'.eltb')); rows.append(row)
        status=dict(state='CONVERTING',completed=i+1,total=57,elapsed_s=time.time()-start)
        (args.output/'conversion_status.json').write_text(json.dumps(status,indent=2))
        print(f'{i+1}/57 {source.stem}: verified',flush=True)
    with (args.output/'tables.index.tsv').open('w',encoding='ascii',newline='\n') as f:
        f.write('# G4ELASTIC_TABLE_INDEX 1 57\n')
        for r in rows:
            pr=r['projectile']; tr=r['target']
            f.write('\t'.join(map(str,[r['file'],pr['Z'],pr['A'],tr['Z'],tr['A']]))+'\n')
    report=dict(format='G4ELT001',version=1,header_bytes=256,byte_order='little',
        units=dict(energy='MeV/nucleon',cross_section='millibarn',mass='MeV/c^2',
                   momentum='MeV/c',x='sin^2(theta_CM/2)'),
        arrays=[dict(name=k,dtype='int64' if k=='offsets' else 'float64') for k in KEYS],
        interpolation='positive linear mixture of unnormalised d_sigma/d_x in kinetic E/A',
        sampling='piecewise linear unscreened R(x), measure C R(x)/(As+x)^2; stored survival tail',
        total_bytes=sum(r['bytes'] for r in rows),projectile_count=57,
        table_count=sum(r['energy_nodes'] for r in rows),projectiles=rows,
        elapsed_s=time.time()-start,passed=True)
    (args.output/'conversion_manifest.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    (args.output/'conversion_status.json').write_text(json.dumps(dict(state='COMPLETE',passed=True,elapsed_s=time.time()-start),indent=2))

if __name__ == '__main__': main()
