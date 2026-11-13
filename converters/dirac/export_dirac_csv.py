"""Export validated CCPACK01 converter output to portable transport CSV files.

No manifest is needed by the Geant4 CSV reader. This converter intentionally
keeps only the stopping coefficients and spectra required during transport.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
if __package__:
    from .converter_import import load_read_pack
else:
    from converter_import import load_read_pack

DELTA_HEADER = 'Zp,Nmax,E_MeV_u,T_MeV,' + ','.join(f'D{n}_cm2_MeV' for n in range(16))
STOP_HEADER = 'Zp,Zt,E_MeV_u,B_MeV_cm2,B_left_MeV_cm2,ATIMA_MeV_cm2,ATIMA_left_MeV_cm2'


def export(source, output, zs, converter_dir=None):
    read_pack = load_read_pack(converter_dir)
    (output / 'delta').mkdir(parents=True, exist_ok=True)
    (output / 'stopping').mkdir(parents=True, exist_ok=True)
    result = []
    for zp in zs:
        start = time.monotonic()
        _, d = read_pack(source / 'delta' / f'Z{zp:03d}.ccbin')
        _, s = read_pack(source / 'stopping' / f'Z{zp:03d}.ccbin')
        nmax = d['D_cm2_MeV'].shape[0] - 1
        delta_path = output / 'delta' / f'Z{zp:03d}.csv'
        stop_path = output / 'stopping' / f'Z{zp:03d}.csv'
        if delta_path.exists() or stop_path.exists():
            raise FileExistsError(f'Refusing to overwrite Z={zp}')
        with delta_path.with_suffix('.csv.partial').open('w', encoding='ascii', newline='') as f:
            f.write(DELTA_HEADER + '\n')
            for ei, energy in enumerate(d['energy_MeV_u']):
                first, last = map(int, d['offsets'][ei:ei+2])
                rows = np.zeros((last-first, 20))
                rows[:, 0] = zp
                rows[:, 1] = nmax
                rows[:, 2] = energy
                rows[:, 3] = d['T_MeV'][first:last]
                rows[:, 4:5+nmax] = d['D_cm2_MeV'][:, first:last].T
                np.savetxt(f, rows, delimiter=',', fmt='%.17g')
        with stop_path.with_suffix('.csv.partial').open('w', encoding='ascii', newline='') as f:
            f.write(STOP_HEADER + '\n')
            for ti, zt in enumerate(s['target_Z']):
                first, last = map(int, s['offsets'][ti:ti+2])
                rows = np.empty((last-first, 7))
                rows[:, 0] = zp
                rows[:, 1] = zt
                for column, name in enumerate(('energy_MeV_u', 'B_MeV_cm2', 'B_left_MeV_cm2',
                                               'atima_MeV_cm2', 'atima_left_MeV_cm2'), 2):
                    rows[:, column] = s[name][first:last]
                np.savetxt(f, rows, delimiter=',', fmt='%.17g')
        delta_path.with_suffix('.csv.partial').replace(delta_path)
        stop_path.with_suffix('.csv.partial').replace(stop_path)
        row = {'Z': zp, 'delta_bytes': delta_path.stat().st_size,
               'stopping_bytes': stop_path.stat().st_size,
               'seconds': time.monotonic()-start}
        result.append(row)
        print(json.dumps(row), flush=True)
    (output / 'csv_export_report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--z', type=int, nargs='+', default=list(range(29, 93)))
    parser.add_argument('--converter-dir', type=Path,
                        help='Source directory containing cc_binary.py; independent of --source data')
    args = parser.parse_args()
    export(args.source, args.output, args.z, args.converter_dir)
