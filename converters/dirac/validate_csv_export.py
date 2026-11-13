"""Check CSV schema, coverage and every exported number against CCPACK data.

Without --z, all source projectiles must be exported. --z explicitly selects
a projectile subset; --target-z explicitly selects complete stopping curves.
This checks faithful export, not the physical validity of the source model.
"""
from pathlib import Path
import argparse
import json
import re
import sys
import tempfile
import numpy as np
if __package__:
    from .converter_import import load_read_pack
else:
    from converter_import import load_read_pack

DELTA_HEADER = 'Zp,Nmax,E_MeV_u,T_MeV,' + ','.join(f'D{n}_cm2_MeV' for n in range(16))
STOP_HEADER = 'Zp,Zt,E_MeV_u,B_MeV_cm2,B_left_MeV_cm2,ATIMA_MeV_cm2,ATIMA_left_MeV_cm2'
STOP_FIELDS = ('energy_MeV_u', 'B_MeV_cm2', 'B_left_MeV_cm2',
               'atima_MeV_cm2', 'atima_left_MeV_cm2')


def write_report(report, data):
    report.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
                                         dir=report.parent, prefix=report.name + '.',
                                         suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(data, stream, indent=2, allow_nan=False)
            stream.write('\n')
        temporary.replace(report)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def inventory(root, extension):
    if not root.is_dir():
        raise ValueError(f'Data directory does not exist: {root}')
    result = {}
    for section in ('delta', 'stopping'):
        directory = root / section
        if not directory.is_dir():
            raise ValueError(f'Required directory does not exist: {directory}')
        files = {}
        for path in sorted(directory.glob('*' + extension)):
            match = re.fullmatch(r'Z(\d{3})' + re.escape(extension), path.name)
            if not match or not 29 <= int(match[1]) <= 92:
                raise ValueError(f'Invalid projectile filename: {path}')
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f'Empty or non-file table: {path}')
            files[int(match[1])] = path
        if not files:
            raise ValueError(f'No {extension} projectile tables in {directory}')
        result[section] = files
    if result['delta'].keys() != result['stopping'].keys():
        raise ValueError(f'Unpaired delta/stopping projectiles in {root}: '
                         f'delta={sorted(result["delta"])}, stopping={sorted(result["stopping"])}')
    return result


def selection(values, name, lower, upper):
    if values is None:
        return None
    chosen = list(values)
    if (not chosen or any(not isinstance(z, (int, np.integer)) or not lower <= z <= upper
                          for z in chosen) or len(set(chosen)) != len(chosen)):
        raise ValueError(f'{name} must be a nonempty unique list of integers in {lower}..{upper}')
    return sorted(chosen)


def vector(arrays, name, label):
    if name not in arrays:
        raise ValueError(f'{label}: missing CCPACK array {name}')
    a = arrays[name]
    if a.ndim != 1 or not a.size or not np.all(np.isfinite(a)):
        raise ValueError(f'{label}: invalid CCPACK vector {name}')
    return a


def offsets(arrays, group_count, row_count, label):
    a = vector(arrays, 'offsets', label)
    if (a.dtype.kind not in 'iu' or len(a) != group_count + 1 or a[0] != 0
            or a[-1] != row_count or np.any(np.diff(a) < 2)):
        raise ValueError(f'{label}: invalid CCPACK offsets')
    return np.diff(a)


def read_csv(path, header):
    with path.open('r', encoding='utf-8-sig') as stream:
        if stream.readline().rstrip('\r\n') != header:
            raise ValueError(f'{path}: CSV header/schema mismatch')
        data = np.loadtxt(stream, delimiter=',', comments=None, ndmin=2)
    if (data.shape[0] == 0 or data.shape[1] != len(header.split(','))
            or not np.all(np.isfinite(data))):
        raise ValueError(f'{path}: empty, nonfinite or wrong-width CSV data')
    return data


def equal(actual, expected, label):
    if not np.array_equal(actual, expected):
        raise ValueError(f'{label}: CSV/CCPACK numerical or row-count mismatch '
                         f'(CSV shape {actual.shape}, expected {expected.shape})')


def validate_projectile(z, source, exported, read_pack, targets):
    dm, d = read_pack(source['delta'][z])
    sm, s = read_pack(source['stopping'][z])
    if dm != {'kind': 1, 'Zp': z} or sm != {'kind': 3, 'Zp': z}:
        raise ValueError(f'Z={z}: CCPACK kind/projectile metadata mismatch')
    de = vector(d, 'energy_MeV_u', f'Z={z} delta')
    dt = vector(d, 'T_MeV', f'Z={z} delta')
    counts = offsets(d, len(de), len(dt), f'Z={z} delta')
    spectra = d.get('D_cm2_MeV')
    if (spectra is None or spectra.ndim != 2 or not 1 <= spectra.shape[0] <= 16
            or spectra.shape[1] != len(dt) or not np.all(np.isfinite(spectra))):
        raise ValueError(f'Z={z}: invalid CCPACK delta spectrum shape/values')
    nmax = spectra.shape[0] - 1
    delta = read_csv(exported['delta'][z], DELTA_HEADER)
    equal(delta[:, 0], np.full(len(dt), z), f'Z={z} delta Zp')
    equal(delta[:, 1], np.full(len(dt), nmax), f'Z={z} delta Nmax')
    equal(delta[:, 2], np.repeat(de, counts), f'Z={z} delta energy')
    equal(delta[:, 3], dt, f'Z={z} delta T')
    equal(delta[:, 4:5+nmax], spectra.T, f'Z={z} delta D')
    if np.any(delta[:, 5+nmax:] != 0):
        raise ValueError(f'Z={z}: nonzero delta padding beyond Nmax')

    st = vector(s, 'target_Z', f'Z={z} stopping')
    if (np.any(st != st.astype(int)) or np.any(st < 1) or np.any(st > 96)
            or np.any(np.diff(st) <= 0)):
        raise ValueError(f'Z={z}: invalid CCPACK stopping target list')
    values = {name: vector(s, name, f'Z={z} stopping') for name in STOP_FIELDS}
    row_count = len(values['energy_MeV_u'])
    if any(len(a) != row_count for a in values.values()):
        raise ValueError(f'Z={z}: inconsistent CCPACK stopping array lengths')
    stop_counts = offsets(s, len(st), row_count, f'Z={z} stopping')
    expected_targets = list(map(int, st)) if targets is None else targets
    missing_targets = sorted(set(expected_targets) - set(st))
    if missing_targets:
        raise ValueError(f'Z={z}: requested targets absent in CCPACK: {missing_targets}')
    row_targets = np.repeat(st, stop_counts)
    mask = np.isin(row_targets, expected_targets)
    stopping = read_csv(exported['stopping'][z], STOP_HEADER)
    equal(stopping[:, 0], np.full(np.count_nonzero(mask), z), f'Z={z} stopping Zp')
    equal(stopping[:, 1], row_targets[mask], f'Z={z} stopping targets')
    for ci, name in enumerate(STOP_FIELDS, 2):
        equal(stopping[:, ci], values[name][mask], f'Z={z} stopping {name}')
    return {'Z': z, 'target_Z': expected_targets,
            'delta_rows': len(delta), 'stopping_rows': len(stopping),
            'delta_bytes': exported['delta'][z].stat().st_size,
            'stopping_bytes': exported['stopping'][z].stat().st_size}


def validate(source, csv_root, report, converter_dir=None, zs=None, targets=None):
    source, csv_root, report = map(Path, (source, csv_root, report))
    data = {'passed': False, 'status': 'running', 'source': str(source.resolve()),
            'csv_root': str(csv_root.resolve()), 'projectiles': 0, 'records': [],
            'total_csv_bytes': 0}
    write_report(report, data)
    try:
        expected = selection(zs, '--z', 29, 92)
        targets = selection(targets, '--target-z', 1, 96)
        original = inventory(source, '.ccbin')
        exported = inventory(csv_root, '.csv')
        if expected is None:
            expected = sorted(original['delta'])
        data['expected_Z'] = expected
        data['selected_target_Z'] = targets
        if set(exported['delta']) != set(expected):
            raise ValueError(f'CSV projectile coverage mismatch: expected {expected}, '
                             f'found {sorted(exported["delta"])}. '
                             'Use --z to explicitly select a smaller export.')
        missing = sorted(set(expected) - original['delta'].keys())
        if missing:
            raise ValueError(f'Requested projectiles absent in source CCPACK: {missing}')
        read_pack = load_read_pack(converter_dir)
        for z in expected:
            record = validate_projectile(z, original, exported, read_pack, targets)
            data['records'].append(record)
            data['projectiles'] = len(data['records'])
            data['total_csv_bytes'] += record['delta_bytes'] + record['stopping_bytes']
            print(f'Z={z}: schema, coverage and all CSV numbers exactly preserved', flush=True)
        data.update(passed=True, status='passed')
        write_report(report, data)
        return data
    except Exception as error:
        data.update(passed=False, status='failed', error_type=type(error).__name__, error=str(error))
        write_report(report, data)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--csv-root', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--converter-dir', type=Path,
                        help='Source directory containing cc_binary.py; independent of --source data')
    parser.add_argument('--z', type=int, nargs='+',
                        help='Expected exported projectiles; default: all projectiles in source')
    parser.add_argument('--target-z', type=int, nargs='+',
                        help='Expected complete stopping target curves; default: all source targets')
    a = parser.parse_args()
    try:
        validate(a.source, a.csv_root, a.report, a.converter_dir, a.z, a.target_z)
    except Exception as error:
        print(f'CSV validation failed: {error}', file=sys.stderr)
        sys.exit(1)

