#!/usr/bin/env python3
"""Refine completed stopping packs without rereading the large GLOBAL CSV.

The source manifest is read-only.  Each replacement is staged completely before
the previous stopping pack is moved into the explicit backup directory.  Per-Z
reports contain the replacement artifact record for a separate manifest update.
Disjoint projectile ranges may be run concurrently.
"""
from __future__ import annotations

import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

import argparse
from datetime import datetime, timezone
import gc
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from cc_binary import KINDS, read_pack, write_pack
from cc_inputs import read_atima
from convert import prepare_stopping


WORKSPACE = Path(__file__).resolve().parents[2]


def emit(value):
    print(json.dumps(value, ensure_ascii=False, allow_nan=False), flush=True)


def checked_destination(path: Path, root: Path) -> Path:
    """Resolve the final path before any move and keep it under its named root."""
    if path.is_symlink() or (hasattr(path, 'is_junction') and path.is_junction()):
        raise ValueError(f'Refusing linked destination: {path}')
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()) or resolved == root.resolve():
        raise ValueError(f'Destination escapes intended directory: {path}')
    return resolved


def verify_record(tables: Path, records: dict, relative: str, kind: str, zp: int):
    if relative not in records:
        raise ValueError(f'Manifest has no artifact {relative}')
    record = records[relative]
    path = checked_destination(tables/relative, tables)
    if not path.is_file() or record.get('kind') != kind or record.get('Zp') != zp:
        raise ValueError(f'Wrong/missing artifact: {relative}')
    if path.stat().st_size != record['bytes']:
        raise ValueError(f'Artifact size mismatch: {relative}')
    return path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, required=True,
                        help='Completed prepared-table directory')
    parser.add_argument('--atima-splines', type=Path, required=True)
    parser.add_argument('--z-min', type=int, required=True)
    parser.add_argument('--z-max', type=int, required=True)
    parser.add_argument('--backup-dir', type=Path, required=True,
                        help='Absolute backup directory under this workspace, outside --tables')
    parser.add_argument('--normalization-tolerance', type=float,
                        help='Default: the tolerance declared in the existing manifest')
    return parser.parse_args()


def main():
    args = parse_args()
    if not 29 <= args.z_min <= args.z_max <= 92:
        raise ValueError('Require 29 <= Z-min <= Z-max <= 92')
    if not args.backup_dir.is_absolute():
        raise ValueError('--backup-dir must be an explicitly supplied absolute path')
    tables = args.tables.resolve()
    backup = args.backup_dir.resolve()
    if not tables.is_relative_to(WORKSPACE) or not backup.is_relative_to(WORKSPACE):
        raise ValueError(f'Table and backup writes must stay under {WORKSPACE}')
    if tables == WORKSPACE or backup == WORKSPACE or backup.is_relative_to(tables) or tables.is_relative_to(backup):
        raise ValueError('Backup must be separate from the tables directory and its ancestors')
    if tables.drive.casefold() != backup.drive.casefold():
        raise ValueError('Backup must be on the same drive for filesystem rename')
    if args.backup_dir.is_symlink() or (hasattr(args.backup_dir, 'is_junction') and args.backup_dir.is_junction()):
        raise ValueError('Backup directory cannot be a symbolic link/junction')
    manifest = json.loads((tables/'manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema') != 'CC_Dirac_AtomicTables_v1' or manifest.get('complete') is not True:
        raise ValueError('Refinement requires a completed compatible source manifest')
    config = manifest['configuration']
    if args.z_min < config['Z_min'] or args.z_max > config['Z_max']:
        raise ValueError('Selected Z range is outside the source manifest')
    tolerance = (args.normalization_tolerance if args.normalization_tolerance is not None
                 else float(config['normalization_relative_tolerance']))
    if not np.isfinite(tolerance) or not 1e-9 <= tolerance <= 1e-2:
        raise ValueError('Tolerance must lie in [1e-9,1e-2]')
    records = {}
    for record in manifest['files']:
        relative = str(record['path']).replace('\\', '/')
        if relative in records:
            raise ValueError(f'Duplicate manifest artifact: {relative}')
        records[relative] = record
    zs = list(range(args.z_min, args.z_max+1))
    started = time.monotonic()
    for zp in zs:
        for kind in ('delta', 'global', 'stopping'):
            verify_record(tables, records, f'{kind}/Z{zp:03d}.ccbin', kind, zp)
        for path, root in [(backup/f'Z{zp:03d}.ccbin', backup),
                           (tables/'refinement_reports'/f'Z{zp:03d}.json', tables),
                           (tables/'refinement_staging'/f'Z{zp:03d}.ccbin', tables),
                           (tables/'refinement_staging'/f'Z{zp:03d}.ccbin.partial', tables)]:
            if checked_destination(path, root).exists():
                raise FileExistsError(f'Refinement will not overwrite an existing artifact: {path}')
    targets = np.asarray(manifest['target_Z'], dtype=np.int32)
    atima = read_atima(args.atima_splines, zs, targets, progress=emit)
    backup.mkdir(parents=True, exist_ok=True)
    (tables/'refinement_reports').mkdir(exist_ok=True)
    emit({'phase':'refinement_inputs_verified', 'Z_min':args.z_min, 'Z_max':args.z_max,
          'tolerance':tolerance, 'backup_directory':str(backup)})
    for zp in zs:
        t0 = time.monotonic()
        relative = f'stopping/Z{zp:03d}.ccbin'
        global_header, global_arrays = read_pack(tables/'global'/f'Z{zp:03d}.ccbin')
        delta_header, delta = read_pack(tables/'delta'/f'Z{zp:03d}.ccbin')
        if global_header != {'kind':KINDS['global'], 'Zp':zp} or delta_header != {'kind':KINDS['delta'], 'Zp':zp}:
            raise ValueError(f'Wrong GLOBAL/Dirac binary identity for Z={zp}')
        if not np.array_equal(global_arrays['target_Z'], targets):
            raise ValueError(f'GLOBAL target list differs from manifest for Z={zp}')
        g = SimpleNamespace(projectile_z=np.array([zp], dtype=np.int32),
                            target_z=global_arrays['target_Z'],
                            target_mass_g_mol=global_arrays['target_mass_g_mol'],
                            energies_MeV_u=global_arrays['energy_MeV_u'],
                            transitions=global_arrays['transitions_N'],
                            sigma_cm2=global_arrays['sigma_cm2'][None, ...])
        stopping, diagnostic = prepare_stopping(g, 0, atima, delta, tolerance)
        staging = checked_destination(tables/'refinement_staging'/f'Z{zp:03d}.ccbin', tables)
        artifact = write_pack(staging, 'stopping', zp, stopping)
        source = verify_record(tables, records, relative, 'stopping', zp)
        previous = checked_destination(backup/f'Z{zp:03d}.ccbin', backup)
        destination = checked_destination(tables/relative, tables)
        if previous.exists():
            raise FileExistsError(previous)
        source.rename(previous)
        try:
            staging.rename(destination)
        except BaseException:
            if not destination.exists():
                previous.rename(destination)
            raise
        artifact['path'] = relative
        report = {'schema':'CC_StoppingRefinement_v1', 'Zp':zp,
                  'created_utc':datetime.now(timezone.utc).isoformat(),
                  'normalization_relative_tolerance':tolerance,
                  'refinement_probe_safety_factor':0.9,
                  'previous_artifact':records[relative], 'backup_path':str(previous),
                  'artifact':artifact, 'normalization':diagnostic,
                  'stopping_energy_nodes':len(stopping['energy_MeV_u']),
                  'elapsed_seconds':time.monotonic()-t0,
                  'manifest_updated':False}
        report_path = checked_destination(tables/'refinement_reports'/f'Z{zp:03d}.json', tables)
        with report_path.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
        emit({'phase':'stopping_refined', 'Zp':zp,
              'stopping_energy_nodes':report['stopping_energy_nodes'],
              'max_probe_mean_error':diagnostic['refinement'][-1]['max_relative_mean_error'],
              'seconds':round(report['elapsed_seconds'], 2), 'report':str(report_path)})
        del global_arrays, delta, g, stopping
        gc.collect()
    emit({'phase':'refinement_complete', 'projectiles':len(zs),
          'seconds':round(time.monotonic()-started, 2), 'manifest_updated':False})


if __name__ == '__main__':
    main()
