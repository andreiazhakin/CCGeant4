"""Read one projectile lazily and print per-atom stopping and delta moments."""
import argparse
import csv
from pathlib import Path
import sys
import numpy as np
from cc_binary import read_pack
from cc_physics import delta_moments, stopping_at_energy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, required=True)
    parser.add_argument('--Z', type=int, default=54)
    parser.add_argument('--target-Z', type=int, default=14)
    parser.add_argument('--N', type=int, nargs='+', default=[0,1,15,27])
    parser.add_argument('--energy', type=float, nargs='+', default=[30,100,400,1000,2000])
    args = parser.parse_args()
    if any(n < 0 or n > 27 for n in args.N):
        parser.error('N must be in 0..27')
    _, delta = read_pack(args.tables/'delta'/f'Z{args.Z:03d}.ccbin')
    _, stop = read_pack(args.tables/'stopping'/f'Z{args.Z:03d}.ccbin')
    targets = np.flatnonzero(stop['target_Z']==args.target_Z)
    if len(targets) != 1:
        parser.error('Target element is not covered')
    full = stopping_at_energy(delta, stop, args.energy, args.Z)[int(targets[0])]
    moments = delta_moments(delta, args.energy, args.Z)
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['Zp','Zt','N','q','E_MeV_u','full_stopping_MeV_cm2_per_atom',
                     'hard_stopping_MeV_cm2_per_atom','delta_sigma_cm2_per_atom_above_1keV'])
    for ei, energy in enumerate(args.energy):
        for n in args.N:
            writer.writerow([args.Z,args.target_Z,n,args.Z-n,energy,full[ei,n],
                             args.target_Z*moments[ei,n,1],args.target_Z*moments[ei,n,0]])


if __name__ == '__main__':
    main()
