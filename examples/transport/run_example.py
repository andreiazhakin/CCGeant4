"""Run the silicon transport example for Cu, Xe, U or Fe."""

import argparse
import csv
import lzma
import os
from pathlib import Path
import shutil
import subprocess


EXAMPLE = Path(__file__).resolve().parent
CASES = {
    "cu": (29, 63, 400, "msc"),
    "xe": (54, 124, 400, "msc"),
    "u": (92, 238, 400, "msc"),
    "ss_cu": (29, 63, 200, "tabulated_ss"),
    "ss_native": (29, 63, 200, "ss"),
    "fe_ss": (26, 56, 200, "tabulated_ss"),
}


def unpack_tables(destination):
    for archive in sorted((EXAMPLE / "data").rglob("*.xz")):
        target = destination / archive.relative_to(EXAMPLE / "data").with_suffix("")
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(target.suffix + ".tmp")
        try:
            with lzma.open(archive, "rb") as source, temporary.open("wb") as output:
                shutil.copyfileobj(source, output)
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--geant4-bin", type=Path, help="DLL directory of the matching Windows installation")
    parser.add_argument("--tables", type=Path, help="Prepared tables directory; otherwise unpack the supplied data")
    parser.add_argument("--elastic-tables", type=Path, help="External elastic tables index for tabulated_ss")
    parser.add_argument("--output", type=Path, default=EXAMPLE / "results")
    parser.add_argument("--events", type=int, default=1000, help="Events per case")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=6110101)
    parser.add_argument("--case", choices=CASES, nargs="+", default=["cu", "xe", "u"])
    parser.add_argument("--model", choices=("dirac", "legacy", "emopt4"), default="dirac")
    args = parser.parse_args()
    if args.events < 1 or args.threads < 1 or args.seed < 1:
        parser.error("--events, --threads and --seed must be positive")

    executable = args.executable.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    tables = args.tables.resolve() if args.tables else None
    if tables is None and args.model != "emopt4":
        tables = unpack_tables(output / "tables")
    elastic = args.elastic_tables.resolve() if args.elastic_tables else None
    environment = os.environ.copy()
    if args.geant4_bin:
        environment["PATH"] = str(args.geant4_bin.resolve()) + os.pathsep + environment.get("PATH", "")

    fields = ["case", "model", "scattering", "Z", "A", "charge", "energy_MeV_u",
              "material", "thickness_um", "cut_um", "events", "threads", "seed",
              "tables", "elastic_tables", "macro", "event_csv", "log", "returncode"]
    with (output / "configurations.csv").open("w", newline="", encoding="utf-8") as record:
        writer = csv.DictWriter(record, fieldnames=fields)
        writer.writeheader()
        for name in dict.fromkeys(args.case):
            z, a, energy, scattering = CASES[name]
            stem = f"{args.model}_{name}"
            macro = output / f"{stem}.mac"
            events = output / f"{stem}.csv"
            log = output / f"{stem}.log"
            macro.write_text(
                "/example/detector/material G4_Si\n"
                "/example/detector/thickness 100 um\n"
                "/run/initialize\n"
                "/gun/particle ion\n"
                f"/gun/ion {z} {a} {z}\n"
                f"/gun/energy {a * energy} MeV\n"
                f"/run/beamOn {args.events}\n", encoding="utf-8")
            command = [str(executable), "--model", args.model, "--scattering", scattering,
                       "--cut-um", "25", "--threads", str(args.threads), "--seed", str(args.seed),
                       "--output", str(events), "--macro", str(macro)]
            if tables is not None:
                command += ["--tables", str(tables)]
            if elastic is not None and scattering == "tabulated_ss":
                command += ["--elastic-tables", str(elastic)]
            print(f"{stem}: {args.events} events, {args.threads} thread(s)", flush=True)
            with log.open("w", encoding="utf-8") as stream:
                result = subprocess.run(command, cwd=output, env=environment,
                                        stdout=stream, stderr=subprocess.STDOUT)
            writer.writerow(dict(case=name, model=args.model, scattering=scattering,
                Z=z, A=a, charge=z, energy_MeV_u=energy, material="G4_Si", thickness_um=100,
                cut_um=25, events=args.events, threads=args.threads, seed=args.seed,
                tables=str(tables) if tables else "", elastic_tables=str(elastic) if elastic else "",
                macro=str(macro), event_csv=str(events), log=str(log), returncode=result.returncode))
            record.flush()
            if result.returncode:
                raise SystemExit(f"{stem} exited with code {result.returncode}; see {log}")
            print(f"Saved {events}", flush=True)


if __name__ == "__main__":
    main()
