"""Extract the versioned runtime tables into a new directory."""
import argparse
import json
import shutil
import sys
import uuid
import zipfile
from pathlib import Path, PurePosixPath, PureWindowsPath


def validate_names(archives):
    seen = set()
    for archive in archives:
        name = archive["name"]
        if Path(name).name != name or PureWindowsPath(name).name != name:
            raise ValueError(f"Invalid archive name: {name}")
        for row in archive["files"]:
            name = row["path"]
            posix = PurePosixPath(name)
            windows = PureWindowsPath(name)
            if (posix.is_absolute() or windows.anchor or not posix.parts
                    or name != posix.as_posix() or "\\" in name or ":" in name
                    or ".." in posix.parts or name.casefold() in seen):
                raise ValueError(f"Invalid or duplicate member path: {name}")
            seen.add(name.casefold())


def read_archive(directory, record, output):
    name = record["name"]
    print(name, flush=True)
    members = [row["path"] for row in record["files"]]
    with zipfile.ZipFile(directory / name) as z:
        if output is not None:
            z.extractall(output, members=members)
            return
        for member in members:
            with z.open(member) as source:
                while source.read(1024 * 1024):
                    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True, help="directory containing the ZIP files")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--verify-only", action="store_true", help="check ZIP contents without extracting")
    action.add_argument("--output", type=Path, help="new directory for the expanded runtime data")
    parser.add_argument("--ionization-only", action="store_true", help="omit the optional silicon elastic bank")
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name("runtime-manifest.json").read_text(encoding="utf-8"))
    if manifest["format"] != "CCGeant4-runtime-data-v1":
        raise ValueError("Unsupported runtime data manifest")
    archives = [r for r in manifest["archives"] if not args.ionization_only or r["group"] == "ionization"]
    legacy_table = Path(__file__).with_name("legacy_q2.csv")
    if not archives:
        raise ValueError("No archives selected")
    validate_names(archives)
    for record in archives:
        if not (args.archives / record["name"]).is_file():
            raise FileNotFoundError(args.archives / record["name"])

    dest = None
    staging = None
    try:
        if args.output is not None:
            dest = args.output.absolute()
            if dest.exists() or dest.is_symlink():
                raise FileExistsError(f"Choose a new output directory: {dest}")
            dest.parent.mkdir(parents=True, exist_ok=True)
            needed = sum(row["bytes"] for r in archives for row in r["files"]) + legacy_table.stat().st_size
            if shutil.disk_usage(dest.parent).free < needed:
                raise OSError(f"Need {needed:,} bytes of free space")
            candidate = dest.parent / f".{dest.name}.partial-{uuid.uuid4().hex}"
            candidate.mkdir()
            staging = candidate
        for record in archives:
            read_archive(args.archives, record, staging)
        if staging is not None:
            shutil.copyfile(legacy_table, staging / "legacy_q2.csv")
            if dest.exists() or dest.is_symlink():
                raise FileExistsError(dest)
            staging.rename(dest)
            staging = None
    finally:
        if staging is not None:
            shutil.rmtree(staging)
    count = sum(len(r["files"]) for r in archives)
    print(f"Read {count} files." if dest is None else f"Installed {count + 1} files in {dest}.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as error:
        sys.exit(str(error))
