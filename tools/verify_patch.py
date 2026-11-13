"""Compare the applied Geant4 patch with the distributed source files."""
import argparse
from contextlib import nullcontext
import json
from pathlib import Path, PurePosixPath
import subprocess
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PATCH = ROOT / "patches" / "ccgeant4-geant4-11.4.0.patch"


def lf(data):
    return data.replace(b"\r\n", b"\n")


def verify(base):
    manifest = json.loads((ROOT / "patches" / "BASE.json").read_text(encoding="utf-8"))
    changes = ROOT / "source_changes"
    entries = manifest["files"]
    paths = {entry["path"] for entry in entries}
    if len(paths) != len(entries) or not paths:
        raise ValueError("Manifest has duplicate paths or no entries")
    for name in paths:
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or ":" in name or "\\" in name:
            raise ValueError(f"Unexpected relative path in manifest: {name}")
    files_present = {f.relative_to(changes).as_posix() for f in changes.rglob("*") if f.is_file()}
    if files_present != paths:
        raise ValueError(f"source_changes and manifest differ: {sorted(files_present ^ paths)}")

    if base.is_dir():
        source = nullcontext(None)
    else:
        source = zipfile.ZipFile(base)
    with source as archive, tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        names = set(archive.namelist()) if archive is not None else set()
        for entry in entries:
            name = entry["path"]
            if archive is not None:
                member = "geant4-v11.4.0/" + name
                original = archive.read(member) if member in names else None
            else:
                member = base / name
                original = member.read_bytes() if member.is_file() else None
            if entry["status"] == "modified":
                if original is None:
                    raise ValueError(f"Geant4 11.4.0 input is missing: {name}")
                target = work / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(lf(original))
            elif entry["status"] != "added" or original is not None:
                raise ValueError(f"New-file status or upstream presence differs: {name}")

        cmd = ["git", "-c", "core.autocrlf=false", "apply", "--no-index", "--whitespace=nowarn"]
        subprocess.run([*cmd, "--check", str(PATCH)], cwd=work, check=True)
        subprocess.run([*cmd, str(PATCH)], cwd=work, check=True)
        produced = {f.relative_to(work).as_posix() for f in work.rglob("*") if f.is_file()}
        if produced != paths:
            raise ValueError(f"Applied patch and manifest differ: {sorted(produced ^ paths)}")
        for entry in entries:
            if lf((work / entry["path"]).read_bytes()) != lf((changes / entry["path"]).read_bytes()):
                raise ValueError(f"Patched source differs: {entry['path']}")
    modified = sum(e["status"] == "modified" for e in entries)
    print(f"PASS: {modified} modified + {len(entries)-modified} added = {len(entries)} files on Geant4 11.4.0")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True, help="Geant4 11.4.0 source root or CERN release ZIP")
    verify(parser.parse_args().base.resolve())
