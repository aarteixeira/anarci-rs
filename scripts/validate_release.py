#!/usr/bin/env python3
"""Inspect release archives, then test isolated wheel and source installations.

This script reads archive member names before reading approved package metadata.
It never opens repository data/workbench files and does not publish artifacts.
"""
import argparse
import configparser
import email
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import venv
import zipfile

REPO = Path(__file__).resolve().parents[1]
SAFE_TEST_FILES = {"safe_reference.py", "safe_expected.json", "safe_consensus.json", "safe_alignment.json", "safe_helpers.json", "safe_hmmer.txt", "test_api.py",
                   "test_compat.py", "test_experimental.py", "test_integration.py", "test_console_entry.py", "test_optional_model.py", "test_properties.py"}
SDIST_ROOT = {"Cargo.toml", "Cargo.lock", "pyproject.toml", "PKG-INFO", "LICENSE", "NOTICE.md", "README.md"}
SDIST_DOCS = {"BUILDING.md", "RELEASING.md", "reference-version.json",
              "API.md", "COMPATIBILITY.md", "FEATURE_COVERAGE.md", "EXPERIMENTAL.md"}
CRATES = {"anarci-core", "anarci-hmm", "anarci-py", "anarci-cli", "anarci-learned", "hmmer-sys"}
WHEEL_PYTHON = {
    "anarci": {"__init__.py", "anarci.py", "schemes.py", "germlines.py"},
    "anarci_rs": {"__init__.py", "__init__.pyi", "py.typed", "_api.py", "compat.py", "io.py", "cli.py", "experimental.py", "schemes.py", "germlines.py"},
}


def allowed_source(path):
    parts = path.parts
    if len(parts) == 1:
        return parts[0] in SDIST_ROOT
    if parts[0] == "licenses":
        return len(parts) == 2 and parts[1].endswith(".txt")
    if parts[0] == "docs":
        return len(parts) == 2 and parts[1] in SDIST_DOCS
    if parts[0] == "python":
        return len(parts) == 3 and parts[1] in WHEEL_PYTHON and parts[2] in WHEEL_PYTHON[parts[1]]
    if parts[0] == "tests":
        return len(parts) == 3 and (
            parts[1] == "python" and parts[2] in SAFE_TEST_FILES
            or parts[1] == "packaging" and parts[2] == "test_source_extraction.py")
    if parts[0] == "scripts":
        return len(parts) == 2 and parts[1] in {"validate_release.py", "benchmark_release.py", "generate_rust_notices.py"}
    if parts[0] == "reference_data":
        return str(path) in {"reference_data/dat/HMMs/ALL.hmm", "reference_data/dat/HMMs/FEW.hmm"}
    if parts[0] == "crates" and len(parts) >= 3 and parts[1] in CRATES:
        rest = parts[2:]
        if len(rest) == 1:
            return rest[0] in {"Cargo.toml", "build.rs", "wrapper.h", "README.md"}
        if rest[0] in {"src", "tests", "examples"}:
            return path.suffix == ".rs"
        return parts[1] == "anarci-core" and rest == ("data", "germlines.json")
    return False


def safe_member(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name:
        raise ValueError("unsafe archive path: " + name)
    return path


def inspect_archive(path):
    path = Path(path).resolve()
    required_licenses = {"LICENSE", "NOTICE.md", "ANARCI-BSD-3-Clause.txt", "HMMER-LICENSE.txt", "Easel-LICENSE.txt", "CC-BY-4.0.txt", "IMGT-DATA-TERMS.txt", "RUST-THIRD-PARTY.txt"}
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            for name in names:
                member = safe_member(name)
                if name.endswith("/"):
                    continue
                first = member.parts[0]
                wheel_metadata = (first.endswith(".dist-info") and (
                    len(member.parts) == 2 and member.name in {"METADATA", "WHEEL", "RECORD", "entry_points.txt"}
                    or len(member.parts) >= 3 and member.parts[1] == "licenses" and member.name in required_licenses
                    or len(member.parts) == 3 and member.parts[1:] == ("sboms", "anarci-py.cyclonedx.json")))
                package_member = (len(member.parts) == 2 and first in WHEEL_PYTHON and
                                  (member.name in WHEEL_PYTHON[first] or first == "anarci_rs" and
                                   member.name.startswith("_native.") and member.suffix in {".so", ".pyd"}))
                permitted = package_member or wheel_metadata
                if not permitted:
                    raise ValueError("unexpected wheel member: " + name)
            metadata_names = [n for n in names if n.endswith(".dist-info/METADATA")]
            if len(metadata_names) != 1:
                raise ValueError("wheel requires one METADATA")
            metadata = email.message_from_bytes(archive.read(metadata_names[0]))
            python_bounds = {part.strip() for part in (metadata["Requires-Python"] or "").split(",")}
            if metadata["Name"] != "anarci-rs" or python_bounds != {">=3.9", "<3.14"}:
                raise ValueError("unexpected distribution metadata")
            if not any(n.startswith("anarci_rs/_native") and n.endswith((".so", ".pyd")) for n in names):
                raise ValueError("missing native extension")
            entry_points = configparser.ConfigParser()
            for name in names:
                if name.endswith("entry_points.txt"):
                    entry_points.read_string(archive.read(name).decode())
            if entry_points.get("console_scripts", "anarci-rs", fallback="") != "anarci_rs.cli:main":
                raise ValueError("missing CLI entry point")
            present = {PurePosixPath(n).name for n in names}
            missing = required_licenses - present
            if missing:
                raise ValueError("wheel missing notices: " + ", ".join(sorted(missing)))
    else:
        with tarfile.open(path, "r:*") as archive:
            members = archive.getmembers()
            roots = {safe_member(m.name).parts[0] for m in members}
            if len(roots) != 1:
                raise ValueError("sdist requires one root directory")
            names = []
            for member in members:
                parts = safe_member(member.name).parts[1:]
                if member.isdir():
                    continue
                if not member.isfile():
                    raise ValueError("sdist links/special files are not allowed")
                relative = PurePosixPath(*parts)
                if not allowed_source(relative):
                    raise ValueError("unexpected sdist member: " + str(relative))
                names.append(str(relative))
            required = {"pyproject.toml", "Cargo.toml", "Cargo.lock", "reference_data/dat/HMMs/ALL.hmm",
                        "reference_data/dat/HMMs/FEW.hmm", "python/anarci/__init__.py", "python/anarci_rs/__init__.py"}
            if not required.issubset(names):
                raise ValueError("sdist is missing required build inputs: " + str(sorted(required - set(names))))
            if not required_licenses.issubset({PurePosixPath(n).name for n in names}):
                raise ValueError("sdist is missing license notices")
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "members": len(names)}



def extract_source_archive(archive_path, destination):
    """Extract an allowlisted source archive as regular files, on Python 3.9+."""
    inspect_archive(archive_path)
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    with tarfile.open(archive_path, "r:*") as archive:
        members = archive.getmembers()
        root = safe_member(members[0].name).parts[0]
        for member in members:
            if not member.isfile():
                continue
            output = destination.joinpath(*safe_member(member.name).parts)
            output.parent.mkdir(parents=True, exist_ok=True)
            with archive.extractfile(member) as source, output.open("wb") as target:
                shutil.copyfileobj(source, target)
    return destination / root


def run(command, cwd, env=None):
    started = time.perf_counter()
    process = subprocess.run(command, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(process.stdout, end="", flush=True)
    if process.returncode:
        raise RuntimeError("command failed (%s): %s" % (process.returncode, command))
    return {"command": [str(c) for c in command], "seconds": time.perf_counter() - started, "returncode": process.returncode}


def installed_checks(artifact, work, tests):
    environment = work / "venv"
    # Match the venv CLI on POSIX: copied standalone Python launchers can lose
    # their relative libpython lookup, while a symlink retains the base runtime.
    venv.EnvBuilder(with_pip=True, symlinks=os.name != "nt").create(environment)
    executable = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    records = [run([str(executable), "-m", "pip", "install", "--no-deps", str(artifact)], work)]
    isolated = os.environ.copy()
    isolated.pop("PYTHONPATH", None)
    isolated.pop("PYTHONHOME", None)
    isolated["PATH"] = str(executable.parent)  # proves no external hmmscan is needed at runtime
    records.append(run([str(executable), "-I", "-m", "unittest", "discover", "-s", str(tests), "-p", "test_*.py", "-v"], work, isolated))
    records.append(run([str(executable), "-I", "-c", "import anarci,anarci_rs,sys; from pathlib import Path; assert all(Path(m.__file__).is_relative_to(sys.prefix) for m in (anarci,anarci_rs)); print(anarci_rs.__file__)"], work, isolated))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path)
    parser.add_argument("--build-sdist", action="store_true", help="also build the supplied sdist in isolation and test the resulting wheel")
    parser.add_argument("--inspect-only", action="store_true")
    parser.add_argument("--output", type=Path, help="write a JSON evidence record")
    args = parser.parse_args()
    if args.build_sdist and not args.sdist:
        parser.error("--build-sdist requires --sdist")
    evidence = {"python": sys.version, "platform": platform.platform(), "artifacts": [], "checks": []}
    evidence["artifacts"].append(inspect_archive(args.wheel))
    if args.sdist:
        evidence["artifacts"].append(inspect_archive(args.sdist))
    if not args.inspect_only:
        with tempfile.TemporaryDirectory(prefix="anarci-release-") as temporary:
            folder = Path(temporary)
            tests = folder / "safe-tests"
            tests.mkdir()
            for name in SAFE_TEST_FILES:
                source = REPO / "tests/python" / name
                if source.is_file():
                    shutil.copy2(source, tests / name)
            wheel_work = folder / "wheel-install"
            wheel_work.mkdir()
            evidence["checks"].extend(installed_checks(args.wheel.resolve(), wheel_work, tests))
            if args.build_sdist:
                build_work = folder / "source-build"
                build_work.mkdir()
                local_sdist = build_work / args.sdist.name
                shutil.copy2(args.sdist, local_sdist)
                wheels = build_work / "wheels"
                wheels.mkdir()
                # pip extracts the validated archive into a separate temporary source tree.
                build_environment = os.environ.copy()
                build_environment.pop("PYTHONPATH", None)
                build_environment["CARGO_TARGET_DIR"] = str(build_work / "cargo-target")
                evidence["checks"].append(run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-cache-dir", "--wheel-dir", str(wheels), str(local_sdist)], build_work, build_environment))
                built = list(wheels.glob("*.whl"))
                if len(built) != 1:
                    raise RuntimeError("expected one wheel from source archive")
                evidence["artifacts"].append(inspect_archive(built[0]))
                install_work = folder / "source-install"
                install_work.mkdir()
                evidence["checks"].extend(installed_checks(built[0], install_work, tests))
    evidence["status"] = "passed"
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
