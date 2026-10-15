#!/usr/bin/env python3
"""Measure installed native/reference APIs with matched permitted inputs/settings.

Each interpreter runs separately; warmups and timing exclude startup/imports.
The JSON includes each duration and dataset hash. Pan and ALL use distinct model
sets, so pan/reference ratios describe throughput, not an accuracy equivalence.
"""
import argparse
import contextlib
import hashlib
import json
import io
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests/python"))
from safe_reference import dataset, dataset_sha256


def cpu_model():
    if sys.platform == "darwin":
        result = subprocess.run(["/usr/sbin/sysctl", "-n", "machdep.cpu.brand_string"], text=True, capture_output=True)
        return result.stdout.strip() if result.returncode == 0 else "unavailable: " + result.stderr.strip()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                return line.partition(":")[2].strip()
    return platform.processor()


def worker(args):
    payload = json.loads(sys.stdin.read())
    records = [tuple(record) for record in payload["records"]]
    if args.backend == "reference":
        import anarci as package
        import importlib.metadata
        distribution = importlib.metadata.version("anarci")
        source_hashes = {}
        base = Path(package.__file__).parent
        for name in ("__init__.py", "anarci.py", "schemes.py", "germlines.py"):
            source_hashes[name] = hashlib.sha256((base / name).read_bytes()).hexdigest()
        if source_hashes != payload["reference"]["source_sha256"] or distribution != payload["reference"]["version"]:
            raise RuntimeError("reference differs from pinned source hashes/version")
        settings = {"scheme": "imgt", "database": "ALL", "assign_germline": True,
                    "allowed_species": ["human", "mouse"], "bit_score_threshold": 80, "ncpu": args.workers}
    else:
        import anarci_rs as package
        distribution = package.__version__
        base = Path(package.__file__).parent
        source_hashes = {name: hashlib.sha256((base / name).read_bytes()).hexdigest()
                         for name in ("__init__.py", "_api.py")}
        settings = {"scheme": "imgt", "database": args.database, "assign_germline": True,
                    "allowed_species": ["human", "mouse"], "bit_score_threshold": 80,
                    "germline_method": "identity", "ncpu": args.workers}
    target = package._native.run_anarci if args.backend == "native-direct" else package.run_anarci
    diagnostics = io.StringIO()
    def invoke(function=target):
        with contextlib.redirect_stdout(diagnostics):
            return function(records, **settings)
    for _ in range(args.warmup):
        invoke()
    durations = []
    paired_direct = []
    for repeat_index in range(args.repeat):
        functions = [("api", target)]
        if args.backend == "native":
            functions.append(("direct", package._native.run_anarci))
            if repeat_index % 2:
                functions.reverse()
        timings = {}
        for label, function in functions:
            started = time.perf_counter()
            result = invoke(function)
            timings[label] = time.perf_counter() - started
        durations.append(timings["api"])
        if "direct" in timings:
            paired_direct.append(timings["direct"])
        if len(result[0]) != len(records) or not all(result[1]):
            raise RuntimeError("benchmark work was not completed")
    median = statistics.median(durations)
    print(json.dumps({"backend": args.backend, "database": args.database, "workers": args.workers,
                      "python": sys.version, "executable": sys.executable, "package_version": distribution,
                      "package_source_sha256": source_hashes,
                      "native_extension_sha256": hashlib.sha256(Path(package._native.__file__).read_bytes()).hexdigest() if args.backend != "reference" else None,
                      "settings": settings, "sequences": len(records), "unique_sequences": len({s for _, s in records}),
                      "dataset_sha256": dataset_sha256(records), "warmup": args.warmup, "repeat": args.repeat,
                      "diagnostic_lines": sorted(set(diagnostics.getvalue().splitlines())),
                      "paired_direct_seconds": paired_direct,
                      "paired_wrapper_to_direct_ratio": statistics.median(a / b for a, b in zip(durations, paired_direct)) if paired_direct else None,
                      "seconds": durations, "median_seconds": median, "sequences_per_second": len(records) / median}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--reference-python", type=Path, help="separate environment with the pinned ANARCI 1.3 and hmmscan")
    parser.add_argument("--count", type=int, default=16)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, min(4, os.cpu_count() or 1)])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cpu-model", help="verified hardware model when local system queries are restricted")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--backend", choices=("native", "native-direct", "reference"), help=argparse.SUPPRESS)
    parser.add_argument("--database", default="ALL", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.count < 2 or args.repeat < 1 or args.warmup < 1 or any(w < 1 for w in args.workers):
        parser.error("count >= 2, repeat/warmup >= 1 and positive workers are required")
    if args.worker:
        args.workers = args.workers[0]
        worker(args)
        return
    if not args.reference_python:
        parser.error("--reference-python is required to compare against the pinned reference")
    reference = json.loads((ROOT / "docs/reference-version.json").read_text())
    evidence = {"platform": platform.platform(), "machine": platform.machine(), "processor": args.cpu_model or cpu_model(),
                "processor_source": "command-line --cpu-model" if args.cpu_model else "local system query",
                "cpu_count": os.cpu_count(), "reference_pin": reference, "runs": [],
                "method": "Separate interpreters; same batch records/threshold/scheme/species/identity assignment; one warmup minimum; per-batch wall time; reference run_anarci uses multiprocessing and native run_anarci uses Rayon; repeated sequences may benefit from native deduplication. Pan uses seven models; ALL/reference use 29. Native paired wrapper/direct calls alternate order within one process; their per-pair ratios assess Python overhead."}
    with tempfile.TemporaryDirectory(prefix="anarci-benchmark-") as folder:
        for unique in (True, False):
            records = dataset(args.count, unique=unique)
            for workers in dict.fromkeys(args.workers):
                for backend, database in (("reference", "ALL"), ("native", "ALL"), ("native-direct", "ALL"), ("native", "pan"), ("native-direct", "pan")):
                    executable = args.reference_python if backend == "reference" else args.native_python
                    executable = executable.absolute()
                    command = [str(executable), str(Path(__file__).resolve()), "--worker", "--backend", backend,
                               "--database", database, "--workers", str(workers), "--count", str(args.count),
                               "--warmup", str(args.warmup), "--repeat", str(args.repeat)]
                    environment = os.environ.copy()
                    environment.pop("PYTHONPATH", None)
                    environment["PATH"] = str(executable.parent) + os.pathsep + environment.get("PATH", "")
                    process = subprocess.run(command, cwd=folder, env=environment, input=json.dumps({"records": records, "reference": reference}), text=True,
                                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)
                    if process.returncode:
                        raise RuntimeError("%s/%s workers=%d failed:\n%s" % (backend, database, workers, process.stderr))
                    run = json.loads(process.stdout)
                    run["dataset_kind"] = "unique-synthetic-variants" if unique else "two-repeated-sequences"
                    evidence["runs"].append(run)
                    print("%s %-13s %-3s workers=%d: %.2f sequences/s" % (run["dataset_kind"], backend, database, workers, run["sequences_per_second"]), flush=True)
    evidence["status"] = "measured"
    evidence["wrapper_to_direct_median_ratios"] = []
    for run in evidence["runs"]:
        if run["backend"] != "native":
            continue
        direct = next(candidate for candidate in evidence["runs"] if candidate["backend"] == "native-direct"
                      and all(candidate[key] == run[key] for key in ("database", "workers", "dataset_sha256")))
        evidence["wrapper_to_direct_median_ratios"].append({
            "database": run["database"], "workers": run["workers"], "dataset_kind": run["dataset_kind"],
            "separate_process_ratio": run["median_seconds"] / direct["median_seconds"],
            "paired_interleaved_ratio": run["paired_wrapper_to_direct_ratio"]})
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    else:
        print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    main()
