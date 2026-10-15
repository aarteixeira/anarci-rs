# Release performance measurements

The release benchmark compared the pinned reference ANARCI, the extended Python
API and direct native bindings on the same permitted inputs. Timings exclude
interpreter startup and the first warmup batch. The exact and pan engines were
measured separately.

Local run: Apple M5 Max, macOS ARM64, CPython 3.10.20, 32 distinct synthetic
variants of the identified anti-lysozyme antibody pair. Settings were IMGT,
human/mouse, threshold 80, germline assignment enabled with identity matching,
one warmup and three timed repeats. Throughput uses the median wall time.

| Workers | Reference ANARCI | Rust exact through Python | Rust pan through Python |
|---|---:|---:|---:|
| 1 | 59.44 sequences/s | 107.55 sequences/s | 535.44 sequences/s |
| 4 | 105.75 sequences/s | 428.03 sequences/s | 2042.31 sequences/s |

These rates describe this small synthetic workload. Pan scans seven profiles,
while exact/reference scan 29; the speed comparison does not establish equivalent
numbering or score distributions. The evalue germline method, which is the pan
default, was outside this matched identity-assignment timing comparison.

A separate repeated-sequence workload measured deduplication benefits. It is
recorded separately from the distinct-input results. Interleaved, paired calls
in the same process compared the public Python wrapper with direct native
bindings. Across eight configurations, the median wrapper/direct time ratios
ranged from 0.963 to 1.002. This run found no material wrapper overhead; separate
process timings were noisier and were not used to estimate that overhead.

The machine-readable record contains each duration, input hashes, exact settings,
worker counts, interpreter versions and any species-filter diagnostics. It is
delivered beside the release artifacts as `performance-results.json`.

Reproduce the measurements with separate installed environments:

```bash
python scripts/benchmark_release.py \
  --native-python /path/to/native-env/bin/python \
  --reference-python /path/to/reference-env/bin/python \
  --count 32 --warmup 1 --repeat 3 --workers 1 4 \
  --output performance-results.json
```

The reference environment must match the hashes in `reference-version.json` and
have its HMMER executable available. The native environment uses the built wheel.
For machines that restrict CPU metadata queries, pass `--cpu-model` with an
independently verified model name. The benchmark records this supplied value.
