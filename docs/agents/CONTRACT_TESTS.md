# Executable contract tests

Run the focused Rust targets separately from suites that open historical fixtures:

```sh
cargo test --offline --locked -p anarci-core --test property_contract
cargo test --offline --locked -p anarci-learned --no-default-features --test property_contract
python -m unittest discover -s tests/python -p test_properties.py -v
```

Use the interpreter where the current package and native extension are installed. Like the other installed Python tests, this suite imports `anarci` and `anarci_rs` normally. It does not add the source tree to `sys.path`, build an extension, or skip a failed import. Record both packages' `__file__` values when reporting results: an old installed extension does not validate a new checkout.

The Rust targets need cached Cargo dependencies, but neither target depends on the native HMMER crate. The learned target needs no Candle feature, checkpoint, GPU, or training run. Stop if offline dependencies are missing; do not fetch them. Coordinate Cargo jobs with the integrator before running them so workers do not duplicate builds or wait on the same target lock. These commands do not authorize a native build, package installation, or download.

## Contracts exercised

[The core target](../../crates/anarci-core/tests/property_contract.rs) contains seven tests:

- All 256 byte values at each of three query positions agree with the scanner alphabet, with one-based error positions. Lengths around the 100,000-residue limit exercise both sides of the boundary.
- Invalid names, bytes, and lengths at every slot of a four-row batch fail before either sequential or parallel orchestration calls the scanner. Duplicate valid rows also exercise validation before deduplication. Invalid numeric settings and schemes fail on empty batches.
- Checked scanner errors retain the engine-error variant; successful empty scans retain the no-hit shape. The fake scanner panics if orchestration uses the unchecked entry point.
- Generated batches of 0, 1, 2, 7, and 31 rows scan each distinct sequence once, restore input order and query names, retain per-query scores and coordinates, and agree with sequential numbering. Mutating one replayed result cannot change its duplicate.
- Every scheme abbreviation resolves in both cases; malformed scheme strings return errors.
- J-rescue rejects inconsistent state/detail counts and out-of-bounds sequence indices without scanning.
- Malformed public alignment states must return errors without panicking. This is a regression contract, not permission to repair input silently.

[The learned target](../../crates/anarci-learned/tests/property_contract.rs) contains seven tests. It checks all 261 vocabulary entries and exhausts 6,561 small canonical domain shapes with matches, deletes, and insertions. Round trips retain nonzero query origins, insertion owners, and residue counts. Boundary-delete normalization is idempotent. All five separator labels and abutting match resets retain absolute query indices.

The raw Rust codec is permissive: its malformed-input test checks bounded behavior over all 256 possible state IDs, including the existing last-match-wins behavior for duplicate matches. It does not claim these inputs are valid biological alignments. Out-of-vocabulary labels must not be passed to the unchecked Rust vocabulary mapper. Python validation is a separate boundary contract.

Decoder tests recover confident valid paths, repeat 24 short matrices from a fixed LCG seed (`0x5eed`), and compare reconstructed residues with the emitted labels in query order. A targeted decreasing-insertion-owner case checks the same order invariant. Merely counting residues would miss this defect.

[The Python target](../../tests/python/test_properties.py) contains seven tests. Mocked native calls check that the compatibility facade forwards `ALL` and human/mouse defaults while the extended facade forwards `pan` and unrestricted species. Invalid facade options must fail before native calls. Generated short sequences check the compatibility facade's silence and the extended facade's warning. Both facades must reject malformed states, invalid scan input, and non-finite thresholds as Python exceptions rather than reporting no hit.

Synthetic alignment geometry checks six schemes, upper-case spellings, three query offsets, preserved residue strings, and inclusive coordinates. Installed native integration compares duplicate-containing batches with singleton calls at one and two workers, including no-hit rows, query metadata, and independently mutable results. Inputs come only from [the existing safe helper](../../tests/python/safe_reference.py): public 12e8 heavy/light examples, synthetic suffixes, and an alanine no-hit sequence. No historical fixtures or additional native data files are opened by the tests.

## Failure policy and limits

Keep a failing invariant visible and report the reproducer. Do not change an expected value to match a panic, reordered sequence, or swallowed scanner error. Production fixes need separate ownership and review. The alignment path now validates states before scheme dispatch and again after insertion smoothing, including cases where smoothing could synthesize an invalid owner. The constrained decoder retains the last insertion owner as part of its hidden state; the segmenter recognizes resets below that owner. These changes preserve query order for the generated contracts instead of accepting a reordered expected result.

These are bounded property tests over explicit finite sets and deterministic generated inputs. They add no property-testing dependency and use no coverage-guided fuzz instrumentation. Passing them does not establish branch coverage, memory safety of native code, biological numbering accuracy, compatibility outside the tested contracts, or trained-model accuracy. The fake scanner cannot validate HMMER itself. The native Python checks require an installed extension with bundled profiles; they exercise exact `ALL` mode explicitly and make no pan-mode accuracy claim. There are no optional correctness skips in this file.
