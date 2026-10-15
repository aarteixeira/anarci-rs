# anarci-rs installation and builds

anarci-rs provides antigen receptor numbering through a Rust backend, an extended
`anarci_rs` Python API, an `anarci` compatibility package and the `anarci-rs` command.
HMMER 3.4 is linked into the extension. The ALL and FEW reference profiles and
germline reference data are bundled.

## Install a wheel

Use CPython 3.9 through 3.13 in an isolated environment. Install a wheel for your
interpreter and platform:

```sh
python -m pip install /path/to/anarci_rs-0.1.0-<wheel-tags>.whl
python -c 'import anarci_rs, anarci; print(anarci_rs.__name__, anarci.__name__)'
anarci-rs --help
```

The package supplies the `anarci` import namespace. Use a separate environment from
the original ANARCI distribution to keep ownership of that namespace clear.
Wheels contain the native backend and need no separate HMMER installation.

## Build from source

Source builds need Rust/Cargo, a C compiler, make, curl, tar, Python development
headers and libclang for bindgen. On Debian/Ubuntu, the native prerequisites are
`build-essential curl clang libclang-dev python3-dev`. On macOS, install the Xcode
command-line tools and LLVM; set `LIBCLANG_PATH` to the LLVM library directory if
bindgen cannot find libclang. Rust is installed separately with rustup.

```sh
python -m pip install 'maturin>=1.13,<2.0'
maturin build --release --interpreter python --out dist
python -m pip install dist/anarci_rs-*.whl
```

The HMMER build downloads the official 3.4 tarball and verifies SHA-256
`ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3`.
For an offline HMMER build, set `HMMER_TARBALL` to a local verified copy. Cargo and
Python build dependencies also need a populated cache or an available package index.

```sh
HMMER_TARBALL=/path/to/hmmer-3.4.tar.gz maturin build --release --interpreter python
```

Build on the target platform: native glibc Linux and macOS on x86_64 or ARM64.
Cross compilation, Windows, musl Linux and free-threaded Python are outside the
configured release matrix. HMMER uses position-independent C code and an x86-64
baseline on Intel/AMD systems. The workflow tests each platform before it qualifies
for a release; a workflow definition alone does not establish platform support.

## License and release status

Source code is BSD-3-Clause. Reference data carries additional attribution and
CC BY 4.0 terms. See `NOTICE.md` and `licenses/` for retained notices.
`docs/RELEASING.md` lists the artifact checks required before upload.

## Optional model research tools

The default wheel includes the experimental codec and data utilities. Build with
`maturin build --release --features learned-model --interpreter python` to enable
the optional CPU model research bindings. That build adds the Candle dependencies.
Model construction and training helpers require explicit user code and model data;
the package does not supply a validated trained numbering model. See
`docs/EXPERIMENTAL.md` for the experimental API and its limits.
