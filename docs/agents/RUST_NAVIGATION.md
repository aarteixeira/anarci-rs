# Rust source navigation

Use the standalone helper to find Rust declarations, identifier occurrences, impl blocks and declared public items without running Cargo:

```sh
python3.12 scripts/rust_symbols.py definition HmmEngine
python3.12 scripts/rust_symbols.py references HmmEngine --limit 20
python3.12 scripts/rust_symbols.py implementations HmmEngine
python3.12 scripts/rust_symbols.py search scan --path crates/anarci-hmm/
python3.12 scripts/rust_symbols.py api --path crates/anarci-core/ --limit 20
python3.12 scripts/rust_symbols.py api HmmEngine --include-restricted
```

Python 3.11 or newer and Git are the only requirements. Run these commands from the repository root, or pass `--root /path/to/repository`. The helper does not import the native Python extension, run Cargo, install tools, download dependencies or regenerate repository maps. Use it alongside [agent_tools.py](../../scripts/agent_tools.py): component context narrows the search, then symbol navigation locates declarations and callers.

## Choose a query

| Command | Lexical behavior |
| --- | --- |
| `search QUERY` | Case-sensitive substring search over declaration names and scope hints. Includes impl headers. |
| `definition QUERY` | Exact declaration name or qualified scope hint. Returns all matches, not a resolved symbol. |
| `references QUERY` | Exact identifier occurrences, excluding indexed declaration-name tokens, comments and literals. A qualified query uses its final component. |
| `implementations QUERY` | Substring search over trait and target text in impl headers, ignoring whitespace. Includes inherent impls when their target matches. |
| `api [QUERY]` | Declarations with `pub` visibility and trait members, optionally filtered by a substring of the scope hint. Public use statements appear as `reexport` records. |

All searches are case-sensitive. `--path PREFIX` filters by repository-relative path prefix. For example, use `--path crates/anarci-core/`, including the trailing slash when you mean a directory. `--include-restricted` adds `pub(crate)`, `pub(super)` and `pub(in ...)` declarations to `api`; it does not change other commands.

A declaration record contains `kind`, `name`, `qualified_name`, `visibility`, a signature excerpt, `path`, `line` and `column`. Impl records also contain `trait` (null for inherent impls) and `target`. Method records identify their `trait` or `impl` container. Restricted visibility is retained as source tokens; spaces inside `pub(in ...)` are not preserved. Excerpts stop at 240 characters.

Inline modules and trait/impl methods keep their containing scope. Out-of-line module scopes are inferred from source filenames, with `src`, `lib.rs`, `main.rs` and `mod.rs` treated conventionally. These hints are not globally unique: different crates can have the same scope hint. Use the returned path to distinguish them.

## Read the JSON honestly

Successful lexical commands exit 0 and emit one JSON object. The envelope includes:

- `schema_version`, `source_digest` and `files_indexed`;
- `mode: "lexical"` and `semantics: "unresolved_candidates"`;
- `command`, `query`, `total`, `offset`, `limit` and `truncated`;
- `results` and explicit `limitations`.

`--limit` defaults to 30 and accepts 1 through 200. `--offset` accepts 0 through 1000000. Results sort by path, line, column and name; repeat the same query with a larger offset to page through them. `total` counts every match before pagination. Invalid arguments, unreadable sources and Git errors exit 2 with a JSON error object. `--help` prints ordinary CLI help.

Locations use one-based lines and Unicode code-point columns, not LSP UTF-16 positions or byte offsets. Reference records contain only the identifier and its location. Read the surrounding source before interpreting a match.

**A lexical occurrence is not a resolved reference.** Identically named types, local variables, imports and methods can all match. Even `other_crate::HmmEngine` searches occurrences of `HmmEngine`, without checking the qualifier. Implementation results likewise match written headers; they do not resolve type aliases, trait identities or blanket-implementation applicability.

`api` is a declaration inventory, not the crate's effective exported API. It includes public declarations inside private modules and members of private traits. It retains public reexports as written, without following them. It does not list struct fields, enum variants, tuple fields, generated declarations or declarations inside foreign `extern` blocks. Review module reachability, reexports and conditional compilation yourself before making an API compatibility claim.

## Source coverage and freshness

Each query obtains Rust paths from `git ls-files --cached --others --exclude-standard`. It includes tracked sources and nonignored, untracked candidate sources. It skips deleted files, symlinks, symlinked parent directories and paths containing `target`, `build`, `dist`, `node_modules`, `.git` or `.venv`. It never recursively scans ignored local research directories or changes their contents.

There is no persistent index. Every query rebuilds an in-memory source index and hashes the schema version, sorted paths and source bytes with SHA-256. Content edits, additions and deletions change `source_digest`; no timestamp-only freshness assumption is used. Concurrent source edits are not an atomic snapshot, so rerun after editing. No generated index needs committing or cleanup.

The lexer skips nested block comments, line comments, ordinary strings, byte/C strings, raw strings with arbitrary hash delimiters and character literals. It distinguishes lifetimes from character literals and normalizes raw identifiers such as `r#type`. The parser recognizes module nesting, generic impl headers, where clauses, trait members, function qualifiers and declaration visibility. It skips function bodies and macro bodies when finding declarations; identifier occurrences in those bodies can still appear as reference candidates.

This is a tolerant source scanner, not a Rust parser or compiler. It does not evaluate `cfg`, expand macros, resolve imports, follow `include!` or `#[path]`, derive crate names, or compute effective visibility. Unsupported or malformed syntax can yield missing or misleading candidates. In particular, local declarations inside function bodies are omitted. These limits apply even when a query returns no results.

## Resolve symbols with rust-analyzer

An optional [LSP adapter](../../scripts/rust_lsp.py) resolves definitions, references and implementations by source position:

```sh
python3.12 scripts/rust_symbols.py definition --mode semantic \
  --at crates/anarci-cli/src/engine.rs:13:6
python3.12 scripts/rust_symbols.py references --mode semantic \
  --at crates/anarci-core/src/orchestrate.rs:44:11 --limit 20
python3.12 scripts/rust_symbols.py implementations --mode semantic \
  --at crates/anarci-core/src/orchestrate.rs:44:11
```

Use an installed `rust-analyzer`, or pass `--analyzer /path/to/executable`. A rustup shim without the component is not an installed analyzer. The helper probes its version, then starts one bounded LSP session. It never installs a tool. `--timeout` defaults to 30 seconds; startup and project loading count toward that deadline.

Semantic mode takes `--at FILE:LINE:COLUMN`, using one-based codepoint positions, rather than a positional symbol name. It reports `semantics: "rust-analyzer"` and `resolution_status: "resolved"` or `"no_results"` only after a loading-to-quiescent server transition. Missing tools, incomplete loading, protocol errors and deadlines produce explicit errors with exit 3. There is no silent lexical fallback.

The adapter sets Cargo offline mode and disables check-on-save, build-script execution, procedural macros and sysroot discovery. Cached dependency metadata is still needed, and standard-library, generated and macro-dependent symbols may be incomplete. External source files, ignored files and symlinked files are not read to convert result positions; such locations retain explicitly labeled UTF-16 ranges. Results include raw LSP ranges and, where source can be read safely, one-based codepoint ranges.

Definitions, references and trait implementations were exercised against a real analyzer during implementation; the protocol also has isolated fake-server tests. Neither that check nor an empty result proves complete workspace coverage. Recheck returned source locations when APIs or conditional features change.

## Offline tests

```sh
python3.12 -m unittest discover -s tests/agent_tools -p 'test_rust_*.py' -v
```

The [tests](../../tests/agent_tools/test_rust_symbols.py) use synthetic Rust and temporary local Git repositories, plus Git-visible repository source. They cover comments and literals, lifetimes, raw identifiers, nested modules, generic trait/inherent impls, qualifiers, visibility, reexports, pagination, inventory exclusions, digest freshness, JSON errors and analyzer availability/dispatch behavior. Repository checks find the `HmmEngine` trait, its written impls, identifier candidates and declared public API. They require no native build, network access or installed Rust component.
