"""Synthetic and repository contracts; stdlib only, no Rust/native builds."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("rust_symbols", ROOT / "scripts/rust_symbols.py")
rust = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rust)


class LexerTests(unittest.TestCase):
    def test_comments_literals_raw_identifiers_and_lifetimes(self):
        source = '''// fn Ghost() {}
        /* outer /* nested */ struct Ghost; */
        const TEXT: &str = "fn Ghost() { \\\" }";
        const RAW: &str = r###"struct Ghost; r#"nested"#"###;
        const BYTE: &[u8] = br##"Ghost"##;
        const C: &CStr = cr#"Ghost"#;
        const B: u8 = b'G';
        const U: char = '\\u{0047}';
        const ESC: char = '\\'';
        fn r#type<'a>(x: &'a str) { let s = "Ghost"; }
        '''
        values = [t["value"] for t in rust.tokenize(source)]
        self.assertNotIn("Ghost", values)
        self.assertNotIn("nested", values)
        self.assertIn("a", values)
        defs, refs = rust.parse_source("src/lib.rs", source)
        self.assertEqual(["TEXT", "RAW", "BYTE", "C", "B", "U", "ESC", "type"], [d["name"] for d in defs])
        self.assertNotIn("Ghost", [r["name"] for r in refs])

    def test_unterminated_comment_or_string_does_not_leak(self):
        for tail in ['/* fn Hidden() {}', 'r##"fn Hidden() {}', '"fn Hidden() {}']:
            defs, _ = rust.parse_source("src/lib.rs", "fn visible() {} " + tail)
            self.assertEqual(["visible"], [d["name"] for d in defs])

    def test_unicode_location_uses_codepoints(self):
        defs, refs = rust.parse_source("src/lib.rs", "// é\npub fn café() { café(); }\n")
        self.assertEqual((2, 8), (defs[0]["line"], defs[0]["column"]))
        self.assertEqual(17, next(r["column"] for r in refs if r["name"] == "café"))


class ParserTests(unittest.TestCase):
    SOURCE = '''
pub mod outer {
    pub(crate) struct Thing<T> { pub field: T }
    pub trait Service<T>: Send {
        type Output;
        const VALUE: usize;
        fn run<'a>(&'a self, arg: T) -> Self::Output;
        fn defaulted(&self) { helper(); }
    }
    impl<T: Clone> Service<T> for Thing<T> where T: Send {
        type Output = T;
        const VALUE: usize = { 1 };
        fn run<'a>(&'a self, arg: T) -> T { arg }
    }
    impl<T> Thing<T> {
        #[must_use]
        pub const unsafe fn new(value: T) -> Self { Self { field: value } }
        pub(super) fn restricted(&self) {}
        fn hidden(&self) { fn local() {} }
    }
    pub use self::Thing as Exported;
    pub use elsewhere::{One, Two};
    mod inner { pub fn same() {} }
}
pub fn same() {}
macro_rules! generated { () => { pub fn fake() {} } }
generated! { struct AlsoFake; }
'''

    def index(self):
        defs, refs = rust.parse_source("crates/demo/src/lib.rs", self.SOURCE)
        return {"definitions": defs, "occurrences": refs, "files": ["crates/demo/src/lib.rs"], "source_digest": "fixture"}

    def test_nested_modules_and_method_containers(self):
        index = self.index()
        same = rust.query_index(index, "definition", "same")["results"]
        self.assertEqual(["outer::inner::same", "same"], [d["qualified_name"] for d in same])
        self.assertEqual(1, rust.query_index(index, "definition", "outer::inner::same")["total"])
        names = [d["name"] for d in index["definitions"]]
        self.assertNotIn("local", names)
        self.assertNotIn("fake", names)
        self.assertNotIn("AlsoFake", names)
        run = rust.query_index(index, "definition", "run")["results"]
        self.assertEqual(["trait", "impl"], [d["container"] for d in run])

    def test_impl_generics_where_clauses_and_trait(self):
        rows = rust.query_index(self.index(), "implementations", "Service")["results"]
        self.assertEqual(1, len(rows))
        self.assertEqual("Service < T >", rows[0]["trait"])
        self.assertEqual("Thing < T >", rows[0]["target"])
        self.assertEqual(2, rust.query_index(self.index(), "implementations", "Thing")["total"])

    def test_const_generics_and_hrtb(self):
        source = '''struct Buffer<const N: usize>;
        impl<const N: usize> Buffer<{ N }> { pub fn new() -> Self { todo!() } }
        impl<T> Visit<for<'a> fn(&'a T)> for Buffer<2> { fn visit(&self) {} }
        '''
        defs, _ = rust.parse_source("src/mod.rs", source)
        impls = [d for d in defs if d["kind"] == "impl"]
        self.assertEqual("Buffer < { N } >", impls[0]["target"])
        self.assertEqual("Buffer < 2 >", impls[1]["target"])
        self.assertIn("Visit", impls[1]["trait"])
        self.assertEqual(["new", "visit"], [d["name"] for d in defs if d["kind"] == "fn"])

    def test_qualified_generic_impl_queries_and_raw_identifier_queries(self):
        defs, refs = rust.parse_source("src/lib.rs", """
        mod r#outer {
            struct Thing<T>(T);
            impl<T> crate::traits::Work<T> for Thing<T> { fn r#type(&self) {} }
            fn use_it(x: Thing<u8>) { x.r#type(); }
        }
        """)
        index = {"definitions": defs, "occurrences": refs, "files": ["src/lib.rs"], "source_digest": "fixture"}
        self.assertEqual(1, rust.query_index(index, "implementations", "crate::traits::Work<T>")["total"])
        self.assertEqual(1, rust.query_index(index, "references", "outer::r#type")["total"])
        self.assertEqual(1, rust.query_index(index, "definition", "r#outer::use_it")["total"])

    def test_visibility_and_api_are_declarations_not_exports(self):
        result = rust.query_index(self.index(), "api")
        self.assertEqual("unresolved_candidates", result["semantics"])
        names = [d["name"] for d in result["results"]]
        self.assertIn("new", names)
        self.assertIn("run", names)
        self.assertIn("Output", names)
        self.assertNotIn("Thing", names)
        self.assertNotIn("hidden", names)
        self.assertNotIn("restricted", names)
        restricted = rust.query_index(self.index(), "api", include_restricted=True)
        self.assertIn("Thing", [d["name"] for d in restricted["results"]])
        self.assertIn("restricted", [d["name"] for d in restricted["results"]])
        self.assertEqual(2, sum(d["kind"] == "reexport" for d in result["results"]))

    def test_references_are_exact_identifier_candidates_not_definitions(self):
        result = rust.query_index(self.index(), "references", "Thing")
        self.assertEqual(3, result["total"])
        self.assertTrue(all(row["name"] == "Thing" for row in result["results"]))
        self.assertEqual("lexical", result["mode"])
        self.assertIn("no name/type resolution", result["limitations"][0])
        self.assertEqual(result["total"], rust.query_index(self.index(), "references", "other::Thing")["total"])

    def test_bounded_deterministic_pagination_and_prefix_filter(self):
        index = self.index()
        full = rust.query_index(index, "search", "same")
        first = rust.query_index(index, "search", "same", limit=1)
        second = rust.query_index(index, "search", "same", limit=1, offset=1)
        self.assertTrue(first["truncated"])
        self.assertFalse(second["truncated"])
        self.assertEqual(full["results"], first["results"] + second["results"])
        self.assertEqual(0, rust.query_index(index, "search", "same", path="different/")["total"])
        for limit, offset in [(0, 0), (201, 0), (1, -1)]:
            with self.assertRaises(ValueError):
                rust.query_index(index, "search", "same", limit=limit, offset=offset)

    def test_out_of_line_module_hints(self):
        defs, _ = rust.parse_source("crates/demo/src/foo/mod.rs", "pub mod bar; pub fn item() {}")
        self.assertEqual(["foo::bar", "foo::item"], [d["qualified_name"] for d in defs])


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        (self.root / ".gitignore").write_text("research/\ntarget/\n", encoding="utf-8")
        (self.root / "src").mkdir()
        (self.root / "src/lib.rs").write_text("pub trait Visible {}", encoding="utf-8")
        subprocess.run(["git", "add", "src/lib.rs", ".gitignore"], cwd=self.root, check=True)

    def test_git_visible_only_and_symlink_exclusion(self):
        (self.root / "src/new.rs").write_text("pub fn candidate() {}")
        for dirname in ["research", "target", "build"]:
            (self.root / dirname).mkdir()
            (self.root / dirname / "local.rs").write_text("struct Secret;")
        (self.root / "src/link.rs").symlink_to(self.root / "research/local.rs")
        self.assertEqual(["src/lib.rs", "src/new.rs"], rust.inventory(self.root))
        self.assertEqual("struct Secret;", (self.root / "research/local.rs").read_text())

    def test_digest_freshness_content_inventory_and_schema(self):
        before = rust.build_index(self.root)
        self.assertEqual(before, rust.build_index(self.root))
        (self.root / "src/lib.rs").write_text("pub trait Changed {}")
        after = rust.build_index(self.root)
        self.assertNotEqual(before["source_digest"], after["source_digest"])
        (self.root / "src/added.rs").write_text("pub fn added() {}")
        added = rust.build_index(self.root)
        self.assertNotEqual(after["source_digest"], added["source_digest"])
        with patch.object(rust, "SCHEMA", rust.SCHEMA + 1):
            self.assertNotEqual(added["source_digest"], rust.build_index(self.root)["source_digest"])
        (self.root / "src/lib.rs").unlink()
        self.assertNotIn("src/lib.rs", rust.inventory(self.root))

    def test_cli_json_and_errors(self):
        def invoke(*args):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                try:
                    code = rust.main([*args, "--root", str(self.root)])
                except SystemExit as exc:
                    code = exc.code
            return code, json.loads(out.getvalue())
        code, output = invoke("definition", "Visible", "--limit", "1")
        self.assertEqual(0, code)
        self.assertEqual("Visible", output["results"][0]["name"])
        for args in [("search",), ("search", "x", "--limit", "0"), ("search", "x", "--unknown"),
                     ("definition", "--at", "src/lib.rs:1:1"), ("definition", "--mode", "semantic")]:
            code, output = invoke(*args)
            self.assertEqual(2, code)
            self.assertEqual("error", output["status"])


class SemanticTests(unittest.TestCase):
    def test_missing_analyzer_is_explicit_and_never_lexical(self):
        with patch.object(rust.shutil, "which", return_value=None), patch.object(rust.subprocess, "run") as run:
            result = rust.semantic_status()
        run.assert_not_called()
        self.assertEqual("unavailable", result["status"])
        self.assertEqual("not_queried", result["semantics"])
        self.assertEqual([], result["results"])

    def test_shim_failure_and_offline_probe(self):
        failure = subprocess.CompletedProcess([], 1, "", "component unavailable")
        with patch.object(rust.shutil, "which", return_value="/bin/rust-analyzer"), patch.object(rust.subprocess, "run", return_value=failure) as run:
            result = rust.semantic_status()
        self.assertEqual("component unavailable", result["reason"])
        self.assertEqual(["/bin/rust-analyzer", "--version"], run.call_args.args[0])
        self.assertEqual("true", run.call_args.kwargs["env"]["CARGO_NET_OFFLINE"])
        self.assertEqual("0", run.call_args.kwargs["env"]["RUSTUP_AUTO_INSTALL"])
        self.assertEqual("1", run.call_args.kwargs["env"]["RUSTUP_NO_UPDATE_CHECK"])
        self.assertEqual(5, run.call_args.kwargs["timeout"])

    def test_installed_analyzer_delegates_to_real_semantic_adapter(self):
        version = subprocess.CompletedProcess([], 0, "rust-analyzer test", "")
        response = {"status": "ok", "mode": "semantic", "semantics": "rust-analyzer", "results": [{"path": "src/lib.rs", "line": 1}]}
        with patch.object(rust.shutil, "which", return_value="/bin/rust-analyzer"), patch.object(rust.subprocess, "run", return_value=version), patch.object(rust, "semantic_query", return_value=response) as adapter:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = rust.main(["references", "--mode", "semantic", "--at", "src/lib.rs:1:2"])
        self.assertEqual(0, code)
        result = json.loads(out.getvalue())
        self.assertEqual("rust-analyzer", result["semantics"])
        self.assertEqual("src/lib.rs:1:2", adapter.call_args.args[2])
        self.assertEqual(1, result["total"])

    def test_probe_timeout(self):
        with patch.object(rust.shutil, "which", return_value="ra"), patch.object(rust.subprocess, "run", side_effect=subprocess.TimeoutExpired("ra", 5)):
            self.assertEqual("unavailable", rust.semantic_status()["status"])


class RepositoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = rust.build_index(ROOT)

    def test_hmm_engine_definition_references_and_implementations(self):
        definition = rust.query_index(self.index, "definition", "HmmEngine")
        self.assertTrue(any(row["kind"] == "trait" and row["path"] == "crates/anarci-core/src/orchestrate.rs" for row in definition["results"]))
        impls = rust.query_index(self.index, "implementations", "HmmEngine")
        self.assertGreaterEqual(impls["total"], 6)
        self.assertTrue(any(row["path"] == "crates/anarci-cli/src/engine.rs" for row in impls["results"]))
        self.assertGreater(rust.query_index(self.index, "references", "HmmEngine")["total"], 6)

    def test_repository_public_api_and_bounds(self):
        result = rust.query_index(self.index, "api", "HmmEngine", path="crates/anarci-core/")
        self.assertTrue(any(row["name"] == "HmmEngine" for row in result["results"]))
        result = rust.query_index(self.index, "api", limit=2)
        self.assertEqual(2, len(result["results"]))
        self.assertTrue(result["truncated"])
        self.assertFalse(any("target" in Path(path).parts for path in self.index["files"]))


if __name__ == "__main__":
    unittest.main()
