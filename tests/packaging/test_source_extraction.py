"""Portable source extraction accepts regular data and rejects unsafe members."""
import io
from pathlib import Path
import runpy
import tarfile
import tempfile
import unittest

TOOLS = runpy.run_path(str(Path(__file__).resolve().parents[2] / "scripts/validate_release.py"))
REQUIRED = [
    "pyproject.toml", "Cargo.toml", "Cargo.lock",
    "reference_data/dat/HMMs/ALL.hmm", "reference_data/dat/HMMs/FEW.hmm",
    "python/anarci/__init__.py", "python/anarci_rs/__init__.py", "LICENSE", "NOTICE.md",
    "licenses/ANARCI-BSD-3-Clause.txt", "licenses/HMMER-LICENSE.txt",
    "licenses/Easel-LICENSE.txt", "licenses/CC-BY-4.0.txt",
    "licenses/IMGT-DATA-TERMS.txt", "licenses/RUST-THIRD-PARTY.txt",
]
ROOT = "anarci_rs-0.1.0"
PAYLOAD = b"synthetic packaging test data\n"


def source_archive(path, extra=None):
    with tarfile.open(path, "w:gz") as archive:
        for name in REQUIRED:
            member = tarfile.TarInfo(ROOT + "/" + name)
            member.size = len(PAYLOAD)
            member.mode = 0o4755
            archive.addfile(member, io.BytesIO(PAYLOAD))
        if extra is not None:
            archive.addfile(extra)


class SourceExtractionTests(unittest.TestCase):
    def test_extracts_data_without_archive_permission_bits(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            source = base / "source.tar.gz"
            source_archive(source)
            root = TOOLS["extract_source_archive"](source, base / "extracted")
            self.assertEqual(root, base.resolve() / "extracted" / ROOT)
            for name in REQUIRED:
                output = root / name
                self.assertEqual(output.read_bytes(), PAYLOAD)
                self.assertEqual(output.stat().st_mode & 0o7111, 0)

    def test_rejects_path_traversal_and_links_before_writing(self):
        for kind in ("traversal", "symlink", "hardlink"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as temporary:
                base = Path(temporary)
                source = base / "source.tar.gz"
                extra = tarfile.TarInfo(ROOT + "/../escape" if kind == "traversal" else ROOT + "/python/anarci_rs/__init__.py")
                if kind != "traversal":
                    extra.type = tarfile.SYMTYPE if kind == "symlink" else tarfile.LNKTYPE
                    extra.linkname = "../../escape"
                source_archive(source, extra)
                with self.assertRaises(ValueError):
                    TOOLS["extract_source_archive"](source, base / "extracted")
                self.assertFalse((base / "extracted").exists())
                self.assertFalse((base / "escape").exists())

    def test_preserves_existing_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            source = base / "source.tar.gz"
            source_archive(source)
            destination = base / "extracted"
            destination.mkdir()
            sentinel = destination / "existing.txt"
            sentinel.write_text("keep")
            with self.assertRaises(FileExistsError):
                TOOLS["extract_source_archive"](source, destination)
            self.assertEqual(sentinel.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
