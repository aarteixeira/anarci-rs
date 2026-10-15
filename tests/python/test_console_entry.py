"""Check that the installed console shim passes arguments and returns failures."""

import contextlib
import io
import sys
import unittest
from unittest.mock import patch

from anarci_rs import _native
from anarci_rs.cli import main


class ConsoleEntryTests(unittest.TestCase):
    def test_passes_arguments_without_program_name(self):
        with patch.object(sys, "argv", ["anarci-rs", "--help"]):
            with patch.object(_native, "cli", return_value=0) as cli:
                self.assertEqual(main(), 0)
                cli.assert_called_once_with(["--help"])

    def test_reports_native_argument_error(self):
        errors = io.StringIO()
        with patch.object(sys, "argv", ["anarci-rs", "--invalid"]):
            with patch.object(_native, "cli", side_effect=ValueError("invalid option")):
                with contextlib.redirect_stderr(errors):
                    self.assertEqual(main(), 1)
        self.assertEqual(errors.getvalue(), "anarci-rs: invalid option\n")


if __name__ == "__main__":
    unittest.main()
