"""Console entry point for the native command-line interface."""

import sys


def main():
    from . import _native

    try:
        return _native.cli(sys.argv[1:])
    except (ValueError, OSError) as error:
        print(f"anarci-rs: {error}", file=sys.stderr)
        return 1
