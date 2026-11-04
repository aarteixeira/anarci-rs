"""Explicit-data CLI for experimental species-reference scoring and adaptation."""
import argparse
import json
from pathlib import Path
import sys

from .species_reference import audit_profile, fit_profile, load_profile, save_profile
from .species_adaptation import profile_sequence, propose_sequence


def _records(path):
    source = Path(path)
    if source.stat().st_size > 256 * 1024 * 1024:
        raise ValueError("Record input exceeds 256 MiB; use a documented bounded cohort")
    rows = []
    with source.open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                rows.append(json.loads(line))
                if len(rows) > 200_000:
                    raise ValueError("At most 200000 records may be fitted in one baseline profile")
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit", help="Validate a supplied model and inspect evidence/provenance")
    audit.add_argument("model", type=Path)
    audit.add_argument("--output", type=Path)
    fit = sub.add_parser("fit", help="Fit one baseline profile from supplied numbered JSONL records")
    fit.add_argument("records", type=Path)
    fit.add_argument("--species", required=True)
    fit.add_argument("--chain", required=True, choices=["H", "K", "L"])
    fit.add_argument("--domain-format", required=True, choices=["VH", "VL", "VHH"])
    fit.add_argument("--atlas-id", required=True)
    fit.add_argument("--version", required=True)
    fit.add_argument("--source-manifest", type=Path, required=True)
    fit.add_argument("--pseudocount", type=float, default=0.5)
    fit.add_argument("--output", type=Path, required=True)
    for name in ("score", "propose"):
        command = sub.add_parser(name)
        command.add_argument("sequence", help="Canonical amino-acid sequence (not a file)")
        command.add_argument("--domain-format", required=True, choices=["VH", "VL", "VHH"])
        command.add_argument("--domain-index", type=int, default=0)
        command.add_argument("--database", choices=["ALL", "pan"], default="ALL")
        command.add_argument("--output", type=Path)
        if name == "score":
            command.add_argument("--model", type=Path, action="append", required=True)
            command.add_argument("--species", action="append")
        else:
            command.add_argument("--target-model", type=Path, required=True)
            command.add_argument("--guard-model", type=Path, action="append", default=[])
            command.add_argument("--protect", nargs="*", default=[], help="Additional immutable IMGT positions")
            command.add_argument("--mutable", nargs="*", help="Optional explicit allowed IMGT positions")
            command.add_argument("--max-mutations", type=int, default=3)
            command.add_argument("--max-candidates", type=int, default=8)
            command.add_argument("--beam-width", type=int, default=16)
            command.add_argument("--max-evaluations", type=int, default=5000)
    args = parser.parse_args(argv)
    try:
        if args.output is not None and args.output.exists():
            raise ValueError("Output already exists; select a new path")
        if args.command == "fit":
            manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
            model = fit_profile(_records(args.records), species=args.species, chain=args.chain,
                                domain_format=args.domain_format, atlas_id=args.atlas_id, version=args.version,
                                pseudocount=args.pseudocount, source_manifest=manifest)
            save_profile(model, args.output)
            print(json.dumps(audit_profile(model), indent=2, sort_keys=True, allow_nan=False))
            return 0
        if args.command == "audit":
            result = audit_profile(load_profile(args.model))
        elif args.command == "score":
            result = profile_sequence(args.sequence, [load_profile(path) for path in args.model],
                                      domain_format=args.domain_format, domain_index=args.domain_index,
                                      database=args.database, species=args.species)
        else:
            result = propose_sequence(args.sequence, load_profile(args.target_model), domain_format=args.domain_format,
                                      domain_index=args.domain_index, database=args.database,
                                      protected_positions=args.protect, mutable_positions=args.mutable,
                                      guard_profiles=[load_profile(path) for path in args.guard_model],
                                      max_mutations=args.max_mutations, max_candidates=args.max_candidates,
                                      beam_width=args.beam_width, max_evaluations=args.max_evaluations)
        text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.output is not None:
            with args.output.open("x", encoding="utf-8") as stream:
                stream.write(text)
        else:
            print(text, end="")
        return 0
    except (OSError, ValueError, AssertionError, KeyError, TypeError) as exc:
        print("anarci-species: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
