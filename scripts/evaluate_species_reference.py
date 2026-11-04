#!/usr/bin/env python3
"""Fit and audit real numbered OAS cohorts; write local artifacts, not validation claims.

Run with an installed current anarci-rs package. Dataset preparation is separate;
this command never downloads data or changes split assignments.
"""
import argparse
import bisect
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def distribution(values):
    values = sorted(values)
    if not values:
        return {"n": 0, "minimum": None, "median": None, "maximum": None}
    return {"n": len(values), "minimum": values[0], "median": statistics.median(values), "maximum": values[-1]}


def rank_auc(positives, negatives):
    """Mann-Whitney AUC with tied scores; descriptive, not independent trials."""
    if not positives or not negatives:
        return None
    negatives = sorted(negatives)
    wins = math.fsum(bisect.bisect_left(negatives, value) + 0.5 * (bisect.bisect_right(negatives, value) - bisect.bisect_left(negatives, value)) for value in positives)
    return wins / (len(positives) * len(negatives))


def raw_score(model, numbering):
    if not all(position in numbering and numbering[position] in model["probabilities"][position] for position in model["positions"]):
        return None
    return math.fsum(math.log(model["probabilities"][pos][numbering[pos]]) for pos in model["positions"]) / len(model["positions"])


def validate_cohort_origins(records, manifest):
    """Require explicit reviewed origin declarations, not inferred OAS species labels."""
    source_units = {unit["url"]: unit for unit in manifest["units"]}
    if len(source_units) != len(manifest["units"]):
        raise ValueError("Duplicate source units in manifest")
    for unit in source_units.values():
        if unit.get("biological_origin") != "natural_repertoire" or unit.get("library_selection") != "unselected":
            raise ValueError("Reference units require reviewed natural, unselected repertoire origin; OAS membership alone is insufficient")
        evidence = unit.get("origin_evidence", {})
        if not isinstance(evidence, dict) or not evidence.get("url") or not evidence.get("review_note"):
            raise ValueError("Biological-origin evidence is required for every unit")
    for row in records:
        unit = source_units.get(row.get("source_url"))
        if unit is None or any(row.get(key) != unit.get(key) for key in ("study", "split", "species", "domain_format")):
            raise ValueError("Record source or cohort assignment differs from the reviewed retrieval manifest")
        if row.get("biological_origin") != "natural_repertoire" or row.get("library_selection") != "unselected":
            raise ValueError("Prepared record lacks reviewed natural-origin eligibility")


def evaluate(records_path, manifest_path, output_dir, atlas_id, version):
    from anarci_rs.species_reference import fit_profile, save_profile, audit_profile, FRAMEWORK_POSITIONS
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise ValueError("Output directory already exists; choose a new versioned artifact directory")
    records = []
    with Path(records_path).open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                row = json.loads(line)
                if row.get("provenance", {}).get("source") == "synthetic":
                    raise ValueError("This evaluator is for real source cohorts, not synthetic references")
                records.append(row)
                if len(records) > 200_000:
                    raise ValueError("Bounded baseline evaluator accepts at most200000 source records")
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or not manifest.get("units"):
        raise ValueError("A source-unit manifest is required")
    validate_cohort_origins(records, manifest)
    groups = defaultdict(list)
    for row in records:
        groups[(row["species"], row["chain"], row["domain_format"])].append(row)
    if not groups:
        raise ValueError("No prepared source records")
    # Matched coordinate panels are selected from training records only, never
    # calibration or test records. This allows descriptive same-format comparisons.
    panels = {}
    for (_, chain, domain_format), rows in groups.items():
        key = (chain, domain_format)
        training = [row for row in rows if row["split"] == "train"]
        if not training:
            raise ValueError("Each profile requires training records")
        counts = Counter(pos for row in training for pos, aa in row["numbering"].items()
                         if pos in FRAMEWORK_POSITIONS and aa in "ACDEFGHIKLMNPQRSTVWY")
        # A single deletion outlier must not erase a position for the whole atlas.
        # Freeze this prevalence rule before calibration/test evaluation; reject
        # incomplete training rows rather than imputing amino acids or weakening
        # the fitter's complete-coordinate contract.
        available = {pos for pos, count in counts.items() if count / len(training) >= 0.99}
        panels[key] = available if key not in panels else panels[key] & available
    if any(len(panel) < 60 for panel in panels.values()):
        raise ValueError("Fewer than60 common training framework positions; review cohort coverage")
    models = {}
    fit_exclusions = {}
    for (species, chain, domain_format), rows in sorted(groups.items()):
        units = [unit for unit in manifest["units"] if unit["study"] in {row["study"] for row in rows}]
        provenance = {"source": "Observed Antibody Space / OAS-nano", "units": units,
                      "retrieval_manifest_sha256": sha256(manifest_path), "prepared_records_sha256": sha256(records_path),
                      "note": "Local downloaded data checksums and recorded metadata, not upstream certification"}
        positions = sorted(panels[(chain, domain_format)], key=int)
        fit_rows = [row for row in rows if row["split"] != "train" or all(pos in row["numbering"] for pos in positions)]
        excluded = len(rows) - len(fit_rows)
        fit_exclusions[(species, chain, domain_format)] = excluded
        provenance["coordinate_policy"] = {"training_position_min_prevalence": 0.99,
                                            "shared_within_chain_format": True,
                                            "incomplete_training_rows_excluded": excluded,
                                            "calibration_test_used_for_position_selection": False}
        model = fit_profile(fit_rows, species=species, chain=chain, domain_format=domain_format,
                            atlas_id=atlas_id, version=version, positions=positions, source_manifest=provenance)
        models[(species, chain, domain_format)] = model
    output_dir.mkdir(parents=True)
    report = {"schema_version": 1, "atlas_id": atlas_id, "version": version,
              "prepared_records_sha256": sha256(records_path), "retrieval_manifest_sha256": sha256(manifest_path),
              "source_record_count": len(records), "profiles": [], "pairwise_test_discrimination": [],
              "limitations": [
                  "Bounded unit selection is a pilot, not a representative species repertoire census.",
                  "Study-held-out sequence ranks do not establish immunogenicity, adaptation success or biological fitness.",
                  "Exact-sequence and CDR3/V-family proxy overlap filters do not prove absence of all clonal or ancestral relatedness.",
                  "Unknown subject identities are unit proxies, not verified independent donors.",
                  "AUC is descriptive on the selected test records; repeated sequences within donors are not independent experiments.",
                  "VH and VHH format models are not interchangeable; absent species/chain/format combinations are unsupported.",
                  "No experimental antibody adaptation or independent structure evaluation has been performed."]}
    for key, model in models.items():
        species, chain, domain_format = key
        filename = "%s-%s-%s.json" % key
        save_profile(model, output_dir / filename)
        per_split = {}
        for split in ("train", "calibration", "test"):
            rows = [row for row in groups[key] if row["split"] == split]
            scores = [raw_score(model, row["numbering"]) for row in rows]
            valid = [value for value in scores if value is not None]
            calibration = model["calibration_scores"]
            percentiles = [100.0 * (bisect.bisect_left(calibration, value) + 0.5 * (bisect.bisect_right(calibration, value) - bisect.bisect_left(calibration, value))) / len(calibration) for value in valid] if calibration else []
            per_split[split] = {"records": len(rows), "studies": sorted({row["study"] for row in rows}),
                                "full_coordinate_coverage": len(valid), "raw_score": distribution(valid),
                                "reference_percentile": distribution(percentiles)}
        report["profiles"].append({"file": filename, "species": species, "chain": chain,
                                   "domain_format": domain_format, "audit": audit_profile(model),
                                   "incomplete_training_rows_excluded": fit_exclusions[key], "split_results": per_split})
        positive = [raw_score(model, row["numbering"]) for row in groups[key] if row["split"] == "test"]
        for other_key, other_rows in groups.items():
            if other_key == key or other_key[1:] != key[1:]:
                continue
            negative = [raw_score(model, row["numbering"]) for row in other_rows if row["split"] == "test"]
            positive_valid = [value for value in positive if value is not None]
            negative_valid = [value for value in negative if value is not None]
            report["pairwise_test_discrimination"].append({"model_species": species, "other_species": other_key[0],
                                                          "chain": chain, "domain_format": domain_format,
                                                          "positive_records": len(positive_valid), "negative_records": len(negative_valid),
                                                          "auc": rank_auc(positive_valid, negative_valid)})
    with (output_dir / "evaluation.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--atlas-id", required=True)
    parser.add_argument("--version", required=True)
    args = parser.parse_args(argv)
    try:
        report = evaluate(args.records, args.manifest, args.output_dir, args.atlas_id, args.version)
        print(json.dumps({"source_record_count": report["source_record_count"],
                          "profiles": [{"species": row["species"], "format": row["domain_format"], "evidence": row["audit"]["evidence"]} for row in report["profiles"]],
                          "output_dir": str(args.output_dir)}, indent=2))
        return 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print("evaluate_species_reference: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
