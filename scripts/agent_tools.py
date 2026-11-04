#!/usr/bin/env python3
"""Offline repository maps and conservative change routing (Python 3.11+)."""
from __future__ import annotations

import argparse
import ast
import fnmatch
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

if sys.version_info < (3, 11):
    raise SystemExit("agent_tools requires Python 3.11+ (development tooling only). Try python3.11 or python3.12.")
import tomllib

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ("docs/agents/repo-map.json", "docs/agents/ARCHITECTURE.md")
SKIP = {"target", "dist", ".git", ".venv", "venv", "__pycache__", ".ruff_cache", "node_modules", "outputs", "output", "reports"}
SUFFIXES = {".rs", ".py", ".toml", ".md", ".json", ".yml", ".yaml", ".c", ".h", ".sh", ".js"}


def load_policy(root: Path) -> dict:
    policy = json.loads((root / ".agents/repo.json").read_text())
    if policy["schema_version"] != 1:
        raise ValueError("Unsupported policy schema")
    return policy


def normalize_path(root: Path, value: str) -> str:
    path = Path(value)
    # Resolve lexically, including deleted files, without following source symlinks.
    absolute = Path(os.path.abspath(root / path))
    try:
        return absolute.relative_to(root.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError(f"Path outside repository: {value}") from exc


def components(root: Path, policy: dict) -> list[dict]:
    workspace = tomllib.loads((root / "Cargo.toml").read_text())
    result = []
    for member in workspace["workspace"]["members"]:
        matches = sorted(root.glob(member))
        if not matches:
            raise ValueError(f"Workspace member missing: {member}")
        for directory in matches:
            manifest = directory / "Cargo.toml"
            data = tomllib.loads(manifest.read_text())
            result.append({"id": data["package"]["name"], "path": directory.relative_to(root).as_posix(),
                           "description": data["package"].get("description", data["package"]["name"]),
                           "manifest": manifest.relative_to(root).as_posix(),
                           "features": data.get("features", {}),
                           "recipes": policy["crate_recipes"].get(data["package"]["name"], [])})
    result.extend(dict(component) for component in policy["components"])
    ids = [component["id"] for component in result]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate component ids")
    for component in result:
        if not (root / component["path"]).is_dir():
            raise ValueError(f"Missing component directory: {component['path']}")
        for recipe in component["recipes"]:
            if recipe not in policy["recipes"]:
                raise ValueError(f"Unknown recipe: {recipe}")
    return sorted(result, key=lambda item: item["id"])


def owners(path: str, comps: list[dict], policy: dict) -> list[str]:
    routed = {name for route in policy["routes"] if any(fnmatch.fnmatchcase(path, pattern) for pattern in route["patterns"]) for name in route["components"]}
    if routed:
        return sorted(routed)
    return sorted(c["id"] for c in comps if path == c["path"] or path.startswith(c["path"] + "/"))


def git_inventory(root: Path) -> set[str] | None:
    if not (root / ".git").exists():
        return None  # Synthetic repositories can exercise the generator without Git.
    result = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=root, check=True, stdout=subprocess.PIPE)
    return {os.fsdecode(path) for path in result.stdout.split(b"\0") if path}


def source_files(root: Path):
    # Explicit roots plus Git's ignore rules avoid local artifacts and native trees.
    eligible = git_inventory(root)
    for base in ("crates", "python", "scripts", "tests", "docs", ".agents", ".github", ".cargo"):
        directory = root / base
        if not directory.exists() or directory.is_symlink():
            continue
        for current, dirs, files in os.walk(directory, followlinks=False):
            dirs[:] = sorted(d for d in dirs if d not in SKIP and not (Path(current) / d).is_symlink())
            for name in sorted(files):
                file = Path(current) / name
                relative = file.relative_to(root).as_posix()
                if file.suffix in SUFFIXES and not file.is_symlink() and relative not in OUTPUTS and (eligible is None or relative in eligible):
                    yield relative
    root_names = {".gitignore", "Cargo.toml", "Cargo.lock", "pyproject.toml", "AGENTS.md", "README.md", "NOTICE.md", "LICENSE"}
    for file in sorted(root.iterdir()):
        if file.name in root_names and file.is_file() and not file.is_symlink() and (eligible is None or file.name in eligible):
            yield file.name


def build_graph(root: Path) -> dict:
    policy = load_policy(root)
    comps = components(root, policy)
    by_id = {c["id"]: c for c in comps}
    by_path = {c["path"]: c["id"] for c in comps}
    for route in policy["routes"]:
        for name in route["components"]:
            if name not in by_id:
                raise ValueError(f"Unknown route component: {name}")
    nodes = [{"id": c["id"], "kind": "component", **c} for c in comps]
    edges = []
    for component in comps:
        for dependency in component.get("depends_on", []):
            if dependency not in by_id:
                raise ValueError(f"Unknown dependency: {dependency}")
            edges.append({"from": component["id"], "to": dependency, "kind": "depends_on", "evidence": ".agents/repo.json"})
        if "manifest" not in component:
            continue
        manifest = tomllib.loads((root / component["manifest"]).read_text())
        sections = [(kind, manifest.get(kind, {})) for kind in ("dependencies", "dev-dependencies", "build-dependencies")]
        for target, sections_data in manifest.get("target", {}).items():
            sections.extend((f"target.{target}.{kind}", sections_data.get(kind, {})) for kind in ("dependencies", "dev-dependencies", "build-dependencies"))
        for section, dependencies in sections:
            for alias, specification in dependencies.items():
                if not isinstance(specification, dict):
                    continue
                if specification.get("workspace"):
                    workspace = tomllib.loads((root / "Cargo.toml").read_text())
                    specification = workspace.get("workspace", {}).get("dependencies", {}).get(alias, {})
                    base = root
                else:
                    base = root / component["path"]
                if not isinstance(specification, dict) or "path" not in specification:
                    continue
                dependency_path = normalize_path(root, str(base / specification["path"]))
                if dependency_path not in by_path:
                    raise ValueError(f"Unmapped local dependency: {dependency_path}")
                edges.append({"from": component["id"], "to": by_path[dependency_path], "kind": "depends_on", "scope": section,
                              "optional": specification.get("optional", False), "evidence": component["manifest"]})
    for dependency in policy.get("data_dependencies", []):
        if dependency["from"] not in by_id or dependency["to"] not in by_id or not (root / dependency["evidence"]).is_file():
            raise ValueError(f"Invalid data dependency: {dependency}")
        edges.append({**dependency, "kind": "depends_on", "scope": "reviewed-data"})
    files = sorted(set(source_files(root)))
    for file in files:
        nodes.append({"id": "file:" + file, "kind": "file", "path": file, "components": owners(file, comps, policy)})
        for owner in owners(file, comps, policy):
            edges.append({"from": owner, "to": "file:" + file, "kind": "contains", "evidence": file})
    modules = {}
    for file in files:
        if file.startswith("python/") and file.endswith(".py"):
            module = file[len("python/"):-3].replace("/", ".")
            if module.endswith(".__init__"):
                module = module[:-9]
            modules[module] = file
    for module, file in sorted(modules.items()):
        tree = ast.parse((root / file).read_text(), filename=file)
        package = module if file.endswith("/__init__.py") else module.rpartition(".")[0]
        for node in ast.walk(tree):
            targets = []
            if isinstance(node, ast.Import):
                targets = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                base = node.module or ""
                if node.level:
                    base = importlib.util.resolve_name("." * node.level + base, package)
                targets = [base] + [base + "." + alias.name for alias in node.names]
            for target in sorted(set(targets)):
                if target in modules and modules[target] != file:
                    edges.append({"from": "file:" + file, "to": "file:" + modules[target], "kind": "imports", "evidence": file, "line": node.lineno})
    for flow in policy["flows"]:
        if flow["from"] not in by_id or flow["to"] not in by_id or not (root / flow["evidence"]).is_file():
            raise ValueError(f"Invalid curated flow: {flow}")
        edges.append({**flow, "kind": "flow"})
    edges = sorted({json.dumps(edge, sort_keys=True): edge for edge in edges}.values(), key=lambda edge: (edge["kind"], edge["from"], edge["to"], json.dumps(edge, sort_keys=True)))
    return {"schema_version": 1, "notes": ["Cargo edges include all declared local optional, target, build and dev dependencies, not a resolved feature configuration. Reviewed data dependencies are separately labeled.", "Python import edges are static AST references to local source modules; dynamic imports and the native extension are not resolved.", "Flow edges are reviewed policy annotations, not extracted call graphs. File containment does not imply test coverage.", "Inventory covers Git-visible source/config/docs in explicit roots plus named root entry points, including untracked candidate source. Reference datasets, binary fixtures, ignored untracked files, symlinks and build artifacts are excluded; commit intended new inputs with generated maps."],
            "nodes": sorted(nodes, key=lambda node: node["id"]), "edges": edges,
            "recipes": policy["recipes"]}


def architecture(graph: dict) -> str:
    comps = [n for n in graph["nodes"] if n["kind"] == "component"]
    ids = {c["id"]: "n" + str(index) for index, c in enumerate(comps)}
    lines = ["# Repository architecture", "", "Generated by `python3 scripts/agent_tools.py generate`. Do not edit this file.", "", "## Declared dependencies", "", "Arrows point from consumer to dependency. Local Cargo dependency kinds, optional features and reviewed data dependencies are included.", "", "```mermaid", "flowchart TD"]
    for c in comps:
        lines.append(f'  {ids[c["id"]]}["{c["id"]}"]')
    pairs = sorted({(e["from"], e["to"]) for e in graph["edges"] if e["kind"] == "depends_on"})
    for source, target in pairs:
        lines.append(f"  {ids[source]} --> {ids[target]}")
    lines.extend(["```", "", "## Annotation flow", "", "Reviewed entry-point relationships; these arrows are not a complete execution trace.", "", "```mermaid", "flowchart LR"])
    for e in graph["edges"]:
        if e["kind"] == "flow":
            lines.append(f'  {ids[e["from"]]}["{e["from"]}"] -->|"{e["label"]}"| {ids[e["to"]]}["{e["to"]}"]')
    lines.extend(["```", "", "## Components", "", "| Component | Source | Purpose | Validation recipes |", "|---|---|---|---|"])
    for c in comps:
        lines.append(f'| {c["id"]} | [{c["path"]}](../../{c["path"]}) | {c["description"]} | {", ".join(c["recipes"]) or "Review task-specific checks"} |')
    lines.extend(["", "## Evidence and limits", ""])
    for note in graph["notes"]:
        lines.append("- " + note)
    lines.extend(["", "Flow evidence:", ""])
    for e in graph["edges"]:
        if e["kind"] == "flow":
            lines.append(f'- {e["from"]} to {e["to"]}: [{e["evidence"]}](../../{e["evidence"]})')
    lines.extend(["", "The [machine-readable map](repo-map.json) includes file ownership, local Python imports, Cargo features, dependency evidence and validation recipes. Use `impact` for change routing. Unmapped files and missing recipes require a task-specific review.", ""])
    return "\n".join(lines)


def rendered(root: Path) -> dict[str, str]:
    graph = build_graph(root)
    return {OUTPUTS[0]: json.dumps(graph, indent=2, sort_keys=True) + "\n", OUTPUTS[1]: architecture(graph)}


def impact(root: Path, paths: list[str]) -> dict:
    policy = load_policy(root)
    graph = build_graph(root)
    comps = [n for n in graph["nodes"] if n["kind"] == "component"]
    normalized = sorted({normalize_path(root, path) for path in paths})
    matched = {path: owners(path, comps, policy) for path in normalized}
    affected = {owner for value in matched.values() for owner in value}
    # A manifest/lock/build configuration change can affect every compiled consumer.
    if any(p in {"Cargo.toml", "Cargo.lock", "pyproject.toml"} or p.startswith(".cargo/") for p in normalized):
        affected.update(c["id"] for c in comps if "manifest" in c)
    changed = True
    while changed:
        changed = False
        for edge in graph["edges"]:
            if edge["kind"] == "depends_on" and edge["to"] in affected and edge["from"] not in affected:
                affected.add(edge["from"])
                changed = True
    recipe_ids = sorted({recipe for c in comps if c["id"] in affected for recipe in c["recipes"]})
    return {"paths": matched, "affected_components": sorted(affected), "unmapped_paths": [p for p in normalized if not matched[p]],
            "manual_review_components": [c["id"] for c in comps if c["id"] in affected and not c["recipes"]],
            "suggested_recipes": {key: policy["recipes"][key] for key in recipe_ids},
            "warning": "Suggestions only; no validation commands executed. Reverse declared dependencies over-approximate impact, not runtime calls or verified test coverage. Review feature flags, integration and scientific contracts separately."}


def context(root: Path, target: str, limit: int = 25) -> dict:
    graph = build_graph(root)
    comps = [node for node in graph["nodes"] if node["kind"] == "component"]
    policy = load_policy(root)
    component_ids = {c["id"] for c in comps}
    if target in component_ids:
        selected = {target}
        file_id = None
    else:
        path = normalize_path(root, target)
        selected = set(owners(path, comps, policy))
        file_id = "file:" + path
    files = sorted(node["path"] for node in graph["nodes"] if node["kind"] == "file" and selected.intersection(node["components"]))
    return {"target": target, "components": [c for c in comps if c["id"] in selected],
            "dependencies": [e for e in graph["edges"] if e["kind"] == "depends_on" and e["from"] in selected],
            "consumers": [e for e in graph["edges"] if e["kind"] == "depends_on" and e["to"] in selected],
            "file_imports": [e for e in graph["edges"] if e["kind"] == "imports" and file_id in (e["from"], e["to"])],
            "files": files[:limit], "total_files": len(files), "truncated": len(files) > limit,
            "note": "Read these sources to verify behavior. File ownership and static imports are navigation aids, not a call graph."}


def git_paths(root: Path) -> list[str]:
    paths = set()
    for args in (("diff", "--no-renames", "--name-only", "-z", "HEAD", "--"), ("ls-files", "--others", "--exclude-standard", "-z")):
        result = subprocess.run(["git", *args], cwd=root, check=True, stdout=subprocess.PIPE)
        paths.update(os.fsdecode(path) for path in result.stdout.split(b"\0") if path)
    return sorted(paths)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("generate", help="Update deterministic checked-in maps")
    sub.add_parser("check", help="Fail if checked-in maps differ from current source topology")
    sub.add_parser("doctor", help="Report tool availability without building or downloading")
    lookup = sub.add_parser("context", help="Retrieve a bounded file list and neighbors for a component or source path")
    lookup.add_argument("target")
    lookup.add_argument("--limit", type=int, default=25)
    change = sub.add_parser("impact", help="Suggest affected components and validation; never run tests")
    change.add_argument("paths", nargs="*")
    change.add_argument("--git", action="store_true", help="Include staged, unstaged and nonignored untracked paths against HEAD")
    args = parser.parse_args(argv)
    try:
        if args.command in {"generate", "check"}:
            stale = []
            for path, content in rendered(ROOT).items():
                destination = ROOT / path
                if args.command == "generate":
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_text(content)
                    print(f"Generated {path}")
                elif not destination.is_file() or destination.read_text() != content:
                    stale.append(path)
            if stale:
                print("Stale maps: " + ", ".join(stale) + ". Run python3 scripts/agent_tools.py generate.", file=sys.stderr)
                return 1
            if args.command == "check":
                print("Repository maps are current.")
        elif args.command == "context":
            if args.limit < 1:
                parser.error("context --limit must be positive")
            print(json.dumps(context(ROOT, args.target, args.limit), indent=2, sort_keys=True))
        elif args.command == "impact":
            paths = args.paths + (git_paths(ROOT) if args.git else [])
            if not paths:
                parser.error("impact requires paths or --git")
            print(json.dumps(impact(ROOT, paths), indent=2, sort_keys=True))
        else:
            print(json.dumps({"python": sys.version.split()[0], "tools": {tool: shutil.which(tool) for tool in ("git", "cargo", "rustc", "cc", "make", "curl", "tar", "clang", "maturin")},
                              "environment": {key: os.environ.get(key) for key in ("HMMER_TARBALL", "LIBCLANG_PATH", "VIRTUAL_ENV")},
                              "note": "Availability only, not build qualification. See docs/BUILDING.md; no network or build commands executed."}, indent=2))
    except (OSError, ValueError, KeyError, SyntaxError, subprocess.CalledProcessError) as exc:
        print(f"agent_tools: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
