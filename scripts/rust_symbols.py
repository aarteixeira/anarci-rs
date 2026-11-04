#!/usr/bin/env python3
"""Offline Rust source navigation, deliberately lexical rather than a compiler.

Python 3.11+, standard library only. No cargo, builds, downloads or persistent index.
See docs/agents/RUST_NAVIGATION.md for the JSON contract and limitations.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

SCHEMA = 1
LIMITATIONS = [
    "Lexical candidates only: no name/type resolution, call graph or compiler validation.",
    "No macro expansion, cfg evaluation, imports/reexport resolution or effective API reachability.",
    "Module paths are file-derived hints; #[path], include!, workspace/crate aliases are not resolved.",
    "API lists declared visibility, including items inside private modules, not a verified exported API.",
]
IDENT = re.compile(r"(?:r#)?[^\W\d]\w*|_\w*", re.UNICODE)
ITEMS = {"fn", "struct", "enum", "trait", "type", "const", "static", "mod", "union", "macro_rules"}
EXCLUDED = {"target", "build", "dist", "node_modules", ".git", ".venv"}


def tokenize(source: str) -> list[dict]:
    """Skip nested comments and all Rust string/char forms; retain source offsets.

    Lifetimes remain tokens; raw identifiers lose r# only in their normalized value.
    This is a tolerant lexer, not a validator of malformed Rust.
    """
    out = []
    i, n = 0, len(source)
    while i < n:
        start = i
        if source[i].isspace():
            i += 1
            continue
        if source.startswith("//", i):
            end = source.find("\n", i)
            i = n if end < 0 else end
            continue
        if source.startswith("/*", i):
            depth, i = 1, i + 2
            while i < n and depth:
                if source.startswith("/*", i):
                    depth, i = depth + 1, i + 2
                elif source.startswith("*/", i):
                    depth, i = depth - 1, i + 2
                else:
                    i += 1
            continue
        raw = re.match(r'(?:br|cr|r)(#*)"', source[i:])
        if raw:
            ending = '"' + raw[1]
            end = source.find(ending, i + raw.end())
            i = n if end < 0 else end + len(ending)
            continue
        quoted = re.match(r'(?:b|c)?"', source[i:])
        if quoted:
            i += quoted.end()
            while i < n:
                if source[i] == "\\":
                    i += 2
                elif source[i] == '"':
                    i += 1
                    break
                else:
                    i += 1
            continue
        # A char has one character or escape followed by a quote; 'a is a lifetime.
        char = re.match(r"(?:b)?'(?:[^'\\\n]|\\(?:u\{[0-9a-fA-F_]+\}|x[0-9a-fA-F]{2}|[^\n]))'", source[i:])
        if char:
            i += char.end()
            continue
        ident = IDENT.match(source, i)
        if ident:
            i = ident.end()
            text = ident[0]
            out.append({"value": text.removeprefix("r#"), "start": start, "end": i, "identifier": True})
            continue
        pair = source[i:i + 2]
        if pair in {"::", "->", "=>"}:
            i += 2
        else:
            i += 1
        out.append({"value": source[start:i], "start": start, "end": i, "identifier": False})
    return out


def inventory(root: Path) -> list[str]:
    """Only Git-visible .rs files; never traverse local ignored research directories."""
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", "*.rs"],
        cwd=root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    paths = set()
    for raw in result.stdout.split(b"\0"):
        if not raw:
            continue
        name = os.fsdecode(raw)
        path = root / name
        if EXCLUDED.intersection(Path(name).parts) or path.suffix != ".rs":
            continue
        if path.is_symlink() or not path.is_file():
            continue
        if any(parent.is_symlink() for parent in path.parents if parent != root and root in parent.parents):
            continue
        if not path.resolve().is_relative_to(root.resolve()):
            continue
        paths.add(name)
    return sorted(paths)


def file_module(path: str) -> list[str]:
    parts = list(Path(path).with_suffix("").parts)
    if "src" in parts:
        parts = parts[parts.index("src") + 1:]
    if parts and parts[-1] in {"lib", "main", "mod"}:
        parts.pop()
    return parts


def parse_source(path: str, source: str) -> tuple[list[dict], list[dict]]:
    tokens = tokenize(source)
    lines = [0] + [m.end() for m in re.finditer("\n", source)]
    values = [t["value"] for t in tokens]
    pairs, stack = {}, []
    for i, value in enumerate(values):
        if value in {"{", "(", "["}:
            stack.append(i)
        elif value in {"}", ")", "]"} and stack:
            opening = stack.pop()
            if values[opening] == {"}": "{", ")": "(", "]": "["}[value]:
                pairs[opening] = i
    definitions, declared = [], set()

    def location(index):
        token = tokens[index]
        line = bisect.bisect_right(lines, token["start"])
        return {"path": path, "line": line, "column": token["start"] - lines[line - 1] + 1}

    def text(a, b):
        return " ".join(source[tokens[a]["start"]:tokens[b - 1]["end"]].split())[:240] if b > a else ""

    def boundary(i, end):
        # Skip nested generic/argument syntax, including const generic blocks.
        angle = 0
        while i < end:
            value = values[i]
            if value == "<":
                angle += 1
            elif value == ">":
                angle = max(0, angle - 1)
            elif value in {"(", "["} or (value == "{" and angle):
                i = pairs.get(i, i) + 1
                continue
            elif not angle and value in {"{", ";", "="}:
                return i
            i += 1
        return end

    def add(kind, name, index, start, stop, scope, visibility, **extra):
        item = dict(kind=kind, name=name, qualified_name="::".join(scope + [name]),
                    visibility=visibility, signature=text(start, stop), **location(index), **extra)
        definitions.append(item)
        declared.add(index)
        return item

    def walk(begin, end, scope, container=None):
        i = begin
        while i < end:
            start, visibility = i, "private"
            # Attributes and qualifiers belong to the next declaration.
            while i < end and values[i] == "#":
                i += 1
                if i < end and values[i] == "!":
                    i += 1
                if i < end and values[i] == "[":
                    i = pairs.get(i, i) + 1
                else:
                    break
            if i < end and values[i] == "pub":
                visibility, i = "pub", i + 1
                if i < end and values[i] == "(":
                    close = pairs.get(i, i)
                    visibility = "pub(" + "".join(values[i + 1:close]) + ")"
                    i = close + 1
            while i < end and values[i] in {"async", "unsafe", "default", "auto", "extern"}:
                i += 1
            # const fn and const unsafe fn are functions, not const declarations.
            if i < end and values[i] == "const":
                look = i + 1
                while look < end and values[look] in {"unsafe", "async", "extern"}:
                    look += 1
                if look < end and values[look] == "fn":
                    i = look
            if i >= end:
                break
            kind = values[i]
            if kind == "impl":
                stop = boundary(i + 1, end)
                header = values[i + 1:stop]
                # Remove leading impl generic parameters, preserving target generics.
                offset, depth = 0, 0
                if header and header[0] == "<":
                    for offset, value in enumerate(header):
                        depth += (value == "<") - (value == ">")
                        if depth == 0:
                            offset += 1
                            break
                rest = header[offset:]
                if "where" in rest:
                    rest = rest[:rest.index("where")]
                # 'for' inside an HRTB binder is not the trait/target separator.
                separator = next((j for j, v in enumerate(rest) if v == "for" and (j + 1 == len(rest) or rest[j + 1] != "<")), None)
                trait = " ".join(rest[:separator]) if separator is not None else None
                target = " ".join(rest[separator + 1:] if separator is not None else rest)
                name = text(i, stop)
                add("impl", name, i, start, stop, scope, "private", trait=trait, target=target)
                if stop < end and values[stop] == "{":
                    close = pairs.get(stop, end - 1)
                    walk(stop + 1, close, scope + [name], "impl")
                    i = close + 1
                else:
                    i = stop + 1
                continue
            if kind in ITEMS:
                name_index = i + (2 if kind == "macro_rules" and i + 1 < end and values[i + 1] == "!" else 1)
                if kind == "static" and name_index < end and values[name_index] == "mut":
                    name_index += 1
                if name_index >= end or not tokens[name_index]["identifier"]:
                    i += 1
                    continue
                name = values[name_index]
                stop = boundary(name_index + 1, end)
                if container == "trait" and visibility == "private":
                    visibility = "trait_member"
                add(kind, name, name_index, start, stop, scope, visibility, container=container)
                if stop < end and values[stop] == "{":
                    close = pairs.get(stop, end - 1)
                    if kind in {"mod", "trait"}:
                        walk(stop + 1, close, scope + [name], kind)
                    i = close + 1
                else:
                    # Initializers may contain arbitrary braces and identifiers.
                    i = stop + 1
                    if stop < end and values[stop] == "=":
                        while i < end and values[i] != ";":
                            i = pairs.get(i, i) + 1
                        i += 1
                continue
            if kind == "use":
                stop = i + 1
                while stop < end and values[stop] != ";":
                    stop = pairs.get(stop, stop) + 1
                if visibility != "private":
                    add("reexport", text(i + 1, stop), i, start, stop, scope, visibility)
                i = stop + 1
                continue
            # Skip unknown blocks (macro bodies, extern blocks, expressions), never
            # advertise their apparent declarations as ordinary source definitions.
            if values[i] in {"{", "(", "["}:
                i = pairs.get(i, i) + 1
            else:
                i += 1

    walk(0, len(tokens), file_module(path))
    occurrences = []
    for i, token in enumerate(tokens):
        if token["identifier"] and i not in declared:
            occurrences.append(dict(name=token["value"], **location(i)))
    return definitions, occurrences


def build_index(root: Path) -> dict:
    digest = hashlib.sha256(f"rust-symbols:{SCHEMA}\0".encode())
    definitions, occurrences, files = [], [], inventory(root)
    for path in files:
        raw = (root / path).read_bytes()
        digest.update(path.encode("utf-8", "surrogateescape") + b"\0" + raw + b"\0")
        defs, refs = parse_source(path, raw.decode("utf-8"))
        definitions.extend(defs)
        occurrences.extend(refs)
    return {"schema_version": SCHEMA, "source_digest": digest.hexdigest(), "files": files,
            "definitions": definitions, "occurrences": occurrences}


def query_index(index: dict, command: str, query: str = "", *, limit: int = 30,
                offset: int = 0, path: str = "", include_restricted: bool = False) -> dict:
    if not 1 <= limit <= 200 or not 0 <= offset <= 1_000_000:
        raise ValueError("limit must be 1..200; offset must be 0..1000000")
    if command not in {"search", "definition", "references", "implementations", "api"}:
        raise ValueError("unknown navigation command")
    if command != "api" and not query:
        raise ValueError("query must not be empty")
    rows = index["occurrences" if command == "references" else "definitions"]
    matches = []
    for row in rows:
        if path and not row["path"].startswith(path):
            continue
        name = row["name"]
        qualified = row.get("qualified_name", name)
        if command == "references":
            match = name == query.split("::")[-1].removeprefix("r#")
        elif command == "definition":
            normalized = "::".join(part.removeprefix("r#") for part in query.split("::"))
            match = row["kind"] != "impl" and normalized in {name, qualified}
        elif command == "implementations":
            match = row["kind"] == "impl" and any(
                "".join(query.split()) in "".join(part.split())
                for part in [row.get("trait") or "", row.get("target") or ""]
            )
        elif command == "api":
            match = (row["visibility"] in {"pub", "trait_member"} or
                     include_restricted and row["visibility"].startswith("pub(")) and query in qualified
        else:
            match = query in qualified
        if match:
            matches.append(row)
    matches.sort(key=lambda row: (row["path"], row["line"], row["column"], row["name"]))
    return {"schema_version": SCHEMA, "mode": "lexical", "semantics": "unresolved_candidates",
            "command": command, "query": query, "source_digest": index["source_digest"],
            "files_indexed": len(index["files"]), "total": len(matches), "offset": offset,
            "limit": limit, "truncated": offset + limit < len(matches),
            "results": matches[offset:offset + limit], "limitations": LIMITATIONS}


def semantic_status(executable=None) -> dict:
    """Probe installed tooling only. Never silently substitute lexical candidates."""
    executable = executable or shutil.which("rust-analyzer")
    detail = "rust-analyzer was not found on PATH"
    if executable:
        try:
            probe = subprocess.run([executable, "--version"], capture_output=True, text=True,
                                   timeout=5, env={**os.environ, "CARGO_NET_OFFLINE": "true"})
            detail = (probe.stdout + probe.stderr).strip()[:1200]
            if probe.returncode == 0:
                return {"schema_version": SCHEMA, "mode": "semantic", "status": "available",
                        "semantics": "not_queried", "results": [], "analyzer": executable,
                        "reason": detail}
        except (OSError, subprocess.TimeoutExpired) as exc:
            detail = str(exc)[:1200]
    return {"schema_version": SCHEMA, "mode": "semantic", "status": "unavailable",
            "semantics": "not_queried", "results": [], "analyzer": executable,
            "reason": detail, "limitations": ["No lexical fallback is silently substituted for semantic results.",
            "Use lexical mode for source candidates. No LSP, cargo, build scripts, proc macros or check-on-save are started."]}


def semantic_query(root, command, at, executable, timeout):
    path, line, column = at.rsplit(":", 2)
    specification = importlib.util.spec_from_file_location("rust_lsp", Path(__file__).with_name("rust_lsp.py"))
    adapter = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(adapter)
    return adapter.query(root.resolve(), command, path, int(line), int(column), executable=executable, timeout=timeout)


class JsonParser(argparse.ArgumentParser):
    def error(self, message):
        print(json.dumps({"schema_version": SCHEMA, "status": "error", "error": message}))
        raise SystemExit(2)


def main(argv=None) -> int:
    parser = JsonParser(description=__doc__)
    parser.add_argument("command", choices=["search", "definition", "references", "implementations", "api"])
    parser.add_argument("query", nargs="?", default="")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--path", default="", help="repository-relative path prefix")
    parser.add_argument("--limit", type=int, default=30)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--include-restricted", action="store_true", help="include pub(crate)/pub(super)/pub(in ...) in api")
    parser.add_argument("--mode", choices=["lexical", "semantic"], default="lexical")
    parser.add_argument("--at", help="semantic request location file:line:column (1-based codepoints)")
    parser.add_argument("--analyzer", help="explicit installed rust-analyzer executable for semantic mode")
    parser.add_argument("--timeout", type=float, default=30, help="semantic LSP deadline in seconds")
    args = parser.parse_args(argv)
    try:
        if not 1 <= args.limit <= 200 or not 0 <= args.offset <= 1_000_000:
            raise ValueError("limit must be 1..200; offset must be 0..1000000")
        if args.mode == "semantic":
            if args.command not in {"definition", "references", "implementations"}:
                raise ValueError("semantic requests support definition/references/implementations only")
            if not args.at or not re.fullmatch(r".+:[1-9][0-9]*:[1-9][0-9]*", args.at):
                raise ValueError("semantic requests require --at file:line:column (1-based)")
            if args.query:
                raise ValueError("semantic mode resolves --at; omit the positional symbol query")
            result = semantic_status(args.analyzer)
            if result["status"] == "available":
                result = semantic_query(args.root, args.command, args.at, result["analyzer"], args.timeout)
                matches = result["results"]
                result.update(total=len(matches), limit=args.limit, offset=args.offset,
                              truncated=args.offset + args.limit < len(matches),
                              results=matches[args.offset:args.offset + args.limit])
            result.update(command=args.command, at=args.at)
            print(json.dumps(result, sort_keys=True))
            return 0 if result["status"] == "ok" else 3
        if args.at:
            raise ValueError("--at requires --mode semantic; lexical mode cannot resolve locations")
        index = build_index(args.root.resolve())
        result = query_index(index, args.command, args.query, limit=args.limit, offset=args.offset,
                             path=args.path, include_restricted=args.include_restricted)
        print(json.dumps(result, sort_keys=True))
        return 0
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(json.dumps({"schema_version": SCHEMA, "status": "error", "error": str(exc)[:1200]}))
        return 2


if __name__ == "__main__":
    sys.exit(main())
