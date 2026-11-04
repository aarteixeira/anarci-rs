#!/usr/bin/env python3
"""Bounded, stdlib-only OAS acquisition and normalization (Python 3.11+).

Commands: discover --species human --output index.json;
 download manifest.json --output acquired.json; prepare acquired.json --output records.jsonl.
Manifests contain a ``units`` list with url, study, split, species, domain_format,
biological_origin='natural_repertoire', library_selection='unselected', and
origin_evidence={url, review_note} recording primary-source eligibility review.
VHH units also require vhh_evidence. OAS species/health labels alone do not
establish biological origin or clinical cohort; review the primary study.
Downloads are complete units, never shell scripts;
oversized units are rejected. A local SHA256 is not an upstream published hash.
"""
import argparse
import ast
from collections import Counter
import csv
import gzip
import hashlib
import html
import json
from pathlib import Path
import random
import re
import urllib.parse
import urllib.request

AA = frozenset('ACDEFGHIKLMNPQRSTVWY')
FRAMEWORK = frozenset(map(str, [*range(1, 27), *range(39, 56), *range(66, 105), *range(118, 129)]))
HOSTS = {'opig.stats.ox.ac.uk', 'oas-nano.opig.stats.ox.ac.uk'}
UNKNOWN = {'', 'no', 'none', 'unknown', 'undefined', 'nan'}


def checked_url(url):
    p = urllib.parse.urlsplit(url)
    if p.scheme != 'https' or p.hostname not in HOSTS or p.username or p.password:
        raise ValueError('only official HTTPS OAS URLs allowed')
    return url


def discover(species, nano=False):
    url = ('https://oas-nano.opig.stats.ox.ac.uk/search/' if nano else
           'https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/')
    fields = dict.fromkeys(['Species', 'BSource', 'BType', 'Longitudinal', 'Age',
                           'Disease', 'Subject', 'Vaccine', 'Chain', 'Isotype', 'Primer'], '*')
    fields.update(Species=species, Disease='None', Vaccine='None', Chain='Heavy')
    with urllib.request.urlopen(url, urllib.parse.urlencode(fields).encode(), timeout=90) as r:
        text = r.read(20_000_001).decode()
    if len(text) > 20_000_000:
        raise ValueError('metadata page exceeds cap')
    urls = re.findall(r'wget (https://[^\s"\']+)', text)
    match = re.search(r'var downloadPaths = (\[.*?\]);', text)
    if match:
        urls += [urllib.parse.urljoin(url, x) for x in json.loads(match[1])]
    return {'query_url': url, 'query': fields, 'page_sha256': hashlib.sha256(text.encode()).hexdigest(),
            'urls': sorted(set(checked_url(html.unescape(u)) for u in urls)), 'page': text}


def validate_origin(unit):
    """Require an explicit primary-source review, not OAS species labels alone."""
    if unit.get('biological_origin') != 'natural_repertoire':
        raise ValueError('verified natural repertoire origin required')
    if unit.get('library_selection') != 'unselected':
        raise ValueError('unselected repertoire library required')
    evidence = unit.get('origin_evidence')
    if not isinstance(evidence, dict) or not isinstance(evidence.get('review_note'), str) or not evidence['review_note'].strip():
        raise ValueError('origin review evidence required')
    url = evidence.get('url', '')
    if not isinstance(url, str) or urllib.parse.urlsplit(url).scheme != 'https' or not urllib.parse.urlsplit(url).hostname:
        raise ValueError('origin evidence requires HTTPS primary-source URL')


def download_unit(unit, directory, cap):
    validate_origin(unit)
    if cap <= 0:
        raise ValueError('download cap must be positive')
    url = checked_url(unit['url'])
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / (hashlib.sha256(url.encode()).hexdigest()[:16] + '.csv.gz')
    if path.exists():
        raise ValueError(f'refusing to overwrite existing download: {path}')
    with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=60) as r:
        checked_url(r.geturl())
        size = int(r.headers.get('Content-Length', '-1'))
    if size < 0 or size > cap:
        raise ValueError(f'unknown/oversized full unit: {size} bytes, cap {cap}')
    digest = hashlib.sha256()
    count = 0
    with urllib.request.urlopen(url, timeout=90) as r, path.open('xb') as out:
        checked_url(r.geturl())
        while count < size:
            block = r.read(min(1024 * 1024, size - count))
            if not block:
                break
            digest.update(block)
            out.write(block)
            count += len(block)
        extra = r.read(1)
    if count != size or extra:
        raise ValueError(f'incomplete/changed unit, do not import: {path}')
    if unit.get('sha256') and digest.hexdigest() != unit['sha256']:
        raise ValueError(f'upstream changed relative to pinned local checksum: {path}')
    return {**unit, 'path': str(path), 'compressed_bytes': count,
            'sha256': digest.hexdigest(), 'complete': True,
            'hash_scope': 'complete downloaded compressed file; not upstream-published hash'}


def imgt_position_key(position):
    """IMGT 33/61/112 insertions precede base in reverse order."""
    match = re.fullmatch(r'([1-9][0-9]{0,2})([A-Z]*)', position)
    if not match:
        raise ValueError('invalid IMGT position')
    number, insertion = int(match[1]), match[2]
    if number > 128 or len(insertion) > 4:
        raise ValueError('IMGT position or insertion exceeds cap')
    if number in {33, 61, 112} and insertion:
        return number, -1, -len(insertion), tuple(-ord(c) for c in insertion)
    return number, 0, len(insertion), tuple(ord(c) for c in insertion)


def parse_numbering(value):
    if len(value) > 30_000:
        raise ValueError('numbering exceeds size cap')
    tree = ast.parse(value, mode='eval')
    if sum(1 for _ in ast.walk(tree)) > 2500:
        raise ValueError('numbering exceeds AST cap')
    raw = ast.literal_eval(tree)
    if not isinstance(raw, dict) or len(raw) > 10:
        raise ValueError('expected region dictionaries')
    result = {}
    for region in raw.values():
        if not isinstance(region, dict) or len(region) > 200:
            raise ValueError('invalid numbering region')
        for position, aa in region.items():
            key = str(position).replace(' ', '')
            if not re.fullmatch(r'[1-9][0-9]{0,2}[A-Z]*', key) or not isinstance(aa, str) or aa not in AA:
                raise ValueError('invalid position or residue')
            if key in result:
                raise ValueError('duplicate position')
            result[key] = aa
    return dict(sorted(result.items(), key=lambda item: imgt_position_key(item[0])))


def unit_metadata(handle):
    first = handle.readline(100_001)
    if len(first) > 100_000:
        raise ValueError('metadata exceeds size cap')
    # OAS stores its JSON metadata as one CSV-quoted field.
    value = first.strip()
    if not value.startswith('{'):
        value = next(csv.reader([value]))[0]
    result = json.loads(value)
    if not isinstance(result, dict):
        raise ValueError('expected metadata object')
    return result


def normalize(row, meta, unit, index):
    validate_origin(unit)
    if row.get('productive') != 'T' or row.get('stop_codon') != 'F' or row.get('vj_in_frame') != 'T':
        raise ValueError('nonfunctional')
    sequence_field = 'sequence_aa' if row.get('sequence_aa') else 'sequence_alignment_aa'
    seq = row.get(sequence_field, '').replace('-', '').replace('.', '')
    if not 90 <= len(seq) <= 180 or not set(seq) <= AA:
        raise ValueError('noncanonical_or_short')
    numbering = parse_numbering(row.get('ANARCI_numbering') or row.get('ANARCII_numbering', ''))
    numbered_sequence = ''.join(numbering.values())
    if numbered_sequence not in seq:
        raise ValueError('numbering_sequence_mismatch')
    # Use the observed numbered domain, excluding unnumbered terminal overhangs.
    seq = numbered_sequence
    if not 90 <= len(seq) <= 180:
        raise ValueError('numbered_domain_length')
    if len(FRAMEWORK.intersection(numbering)) < 60 or not {'23', '104', '118', '128'} <= numbering.keys():
        raise ValueError('incomplete_framework')
    # Require near-full domain; normal IMGT gaps (e.g. 10 and 73) are allowed.
    if not {'2', '3', '4', '5', '6', '7', '8', '9'} <= numbering.keys():
        raise ValueError('truncated_framework')
    chain = {'IGH': 'H', 'IGK': 'K', 'IGL': 'L'}.get(row.get('locus'), row.get('locus'))
    if chain not in {'H', 'K', 'L'}:
        raise ValueError('unknown_chain')
    domain = unit.get('domain_format', 'VH' if chain == 'H' else 'VL')
    if domain not in ({'VH', 'VHH'} if chain == 'H' else {'VL'}):
        raise ValueError('chain_domain_mismatch')
    if domain == 'VHH' and not (unit.get('vhh_evidence') and chain == 'H'):
        raise ValueError('VHH_requires_evidence')
    species_raw = str(meta.get('Species', ''))
    species = 'mouse' if species_raw.startswith('mouse') else species_raw
    if species != unit['species']:
        raise ValueError('species_metadata_mismatch')
    donor = meta.get('Subject')
    unknown = str(donor).lower() in UNKNOWN
    cdr3 = row.get('cdr3_aa', '')
    # Camelid genes use IGHV3S60 as well as human/mouse IGHV3-23 style.
    family_match = re.match(r'(IG[HKL]V[0-9]+)', row.get('v_call', ''))
    vfamily = family_match[1] if family_match else ''
    if not cdr3 or not vfamily or not set(cdr3) <= AA:
        raise ValueError('missing_lineage_proxy')
    lineage = hashlib.sha256((vfamily + '|' + cdr3).encode()).hexdigest()
    return {'record_id': hashlib.sha256((unit['url'] + ':' + str(index)).encode()).hexdigest(),
            'species': species, 'chain': chain, 'domain_format': domain, 'study': unit['study'],
            'donor': 'unknown' if unknown else donor, 'donor_unknown': unknown, 'lineage': lineage,
            'lineage_method': 'exact_CDR3_AA_plus_V_family_proxy', 'split': unit['split'],
            'sequence': seq, 'numbering': numbering, 'source_url': unit['url'],
            'source_row': index, 'source_metadata': meta, 'vhh_evidence': unit.get('vhh_evidence'),
            'numbering_source': 'ANARCI_numbering' if row.get('ANARCI_numbering') else 'ANARCII_numbering',
            'numbering_scheme': 'IMGT', 'near_full_framework_required': True,
            'sequence_source_field': sequence_field, 'sequence_scope': 'numbered domain; verified source substring',
            'biological_origin': unit['biological_origin'], 'library_selection': unit['library_selection'],
            'origin_evidence': unit['origin_evidence']}


def sample_unit(unit, sample_size=2000, seed=1729):
    validate_origin(unit)
    if sample_size <= 0:
        raise ValueError('sample size must be positive')
    if not unit.get('complete'):
        raise ValueError('only complete downloaded units can be prepared')
    path = Path(unit['path'])
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(block)
    if digest.hexdigest() != unit['sha256']:
        raise ValueError('download checksum mismatch')
    rng = random.Random(str(seed) + unit['url'])
    sample, counts = [], Counter()
    with gzip.open(path, 'rt', newline='') as f:
        meta = unit_metadata(f)
        if str(meta.get('Vaccine', '')).lower() != 'none' or str(meta.get('Disease', '')).lower() != 'none':
            raise ValueError('unit must be healthy and unvaccinated')
        for i, row in enumerate(csv.DictReader(f), 1):
            counts['rows'] += 1
            try:
                record = normalize(row, meta, unit, i)
            except (ValueError, SyntaxError, TypeError, RecursionError) as exc:
                counts['excluded:' + str(exc)] += 1
                continue
            counts['eligible'] += 1
            if len(sample) < sample_size:
                sample.append(record)
            else:
                j = rng.randrange(counts['eligible'])
                if j < sample_size:
                    sample[j] = record
    counts['sampled'] = len(sample)
    return sample, dict(counts), meta


def deduplicate(records):
    """Train-first priority; prune exact and lineage duplicates across splits."""
    sequences, lineages, kept, counts = set(), {}, [], Counter()
    priority = {'train': 0, 'calibration': 1, 'test': 2}
    study_splits = {}
    for r in records:
        key = r['study']
        if key in study_splits and study_splits[key] != r['split']:
            raise ValueError('study occurs in multiple splits')
        study_splits[key] = r['split']
    for r in sorted(records, key=lambda x: (priority[x['split']], x['record_id'])):
        seq = r['sequence']
        if seq in sequences:
            counts['exact_duplicates_removed'] += 1
            continue
        if r['lineage'] in lineages and lineages[r['lineage']] != r['split']:
            counts['cross_split_lineage_removed'] += 1
            continue
        sequences.add(seq)
        lineages[r['lineage']] = r['split']
        kept.append(r)
    return kept, dict(counts)


def validate_manifest(manifest):
    units = manifest.get('units')
    if not isinstance(units, list) or not units:
        raise ValueError('manifest requires nonempty units list')
    studies, urls = {}, set()
    for unit in units:
        if not isinstance(unit, dict):
            raise ValueError('unit must be an object')
        for key in ('url', 'study', 'species', 'split', 'domain_format'):
            if not isinstance(unit.get(key), str) or not unit[key].strip():
                raise ValueError(f'unit requires {key}')
        checked_url(unit['url'])
        validate_origin(unit)
        if unit['split'] not in {'train', 'calibration', 'test'}:
            raise ValueError('invalid split')
        if unit['domain_format'] not in {'VH', 'VL', 'VHH'}:
            raise ValueError('invalid domain format')
        if unit['study'] in studies and studies[unit['study']] != unit['split']:
            raise ValueError('study occurs in multiple splits')
        studies[unit['study']] = unit['split']
        if unit['url'] in urls:
            raise ValueError('duplicate unit URL')
        urls.add(unit['url'])


def preflight_outputs(paths, input_path=None):
    resolved = [Path(path).resolve() for path in paths]
    if len(set(resolved)) != len(resolved):
        raise ValueError('output paths must be distinct')
    if input_path and Path(input_path).resolve() in resolved:
        raise ValueError('input and output paths must differ')
    for path in resolved:
        if path.exists() or path.is_symlink():
            raise ValueError(f'refusing to overwrite output: {path}')
        if not path.parent.is_dir():
            raise ValueError(f'output parent must exist: {path.parent}')


def write_new(path, text):
    with Path(path).open('x', encoding='utf-8') as handle:
        handle.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    d = sub.add_parser('discover'); d.add_argument('--species', default='human'); d.add_argument('--nano', action='store_true'); d.add_argument('--output', required=True)
    d = sub.add_parser('download'); d.add_argument('manifest'); d.add_argument('--output', required=True); d.add_argument('--directory', default='outputs/species-reference/raw'); d.add_argument('--max-total-bytes', type=int, default=500_000_000); d.add_argument('--max-unit-bytes', type=int, default=50_000_000)
    d = sub.add_parser('prepare'); d.add_argument('manifest'); d.add_argument('--output', required=True); d.add_argument('--sample-size', type=int, default=2000); d.add_argument('--seed', type=int, default=1729)
    args = parser.parse_args()
    outputs = [args.output]
    if args.command == 'prepare':
        outputs.append(args.output + '.audit.json')
        if args.sample_size <= 0:
            parser.error('sample size must be positive')
    if args.command == 'download' and (args.max_total_bytes <= 0 or args.max_unit_bytes <= 0):
        parser.error('download caps must be positive')
    preflight_outputs(outputs, getattr(args, 'manifest', None))
    if args.command == 'discover':
        result = discover(args.species, args.nano)
    else:
        manifest = json.loads(Path(args.manifest).read_text())
        validate_manifest(manifest)
        if args.command == 'download':
            acquired, used = [], 0
            for unit in manifest['units']:
                item = download_unit(unit, args.directory, min(args.max_unit_bytes, args.max_total_bytes - used))
                acquired.append(item); used += item['compressed_bytes']
            result = {**manifest, 'units': acquired, 'compressed_bytes': used}
        else:
            if args.sample_size <= 0:
                parser.error('sample size must be positive')
            records, audits = [], []
            for unit in manifest['units']:
                sample, counts, metadata = sample_unit(unit, args.sample_size, args.seed)
                records.extend(sample); audits.append({'url': unit['url'], 'counts': counts, 'metadata': metadata})
            records, dedup = deduplicate(records)
            write_new(args.output, ''.join(json.dumps(r, sort_keys=True) + '\n' for r in records))
            result = {'units': audits, 'deduplication': dedup, 'records': len(records),
                      'sampling': 'deterministic per-unit reservoir over all eligible rows of complete files',
                      'seed': args.seed, 'sample_size': args.sample_size,
                      'limitations': ['small purposively chosen units; not representative of all OAS',
                                      'lineage proxy is not inferred clonal ancestry',
                                      'unknown donors retained and flagged; study-level splits only']}
            write_new(args.output + '.audit.json', json.dumps(result, indent=2))
            return
    write_new(args.output, json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
