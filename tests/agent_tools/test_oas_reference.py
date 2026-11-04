"""Offline synthetic parser contracts, not biological reference evidence."""
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('oas_import', Path(__file__).resolve().parents[2] / 'scripts/build_oas_reference.py')
oas = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(oas)


class OASImportTests(unittest.TestCase):
    def row(self):
        return {'productive': 'T', 'stop_codon': 'F', 'vj_in_frame': 'T',
                'sequence_alignment_aa': 'A' * 128, 'locus': 'H', 'cdr3_aa': 'AAA',
                'v_call': 'IGHV1-2*01',
                'ANARCI_numbering': repr({'fw': {str(i): 'A' for i in range(1, 129)}})}

    def unit(self):
        return {'url': 'https://opig.stats.ox.ac.uk/test.csv.gz', 'study': 'study1',
                'species': 'human', 'domain_format': 'VH', 'split': 'train',
                'biological_origin': 'natural_repertoire', 'library_selection': 'unselected',
                'origin_evidence': {'url': 'https://example.org/offline-test-only',
                                    'review_note': 'Synthetic parser fixture assertion only; NOT real biological evidence.'}}

    def meta(self):
        return {'Species': 'human', 'Subject': 'no', 'Disease': 'None', 'Vaccine': 'None'}

    def test_metadata_csv_quoted_json(self):
        out = io.StringIO(); csv.writer(out).writerow([json.dumps(self.meta())])
        self.assertEqual(oas.unit_metadata(io.StringIO(out.getvalue())), self.meta())
        self.assertEqual(oas.unit_metadata(io.StringIO(json.dumps(self.meta()))), self.meta())

    def test_numbering_safe_and_bounded(self):
        self.assertEqual(oas.parse_numbering("{'fwh1': {'1 ': 'Q', '111 A': 'G'}}"), {'1': 'Q', '111A': 'G'})
        for bad in ["__import__('os').system('false')", "{'fw': {'1':'X'}}", 'x' * 30001, "{'a': {'1':'Q'},'b':{'1':'Q'}}"]:
            with self.assertRaises((ValueError, SyntaxError)):
                oas.parse_numbering(bad)

    def test_normalization_missing_normal_imgt_positions(self):
        row = self.row(); row['ANARCI_numbering'] = repr({'fw': {str(i): 'A' for i in range(1,129) if i not in {10,73}}})
        rec = oas.normalize(row, self.meta(), self.unit(), 1)
        self.assertTrue(rec['donor_unknown']); self.assertEqual(rec['donor'], 'unknown')
        self.assertEqual(rec['domain_format'], 'VH')
        row['stop_codon'] = 'T'
        with self.assertRaisesRegex(ValueError, 'nonfunctional'):
            oas.normalize(row, self.meta(), self.unit(), 1)

    def test_nano_anarcii_column(self):
        row = self.row(); row['ANARCII_numbering'] = row.pop('ANARCI_numbering')
        self.assertEqual(len(oas.normalize(row, self.meta(), self.unit(), 1)['numbering']), 128)

    def test_imgt_insertion_order_and_camelid_v_family(self):
        keys = ['112','111B','112A','111','112B','111A','113']
        self.assertEqual(sorted(keys,key=oas.imgt_position_key), ['111','111A','111B','112B','112A','112','113'])
        for center in (33,61):
            keys = [str(center),f'{center}A',f'{center}B']
            self.assertEqual(sorted(keys,key=oas.imgt_position_key), [f'{center}B',f'{center}A',str(center)])
        a = self.row(); b = self.row(); a['v_call']='IGHV3S60*01'; b['v_call']='IGHV3S61*01'
        self.assertEqual(oas.normalize(a,self.meta(),self.unit(),1)['lineage'], oas.normalize(b,self.meta(),self.unit(),2)['lineage'])

    def test_nano_uses_full_sequence_not_shortened_alignment(self):
        row = self.row(); row['sequence_aa'] = 'A'*128; row['sequence_alignment_aa'] = 'A'*120
        result = oas.normalize(row, self.meta(), self.unit(), 1)
        self.assertEqual(result['sequence'], 'A'*128)
        self.assertEqual(result['sequence_source_field'], 'sequence_aa')
        row['sequence_aa'] = 'C'*128
        with self.assertRaisesRegex(ValueError, 'numbering_sequence_mismatch'):
            oas.normalize(row, self.meta(), self.unit(), 1)

    def test_vhh_not_inferred_from_species(self):
        unit = self.unit(); unit.update(species='alpaca', domain_format='VHH')
        meta = self.meta(); meta['Species'] = 'alpaca'
        with self.assertRaisesRegex(ValueError, 'VHH_requires_evidence'):
            oas.normalize(self.row(), meta, unit, 1)
        unit['vhh_evidence'] = 'https://oas-nano.opig.stats.ox.ac.uk/documentation/'
        self.assertEqual(oas.normalize(self.row(), meta, unit, 1)['domain_format'], 'VHH')

    def test_dedup_priority_and_study_isolation(self):
        r = oas.normalize(self.row(), self.meta(), self.unit(), 1)
        other = dict(r, study='study2', split='test', record_id='other')
        kept, audit = oas.deduplicate([other,r])
        self.assertEqual(kept, [r]); self.assertEqual(audit['exact_duplicates_removed'], 1)
        other['sequence'] = 'C'*128
        kept, audit = oas.deduplicate([other,r])
        self.assertEqual(kept,[r]); self.assertEqual(audit['cross_split_lineage_removed'],1)
        other['study'] = r['study']
        with self.assertRaisesRegex(ValueError, 'study occurs'):
            oas.deduplicate([other,r])

    def test_complete_unit_reservoir_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'unit.gz'; row = self.row()
            with gzip.open(path,'wt',newline='') as f:
                csv.writer(f).writerow([json.dumps(self.meta())]); writer=csv.DictWriter(f,fieldnames=list(row));writer.writeheader()
                for _ in range(80): writer.writerow(row)
            unit = dict(self.unit(), path=str(path), complete=True, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            a,counts,_ = oas.sample_unit(unit,5); b,_,_=oas.sample_unit(unit,5)
            self.assertEqual(a,b);self.assertEqual(counts['rows'],80)
            self.assertTrue(any(r['source_row']>5 for r in a))
            unit['sha256']='0'*64
            with self.assertRaisesRegex(ValueError,'checksum'):
                oas.sample_unit(unit)
            unit['complete']=False
            with self.assertRaisesRegex(ValueError,'complete'):
                oas.sample_unit(unit)

    def test_download_full_size_cap_and_pinned_hash(self):
        class Response(io.BytesIO):
            headers = {'Content-Length': '3'}
            def geturl(self): return 'https://opig.stats.ox.ac.uk/test.csv.gz'
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(oas.urllib.request, 'urlopen', return_value=Response(b'')):
                with self.assertRaisesRegex(ValueError, 'oversized'):
                    oas.download_unit(self.unit(), tmp, 2)
            unit = dict(self.unit(), sha256=hashlib.sha256(b'abc').hexdigest())
            with patch.object(oas.urllib.request, 'urlopen', side_effect=[Response(b''),Response(b'abc')]):
                acquired = oas.download_unit(unit, tmp, 3)
            self.assertTrue(acquired['complete']); self.assertEqual(acquired['compressed_bytes'],3)
        with tempfile.TemporaryDirectory() as tmp:
            unit['sha256'] = '0'*64
            with patch.object(oas.urllib.request, 'urlopen', side_effect=[Response(b''),Response(b'abc')]):
                with self.assertRaisesRegex(ValueError, 'pinned local checksum'):
                    oas.download_unit(unit,tmp,3)

    def test_output_preflight_and_exclusive_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/'records.jsonl'; audit = Path(str(output)+'.audit.json')
            audit.write_text('preserve')
            with patch.object(oas, 'sample_unit') as sample:
                with patch('sys.argv',['build_oas_reference.py','prepare','missing.json','--output',str(output)]):
                    with self.assertRaisesRegex(ValueError,'overwrite'):
                        oas.main()
                sample.assert_not_called()
            self.assertEqual(audit.read_text(),'preserve'); self.assertFalse(output.exists())
            with self.assertRaisesRegex(ValueError,'input and output'):
                oas.preflight_outputs([output], output)
            oas.write_new(output,'first')
            with self.assertRaises(FileExistsError): oas.write_new(output,'second')
            self.assertEqual(output.read_text(),'first')

    def test_manifest_split_validation_and_position_caps(self):
        unit = self.unit(); unit['split']='validation'
        with self.assertRaisesRegex(ValueError,'invalid split'):
            oas.validate_manifest({'units':[unit]})
        unit=self.unit(); second=dict(unit,url='https://opig.stats.ox.ac.uk/other.csv.gz',split='test')
        with self.assertRaisesRegex(ValueError,'multiple splits'):
            oas.validate_manifest({'units':[unit,second]})
        for position in ('129','999','112ABCDE'):
            with self.assertRaisesRegex(ValueError,'cap'):
                oas.parse_numbering(repr({'fw':{position:'A'}}))
        with self.assertRaisesRegex(ValueError,'positive'):
            oas.download_unit(self.unit(),'.',0)
        with self.assertRaisesRegex(ValueError,'positive'):
            oas.sample_unit(self.unit(),0)

    def test_origin_gate_before_network_or_local_data(self):
        for field, value in [('biological_origin','synthetic_library'),('biological_origin','unresolved'),('library_selection','antigen_panned'),('origin_evidence',None)]:
            unit = dict(self.unit()); unit[field]=value
            with patch.object(oas.urllib.request,'urlopen') as request:
                with self.assertRaises(ValueError): oas.download_unit(unit,'.',100)
                request.assert_not_called()
            with self.assertRaises(ValueError): oas.sample_unit(unit)
            with self.assertRaises(ValueError): oas.validate_manifest({'units':[unit]})
        unit=self.unit();unit.pop('biological_origin')
        with self.assertRaises(ValueError): oas.normalize(self.row(),self.meta(),unit,1)
        rec=oas.normalize(self.row(),self.meta(),self.unit(),1)
        self.assertEqual(rec['origin_evidence'],self.unit()['origin_evidence'])

    def test_restrict_download_hosts(self):
        for url in ['http://opig.stats.ox.ac.uk/a', 'https://example.org/a', 'file:///tmp/a']:
            with self.assertRaises(ValueError): oas.checked_url(url)


if __name__ == '__main__':
    unittest.main()
