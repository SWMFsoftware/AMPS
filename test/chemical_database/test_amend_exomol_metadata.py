"""Native-evidence amendment controls; no constructed science observations."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'utility/chemical_database'))
spec = importlib.util.spec_from_file_location('tested_exomol_amendment',
    ROOT/'utility/chemical_database/amend_exomol_metadata.py')
amendment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(amendment)
LABEL = 'Higher energy with complete set of transitions (in cm-1)'


class AmendmentSchemaControls(unittest.TestCase):
    def test_native_label_is_distinct_from_maximum_wavenumber(self):
        # Administrative format tokens only; these are not molecular data.
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch)/'administrative.def'
            path.write_text('12 # Maximum wavenumber (in cm-1)\n8 # '+LABEL+'\n')
            self.assertEqual(amendment.declaration(path, LABEL)['native_value'], '8')

    def test_missing_or_ambiguous_label_is_refused(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch)/'administrative.def'
            for contents in ('12 # Maximum wavenumber (in cm-1)\n',
                             '8 # '+LABEL+'\n8 # '+LABEL+'\n'):
                path.write_text(contents)
                with self.assertRaisesRegex(ValueError, 'one exact native'):
                    amendment.declaration(path, LABEL)


class ActualExoMolMetadataAmendment(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        old = os.environ.get('AMPS_EXOMOL_AMENDMENT_OLD_PACKET')
        evidence = os.environ.get('AMPS_EXOMOL_AMENDMENT_EVIDENCE_PACKET')
        if not old or not evidence:
            raise unittest.SkipTest('explicit old/evidence ExoMol packets are required')
        cls.old, cls.evidence = Path(old).resolve(), Path(evidence).resolve()
        cls.old_sha = 'f777959c00f8bfda6c5297dcb5e0255469b423ef1ca25e994756d78a9380139d'
        cls.evidence_sha = '88d9d5095c687d9109597d67040cff7eb1c1ec08b7e63606b1bd855aff6f15b1'
        cls.source_id = 'molecular-20261119-r4:oh-mythos-completeness-definition'

    def test_real_amendment_preserves_all_science_members(self):
        original = json.loads(self.old.read_text())
        with tempfile.TemporaryDirectory() as scratch:
            result = amendment.amend(self.old, self.old_sha, self.evidence,
                self.evidence_sha, self.source_id, Path(scratch)/'amended')
            published = Path(result['packet'])
            corrected = json.loads(published.read_text())
            self.assertEqual([row['id'] for row in original['sources']],
                             [row['id'] for row in corrected['sources']])
            for before, after in zip(original['sources'], corrected['sources']):
                self.assertEqual(before['sha256'], after['sha256'])
                self.assertEqual((self.old.parent/before['path']).read_bytes(),
                                 (published.parent/after['path']).read_bytes())
            proof = corrected['metadata_amendment']
            self.assertEqual(proof['version_declaration']['native_value'], '20240526')
            self.assertEqual(proof['completeness_declaration']['native_value'], '131898.260953')
            checks = next(row for row in corrected['sources']
                if row['id'] == 'exomol-20261119:oh-mythos-states-bz2:lexical-csv')['quantity_metadata']['source_checks']
            self.assertNotIn('maximum_energy_matches_definition', checks)
            self.assertEqual(checks['maximum_energy_native_cm-1'], '132196.65537')
            self.assertFalse(result['physical_activation'])
        self.assertEqual(hashlib.sha256(self.old.read_bytes()).hexdigest(), self.old_sha)
        self.assertEqual(hashlib.sha256(self.evidence.read_bytes()).hexdigest(), self.evidence_sha)

    def test_changed_original_identity_is_refused_before_copy(self):
        with tempfile.TemporaryDirectory() as scratch:
            output = Path(scratch)/'refused'
            with self.assertRaisesRegex(ValueError, 'identity changed'):
                amendment.amend(self.old, '0'*64, self.evidence,
                    self.evidence_sha, self.source_id, output)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
