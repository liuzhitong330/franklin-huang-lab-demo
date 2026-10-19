"""Standard-library regression tests for deposited-data summaries and decisions.

These tests do not replace independent raw-count validation or donor-level
validation in a new cohort. Run: python analysis/test_analysis.py
"""
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import unittest
from summarize_markers import comparison, COLLECTIONS, CONTRASTS, MIN_CELLS, THRESHOLDS

ROOT=Path(__file__).resolve().parents[1]
DATA=json.loads((ROOT/'patient-marker-data.json').read_text())

class AnalysisTests(unittest.TestCase):
    def test_01_source_and_release_checksum(self):
        manifest=json.loads((ROOT/'analysis/source_manifest.json').read_text())
        self.assertEqual(DATA['meta']['sourceSha256'],manifest['deposited_object']['sha256'])
        self.assertEqual(hashlib.sha256((ROOT/'patient-marker-data.json').read_bytes()).hexdigest(),manifest['derived_json']['sha256'])
        self.assertEqual(DATA['meta']['normalizationQA']['maxAbsoluteCountTotalDifferenceVsNCountRNA'],0)
        self.assertLess(DATA['meta']['normalizationQA']['maxAbsoluteErrorVsDepositedData'],1e-10)

    def test_02_dimensions_and_coverage(self):
        self.assertEqual(len(DATA['genes']),20);self.assertEqual(len(set(DATA['genes'])),20)
        self.assertEqual(len(DATA['groups']),73);self.assertEqual(sum(g['n'] for g in DATA['groups']),13322)
        self.assertEqual(len({g['donor'] for g in DATA['groups']}),11)
        self.assertEqual(len({s for g in DATA['groups'] for s in g['sourceSamples']}),15)
        totals=Counter()
        for g in DATA['groups']: totals[g['state']]+=g['n']
        self.assertEqual(dict(totals),DATA['meta']['stateCounts'])

    def test_03_finite_expression_and_detection(self):
        keys=[]
        for g in DATA['groups']:
            keys.append((g['donor'],g['tissue'],g['state']))
            self.assertEqual(len(g['mean']),20);self.assertEqual(len(g['detect']),20)
            self.assertTrue(all(math.isfinite(x) and x>=0 for x in g['mean']))
            for f in g['detect']:
                self.assertTrue(math.isfinite(f) and 0<=f<=1)
                self.assertAlmostEqual(f*g['n'],round(f*g['n']),places=9)
        self.assertEqual(len(keys),len(set(keys)))

    def test_04_donor_and_tissue_provenance(self):
        for g in DATA['groups']:
            self.assertEqual(len(g['sourceSamples']),1)
            sample=g['sourceSamples'][0]
            self.assertEqual(DATA['meta']['donorMap'][sample],g['donor'])
            if sample.endswith('_N'): self.assertEqual(g['tissue'],'Paired_N')
            if sample.endswith('_T'): self.assertEqual(g['tissue'],'Paired_T')
        self.assertEqual(DATA['meta']['donorMap']['PR5261_N'],'PR5261')
        self.assertEqual(DATA['meta']['donorMap']['PR5254_N'],'PR5254')

    def test_05_default_eligibility(self):
        x=comparison(DATA)
        self.assertEqual(x['eligibleDonors'],['PR5196','PR5249','PR5251','PR5269'])
        self.assertEqual((x['targetCells'],x['referenceCells']),(268,1607))
        self.assertEqual(x['passingGenes'],['KRT15','SERPINB1','PIGR','MMP7','CP','CD74','LTF','IL8'])

    def test_06_checked_default_metrics(self):
        metrics={m['gene']:m for m in comparison(DATA)['metrics']}
        for gene,expected in {'MMP7':(2.6413943794153125,2.5770627265898844,2.364718880839915),
                              'PIGR':(2.3843852610210163,2.5956945068625012,2.2630429280114046),
                              'SCGB1A1':(0.23583460960116562,0.1639958016090296,0.12216399643206736)}.items():
            for key,value in zip(('equalDonorEffect','pooledCellEffect','looMin'),expected):
                self.assertAlmostEqual(metrics[gene][key],value,places=12)
        self.assertFalse(metrics['SCGB1A1']['passesReviewRule'])
        self.assertEqual(metrics['SCGB1A1']['positiveDonorFraction'],0.5)

    def test_07_threshold_monotonicity(self):
        for c in COLLECTIONS:
            for target,reference in CONTRASTS:
                for n in MIN_CELLS:
                    sets=[set(comparison(DATA,c,reference,n,t,target)['passingGenes']) for t in THRESHOLDS]
                    self.assertTrue(sets[2]<=sets[1]<=sets[0])

    def test_08_insufficient_donors_and_loo(self):
        for c in COLLECTIONS:
            for target,reference in CONTRASTS:
                for n in MIN_CELLS:
                    x=comparison(DATA,c,reference,n,0.25,target)
                    for m in x['metrics']:
                        if m['nDonors']<4: self.assertFalse(m['passesReviewRule'])
                        if m['nDonors']<2: self.assertIsNone(m['looMin'])
                        if m['passesReviewRule']:
                            self.assertGreater(m['looMin'],0)
                            self.assertGreaterEqual(m['positiveDonorFraction'],0.75)

    def test_09_input_guards(self):
        for kwargs in ({'collection':'All'},{'reference':'Club'},{'target':'ERGneg_Tumor'},
                       {'min_cells':0},{'threshold':float('nan')}):
            with self.assertRaises(ValueError): comparison(DATA,**kwargs)

    def test_10_all_108_configurations(self):
        scenarios=[comparison(DATA,c,r,n,t,target) for c in COLLECTIONS for target,r in CONTRASTS for n in MIN_CELLS for t in THRESHOLDS]
        self.assertEqual(len(scenarios),108)
        self.assertEqual(sum(len(s['metrics']) for s in scenarios),2160)
        # Serialization rejects any accidental NaN/Infinity in no-donor cases.
        json.dumps(scenarios,allow_nan=False)

if __name__=='__main__': unittest.main(verbosity=2)
