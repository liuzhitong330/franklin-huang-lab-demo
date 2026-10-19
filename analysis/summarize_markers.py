"""Transparent donor-level descriptive sensitivity analysis; no significance test.

Run from the repository root:
python analysis/summarize_markers.py --data patient-marker-data.json --output summary.json
All comparisons are within one collection context. A donor is included only
when both target and comparator have the requested minimum number of cells.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from statistics import mean

COLLECTIONS={'RP_Tumor':{'Paired_T','Unpaired_T'},'Biopsy':{'Biopsy'},'Grossly_normal':{'Paired_N'}}
CONTRASTS=(('Club','LE'),('Club','BE'),('LE','BE'))
MIN_CELLS=(10,20,30,50)
THRESHOLDS=(0.1,0.25,0.5)

def comparison(data,collection='RP_Tumor',reference='LE',min_cells=20,threshold=0.25,target='Club'):
    if collection not in COLLECTIONS: raise ValueError('Unknown collection')
    if (target,reference) not in CONTRASTS: raise ValueError('Unsupported target/reference pair')
    if min_cells not in MIN_CELLS: raise ValueError('Unsupported coverage threshold')
    if threshold not in THRESHOLDS: raise ValueError('Unsupported effect threshold')
    selected=[g for g in data['groups'] if g['tissue'] in COLLECTIONS[collection] and g['state'] in {target,reference}]
    grouped=defaultdict(list)
    for g in selected: grouped[(g['donor'],g['state'])].append(g)
    merged={}
    for key,gs in grouped.items():
        n=sum(g['n'] for g in gs)
        merged[key]={'n':n,'mean':[sum(g['mean'][j]*g['n'] for g in gs)/n for j in range(len(data['genes']))],
                     'detect':[sum(g['detect'][j]*g['n'] for g in gs)/n for j in range(len(data['genes']))]}
    donors=sorted(d for d in {g['donor'] for g in selected} if
        merged.get((d,target),{'n':0})['n']>=min_cells and merged.get((d,reference),{'n':0})['n']>=min_cells)
    metrics=[]
    for j,gene in enumerate(data['genes']):
        pairs=[{'donor':d,'targetN':merged[d,target]['n'],'referenceN':merged[d,reference]['n'],
                'targetMean':merged[d,target]['mean'][j],'referenceMean':merged[d,reference]['mean'][j],
                'effect':merged[d,target]['mean'][j]-merged[d,reference]['mean'][j],
                'targetDetection':merged[d,target]['detect'][j],'referenceDetection':merged[d,reference]['detect'][j]} for d in donors]
        effects=[p['effect'] for p in pairs]
        n=len(effects)
        equal=mean(effects) if n else None
        pooled=(sum(p['targetMean']*p['targetN'] for p in pairs)/sum(p['targetN'] for p in pairs)
                -sum(p['referenceMean']*p['referenceN'] for p in pairs)/sum(p['referenceN'] for p in pairs)) if n else None
        loo=[{'omittedDonor':d,'effect':mean(effects[:i]+effects[i+1:])} for i,d in enumerate(donors)] if n>=2 else []
        positive_fraction=sum(e>0 for e in effects)/n if n else None
        passing=bool(n>=4 and equal>=threshold and positive_fraction>=0.75 and all(x['effect']>0 for x in loo))
        metrics.append({'gene':gene,'nDonors':n,'equalDonorEffect':equal,'pooledCellEffect':pooled,
          'positiveDonorFraction':positive_fraction,'looMin':min(x['effect'] for x in loo) if loo else None,
          'looMax':max(x['effect'] for x in loo) if loo else None,'passesReviewRule':passing,'donors':pairs,'leaveOneOut':loo})
    return {'collection':collection,'target':target,'reference':reference,'minCells':min_cells,'threshold':threshold,
            'eligibleDonors':donors,'targetCells':sum(merged[d,target]['n'] for d in donors),
            'referenceCells':sum(merged[d,reference]['n'] for d in donors),
            'metrics':metrics,'passingGenes':[m['gene'] for m in metrics if m['passesReviewRule']]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();data=json.loads(args.data.read_text())
    scenarios=[comparison(data,c,r,n,t,target) for c in COLLECTIONS for target,r in CONTRASTS for n in MIN_CELLS for t in THRESHOLDS]
    result={'method':'Descriptive within-donor target-minus-reference differences; no population inference or independent validation.',
      'reviewRule':'At least 4 eligible donors, mean within-donor effect >= selected threshold, at least 75% donor effects > 0, every leave-one-donor-out mean > 0.',
      'scenarioCount':len(scenarios),'scenarios':scenarios}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
    baseline=comparison(data)
    print(json.dumps({'scenarioCount':len(scenarios),'defaultEligibleDonors':baseline['eligibleDonors'],
      'defaultPassingGenes':baseline['passingGenes'],'defaultTargetCells':baseline['targetCells'],
      'defaultReferenceCells':baseline['referenceCells']},indent=2))

if __name__=='__main__': main()
