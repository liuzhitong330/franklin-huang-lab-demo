"""Extract patient/state summaries from the deposited 2022 Seurat object.

This is a processed-count reanalysis, not a rerun of alignment, QC, clustering,
or cell annotation. Only metadata, active identities and RNA counts/data are
converted; stored R commands are never evaluated. Requires rdata 1.1.0,
numpy, scipy, pandas (rdata also installs xarray).
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import rdata
from scipy.sparse import csc_matrix

SOURCE_SHA='c0cbde35277b2f250ee13160afe9d163b62c720dca8c9f43f121a47e98b29430'
REPO='https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples'
COMMIT='71dff9b7dc3d096b6ea2da1ed2cc5bae22876982'
GENES=['KRT5','KRT15','KRT17','TP63','KLK3','KRT8','KRT18','KRT13',
       'APOBEC3A','CSTB','LYPD3','SERPINB1','KLK2','PIGR','MMP7','CP',
       'SCGB1A1','CD74','LTF','IL8']
PAIRED_DONORS=['PR5249','PR5251','PR5254','PR5261']
DONOR_MAP={f'{d}_{s}':d for d in PAIRED_DONORS for s in ['N','T']}
DONOR_MAP.update({s:s for s in ['AUG_PB1','MAY_PB1','MAY_PB2','PR5186','PR5196','PR5199','PR5269']})

def deref(o):
    while o.info.type==rdata.parser.RObjectType.REF: o=o.referenced_object
    return o

def slots(o):
    node=deref(o).attributes
    out={}
    while node is not None and node.info.type==rdata.parser.RObjectType.LIST:
        tag=deref(node.tag)
        if tag.info.type==rdata.parser.RObjectType.SYM: tag=deref(tag.value)
        out[tag.value.decode('utf-8')]=node.value[0]
        node=deref(node.value[1])
    return out

def matrix(x):
    assert len(x.p)==int(x.Dim[1])+1
    return csc_matrix((x.x,x.i,x.p),shape=tuple(x.Dim)), list(map(str,x.Dimnames[0])),list(map(str,x.Dimnames[1]))

def extract(path,output_path):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==SOURCE_SHA,'Source SHA256 mismatch'
    parsed=rdata.parser.parse_file(path)
    nodes=slots(parsed.object)
    metadata=rdata.conversion.convert(nodes['meta.data'])
    ident=rdata.conversion.convert(nodes['active.ident'])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Missing constructor for R class')
        assays=rdata.conversion.convert(nodes['assays'])
    assert set(assays)=={'RNA'}
    counts,features,cells=matrix(assays['RNA'].counts)
    stored,stored_features,stored_cells=matrix(assays['RNA'].data)
    assert list(map(str,metadata.index))==cells==stored_cells
    assert features==stored_features
    assert list(map(str,metadata['ID']))==list(map(str,ident))
    assert counts.shape==(21877,13322)
    assert np.all(np.isfinite(counts.data)) and np.all(counts.data>=0)
    assert np.all(counts.data==np.floor(counts.data))
    totals=np.asarray(counts.sum(axis=0)).ravel()
    assert np.all(totals>0)
    meta_totals=np.asarray(metadata['nCount_RNA'],dtype=float)
    original_genes=GENES.copy()
    aliases={}
    if 'IL8' not in features and 'CXCL8' in features: aliases['IL8']='CXCL8'
    resolved=[aliases.get(g,g) for g in GENES]
    missing=[g for g in resolved if g not in features]
    assert not missing,missing
    selected=counts[[features.index(g) for g in resolved],:].toarray()
    log_expression=np.log1p(selected/totals[None,:]*10000.)
    stored_selected=stored[[features.index(g) for g in resolved],:].toarray()
    max_stored_error=float(np.max(np.abs(stored_selected-log_expression)))
    # Stored RNA data is checked, not assumed to be identically normalized.
    assert max_stored_error<1e-10,('Deposited data differs from recomputed LogNormalize',max_stored_error)
    metadata=metadata.copy()
    metadata['donor']=[DONOR_MAP[str(x)] for x in metadata['orig.ident']]
    metadata['state']=metadata['ID'].astype(str)
    metadata['tissue']=metadata['sampletype'].astype(str)
    groups=[]
    for (donor,tissue,state),indices in metadata.groupby(['donor','tissue','state'],observed=True).indices.items():
        indices=np.asarray(indices,dtype=int)
        groups.append({'donor':str(donor),'tissue':str(tissue),'state':str(state),
                       'sourceSamples':sorted(set(map(str,metadata.iloc[indices]['orig.ident']))),
                       'n':int(len(indices)),
                       'mean':np.mean(log_expression[:,indices],axis=1).tolist(),
                       'detect':np.mean(selected[:,indices]>0,axis=1).tolist()})
    coverage=[]
    for (sample,donor,tissue,state),indices in metadata.groupby(['orig.ident','donor','tissue','state'],observed=True).indices.items():
        coverage.append({'sample':str(sample),'donor':str(donor),'tissue':str(tissue),'state':str(state),'n':int(len(indices))})
    source_meta={'cells':int(counts.shape[1]),'rnaFeatures':int(counts.shape[0]),
      'donors':int(metadata['donor'].nunique()),'sourceSamples':int(metadata['orig.ident'].nunique()),
      'sourceFile':path.name,'sourceSha256':SOURCE_SHA,'repository':REPO,'commit':COMMIT,
      'paper':'https://doi.org/10.1038/s41467-021-27322-4','geo':'GSE176031',
      'archiveUrl':'https://www.dropbox.com/sh/ga770h1ik75l1o2/AADfmpEsHp8OLTtfkZaCxPd8a?dl=1',
      'normalization':'natural log(1 + 10000 * gene count / total counts across all 21877 deposited RNA features)',
      'expressionUnit':'mean log1p counts per 10000 deposited RNA counts',
      'differenceUnit':'difference of mean log1p-normalized expression; not log fold change',
      'detectionUnit':'fraction of cells with deposited count > 0',
      'sourceLabelPolicy':'Original author cell-state labels retained; no new malignancy diagnosis or lineage inference.',
      'genePanelSource':'Fixed 20-gene union from Seqwell_Combined_E.R lines 232/252 and 836-847 at the pinned commit; selected before derived results, not an unbiased transcriptome-wide target screen.',
      'originalGeneSymbols':original_genes,'aliases':aliases,'donorMap':DONOR_MAP,
      'il8FeatureAudit':{'IL8_present':'IL8' in features,'CXCL8_present':'CXCL8' in features,'policy':'Use deposited IL8 only; do not silently combine feature rows.'},
      'tissues':sorted(map(str,metadata['tissue'].unique())),
      'stateCounts':{str(k):int(v) for k,v in metadata['state'].value_counts().items()},
      'tissueCounts':{str(k):int(v) for k,v in metadata['tissue'].value_counts().items()},
      'normalizationQA':{'maxAbsoluteErrorVsDepositedData':max_stored_error,
        'maxAbsoluteCountTotalDifferenceVsNCountRNA':float(np.max(np.abs(totals-meta_totals))),
        'integerNonnegativeCounts':True,'cellOrderVerified':True,'identitiesMatchMetadata':True}}
    output={'meta':source_meta,'genes':resolved,'groups':groups,'coverage':coverage}
    output_path.parent.mkdir(parents=True,exist_ok=True)
    output_path.write_text(json.dumps(output,ensure_ascii=False,separators=(',',':'))+'\n')
    print(json.dumps({'cells':len(cells),'features':len(features),'donors':source_meta['donors'],'groups':len(groups),'genes':len(resolved),'normalizationQA':source_meta['normalizationQA'],'il8FeatureAudit':source_meta['il8FeatureAudit'],'output':str(output_path)},indent=2))
    return output

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--rds',type=Path,required=True,help='Path to the checksum-matched deposited dge_E.rds')
    ap.add_argument('--output',type=Path,required=True,help='Output patient-marker-data.json')
    args=ap.parse_args()
    extract(args.rds,args.output)

