#!/usr/bin/env python3
"""Index and hash AUTHOR-PROVIDED trajectory files; never infer their generator.

Requires an explicit metadata CSV: source,scenario_id,line_row,lambda_code,
Top_s,Tev_s,file,split,paraemt_commit. Files are relative to --data-root.
This prevents a surrogate file with a 'ParaEMT' prefix becoming an EMT reference.
"""
import argparse,csv,hashlib,pathlib,re,json

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata',required=True,type=pathlib.Path)
    ap.add_argument('--data-root',required=True,type=pathlib.Path)
    ap.add_argument('--output',required=True,type=pathlib.Path)
    ap.add_argument('--require-1150',action='store_true')
    ap.add_argument('--allow-unknown-commit',action='store_true',
        help='Allow an explicitly unknown historical commit only with source_archive_sha256 for the saved source bundle.')
    a=ap.parse_args();root=a.data_root.resolve()
    rows=list(csv.DictReader(a.metadata.open(newline='',encoding='utf-8-sig')))
    if not rows:ap.error('Metadata must contain actual trajectory records.')
    required={'source','scenario_id','line_row','lambda_code','Top_s','Tev_s','file','split','paraemt_commit'}
    if not required<=set(rows[0]):ap.error('Missing columns: '+','.join(required-set(rows[0])))
    ids=set();combos=set();paths=set();out=[]
    for r in rows:
        if r['source'] not in ('ParaEMT','Surrogate'):ap.error('Explicit source must be ParaEMT or Surrogate.')
        key=r['source'],r['scenario_id']
        if key in ids:ap.error('Duplicate source/scenario_id: '+str(key))
        ids.add(key)
        line,lc=int(r['line_row']),int(r['lambda_code']);top,tev=float(r['Top_s']),float(r['Tev_s'])
        if line<1 or lc<1:ap.error('line_row and lambda_code must be positive.')
        if not(0<=top<=tev):ap.error('Invalid time convention for '+r['scenario_id'])
        if r['split'] not in ('train','validation','test','not_split'):ap.error('Set an explicit split status.')
        if r['source']=='ParaEMT' and not re.fullmatch('[0-9a-fA-F]{40}',r['paraemt_commit']):
            unknown=r['paraemt_commit'].strip().lower() in ('','unknown','unavailable')
            hashed=bool(re.fullmatch('[0-9a-fA-F]{64}',r.get('source_archive_sha256','')))
            if not(a.allow_unknown_commit and unknown and hashed):
                ap.error('Provide the actual full ParaEMT commit, or explicit unknown status with a saved source-bundle hash and --allow-unknown-commit. Never guess a SHA.')
        path=(root/r['file']).resolve()
        if not path.is_relative_to(root) or not path.is_file():ap.error('Missing/out-of-root data file: '+r['file'])
        if path in paths:ap.error('One trajectory file is assigned to multiple scenarios: '+r['file'])
        paths.add(path)
        record=dict(r);record.update(bytes=path.stat().st_size,sha256=sha(path));out.append(record)
        if r['source']=='ParaEMT':combos.add((line,lc,top,tev))
    if a.require_1150:
        emt=[r for r in rows if r['source']=='ParaEMT']
        if len(emt)!=1150 or len(combos)!=1150:ap.error('Expected exactly 1,150 unique EMT references.')
        lines={int(r['line_row']) for r in emt};lambdas={int(r['lambda_code']) for r in emt}
        if len(lines)!=25 or lambdas!=set(range(75,121)):ap.error('Expected 25 lines x lambda codes 75:120.')
        if len({(r['Top_s'],r['Tev_s']) for r in emt})!=1:ap.error('The base 1,150 set must use one fixed pair of times.')
        if {(l,c) for l,c,_,_ in combos}!={(l,c) for l in lines for c in lambdas}:ap.error('Incomplete 25 x 46 Cartesian grid.')
    if a.output.exists():ap.error('Output exists; do not overwrite a published index.')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0]));writer.writeheader();writer.writerows(out)
    print('Indexed',len(out),'files with source declarations and SHA-256 checksums.')
if __name__=='__main__':main()
