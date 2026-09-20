#!/usr/bin/env python3
"""Check and document supplied input bytes, row mappings and graph connectivity."""
import csv,json,hashlib,pathlib,collections
ROOT=pathlib.Path(__file__).resolve().parents[1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def audit():
    directory=ROOT/'data/wecc240/processed'
    branches=list(csv.DictReader((directory/'branch240E3.csv').open()))
    buses=list(csv.DictReader((directory/'bus240E3.csv').open()))
    ids=[int(r['index']) for r in buses];pos={v:k for k,v in enumerate(ids)}
    assert len(pos)==len(ids),'Duplicate bus IDs'
    adj=[set() for _ in ids]
    edges=[]
    for r in branches:
        i,j=pos[int(r['f_bus'])],pos[int(r['t_bus'])]
        assert float(r['br_x'])>0
        adj[i].add(j);adj[j].add(i);edges.append((i,j))
    unseen=set(range(len(ids)));components=[]
    while unseen:
        stack=[min(unseen)];component=[];unseen.remove(stack[0])
        while stack:
            v=stack.pop();component.append(v)
            for w in sorted(adj[v]):
                if w in unseen:unseen.remove(w);stack.append(w)
        components.append(sorted(component))
    p=[float(r['p']) for r in buses]
    upstream=json.loads((ROOT/'data/wecc240/upstream/pfd_240_1_1.json').read_text())
    report={'processed_buses':len(buses),'processed_branches':len(branches),
            'connected_components':len(components),
            'component_sizes':[len(c) for c in components],
            'component_injection_sums':[sum(p[i] for i in c) for c in components],
            'isolated_buses':[ids[i] for i,a in enumerate(adj) if not a],
            'inertia_values':sorted({float(r['m']) for r in buses}),
            'damping_values':sorted({float(r['d']) for r in buses}),
            'limit_values':sorted({float(r['rate_a']) for r in branches}),
            'parallel_same_direction':[list(k)+[v] for k,v in collections.Counter(edges).items() if v>1],
            'json_buses':len(upstream['bus_num']),'json_lines':len(upstream['line_from']),
            'json_transformers':len(upstream['xfmr_from']),
            'json_bus_ids_match_processed':set(ids)==set(upstream['bus_num']),
            'raw_to_processed_transformation_script_supplied':False}
    out=ROOT/'data/wecc240/metadata';out.mkdir(parents=True,exist_ok=True)
    (out/'input_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    with (out/'bus_map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['state_index','bus_id','component'])
        membership={v:c+1 for c,nodes in enumerate(components) for v in nodes}
        for k,id in enumerate(ids):w.writerow([k+1,id,membership[k]])
    with (out/'branch_map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['line_row','original_index','from_bus','to_bus','From','To','state_label'])
        for k,r in enumerate(branches):
            i,j=pos[int(r['f_bus'])]+1,pos[int(r['t_bus'])]+1
            w.writerow([k+1,r['index'],r['f_bus'],r['t_bus'],i,j,f'X{i}#{j}'])
    inputs=list((ROOT/'data/wecc240/upstream').glob('*'))+list(directory.glob('*.csv'))+list((ROOT/'legacy/source').glob('*'))
    with (ROOT/'provenance/input_checksums.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['path','bytes','sha256'])
        for item in sorted(inputs):
            if item.is_file():w.writerow([item.relative_to(ROOT).as_posix(),item.stat().st_size,digest(item)])
    return report
if __name__=='__main__':print(json.dumps(audit(),indent=2))
