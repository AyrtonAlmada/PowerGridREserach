#!/usr/bin/env python3
"""Record a LOCAL ParaEMT checkout without modifying it. Python standard library only.

Run in the actual ParaEMT Python/Conda environment. This records today's checkout;
it cannot establish which revision generated older trajectories. The author must
associate this record with the relevant run manifests. Inspect diffs for private
paths/credentials before publishing them.
"""
import argparse, datetime, hashlib, json, os, pathlib, platform, shutil, subprocess, sys

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()

def command(args,cwd):
    p=subprocess.run(args,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
    return p.returncode,p.stdout,p.stderr

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('checkout',type=pathlib.Path)
    ap.add_argument('output',type=pathlib.Path)
    ap.add_argument('--run-label',required=True,help='E.g. calibration-rerun-2026-09-20; not a historical inference.')
    ap.add_argument('--extra-file',action='append',default=[],help='Additional local snapshot/config path relative to checkout.')
    args=ap.parse_args();root=args.checkout.resolve();out=args.output.resolve()
    if not root.is_dir():ap.error('Checkout directory does not exist.')
    if out.exists():ap.error('Output already exists; do not overwrite provenance records.')
    out.mkdir(parents=True)
    meta={'capture_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'run_label':args.run_label,'upstream_url':'https://github.com/NatLabRockies/ParaEMT_public',
          'python':sys.version,'python_executable_name':pathlib.Path(sys.executable).name,
          'platform':platform.platform(),'processor':platform.processor(),'logical_cpus':os.cpu_count(),
          'commit_status':'unknown','author_confirmed_historical_run':False}
    if shutil.which('git'):
        rc,head,err=command(['git','rev-parse','HEAD'],root)
        if rc==0:
            sha=head.decode().strip();meta['commit']=sha;meta['commit_status']='captured_local_checkout'
            meta['permalink']='https://github.com/NatLabRockies/ParaEMT_public/tree/'+sha
            for name,cmd in [('status.txt',['git','status','--porcelain=v1']),
                             ('tracked_changes.patch',['git','diff','--binary','HEAD']),
                             ('untracked_files.txt',['git','ls-files','--others','--exclude-standard']),
                             ('remotes.txt',['git','remote','-v'])]:
                rc,b,err=command(cmd,root);(out/name).write_bytes(b)
            meta['dirty']=bool((out/'status.txt').read_text().strip())
            rc,b,err=command(['git','show',sha+':LICENSE.md'],root)
            if rc==0:(out/'LICENSE_paraemt_at_commit.md').write_bytes(b)
        else: meta['git_message']=err.decode(errors='replace').strip()
    # Hash Python source, environment specifications, and named case inputs.
    targets=set(root.glob('*.py'))
    for p in root.rglob('*'):
        if not p.is_file() or '.git' in p.parts:continue
        if p.name.lower() in {'240bus.xls','wecc240.raw','pfd_240_1_1.json','requirements.txt','environment.yml','license.md'}:
            targets.add(p)
    for item in args.extra_file:
        p=(root/item).resolve()
        if not p.is_relative_to(root) or not p.is_file():ap.error('extra-file must name a file inside the checkout: '+item)
        targets.add(p)
        destination=out/'extra_files'/p.relative_to(root)
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,destination)
    meta['files']=[{'path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':digest(p)}
                   for p in sorted(targets)]
    rc,freeze,err=command([sys.executable,'-m','pip','freeze'],root)
    if rc==0:(out/'pip-freeze.txt').write_bytes(freeze)
    else:meta['pip_freeze_error']=err.decode(errors='replace')
    (out/'paraemt_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Saved checkout provenance:',out)
    print('Commit:',meta.get('commit','UNKNOWN -- no Git history'))
    print('Confirm separately that this checkout and patch generated the archived runs.')
if __name__=='__main__':main()
