#!/usr/bin/env python3
"""Production CUDA output failure tests; requires known-shift SLC fixtures."""
import argparse,json,os,pathlib,resource,subprocess,sys
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=pathlib.Path,required=True)
p.add_argument('--fixtures',type=pathlib.Path,required=True)
p.add_argument('--work',type=pathlib.Path,required=True)
a=p.parse_args();r=a.work.resolve();r.mkdir(parents=True,exist_ok=True)
records=[]
for name,extra,token in [('normal',[] ,None),('write_limit',[],'correlation output'),('missing_trans',['-geocode'],'trans.dat'),('directory',[],'not a regular file'),('device_symlink',[],'not a regular file'),('unwritable',[],'staging directory')]:
    w=r/name;w.mkdir(exist_ok=True)
    for n in ('master.PRM','aligned.PRM','master.SLC','aligned.SLC'):
        (w/n).symlink_to((a.fixtures/n).resolve())
    out=w/'freq_xcorr.dat';old=b'previous complete result\n';out.write_bytes(old)
    if name=='directory':out.unlink();out.mkdir()
    if name=='device_symlink':out.unlink();out.symlink_to('/dev/full')
    if name=='unwritable':w.chmod(0o555)
    def restrict():resource.setrlimit(resource.RLIMIT_FSIZE,(128,128))
    env=dict(os.environ,CUDA_CACHE_DISABLE='1')
    try:
        done=subprocess.run([str(a.binary.resolve()),'master.PRM','aligned.PRM','-nx','7','-ny','5','-xsearch','32','-ysearch','32','-nointerp',*extra],cwd=w,env=env,capture_output=True,text=True,timeout=90,preexec_fn=restrict if name=='write_limit' else None)
    finally:
        if name=='unwritable':w.chmod(0o755)
    passed=not list(w.glob('.xcorr_cc-*'))
    if token:
        if name=='write_limit':token='flush' if 'flush' in done.stderr else token
        preserved=out.is_dir() if name=='directory' else out.is_symlink() if name=='device_symlink' else out.read_bytes()==old
        passed=passed and done.returncode>0 and token in done.stderr and preserved
    else:
        rows=[line.split() for line in out.read_text().splitlines()]
        passed=passed and done.returncode==0 and len(rows)==35 and all(float(row[1])==3 and float(row[3])==2 for row in rows)
    records.append(dict(name=name,pass_=bool(passed),returncode=done.returncode,stderr=done.stderr))
(r/'results.json').write_text(json.dumps({'pass':all(x['pass_'] for x in records),'cases':records},indent=2))
print(f"{sum(x['pass_'] for x in records)}/{len(records)} production output cases passed")
for x in records:
    if not x['pass_']:print(x)
sys.exit(0 if all(x['pass_'] for x in records) else 1)
