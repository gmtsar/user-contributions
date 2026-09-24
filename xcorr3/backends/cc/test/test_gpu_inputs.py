#!/usr/bin/env python3
"""Exercise production input rejection before CUDA and valid SLC path handling.

Requires a built candidate; valid-path checks additionally require a GPU and
the 1536x1024 known-shift (3,2) fixtures from the Dingri validation suite.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',required=True,type=Path)
p.add_argument('--fixtures',required=True,type=Path)
p.add_argument('--work',required=True,type=Path)
a=p.parse_args();binary=a.binary.resolve();root=a.work.resolve();root.mkdir(parents=True,exist_ok=True)
baseline={name:(a.fixtures/name).read_text() for name in ('master.PRM','aligned.PRM')}
records=[]
def run(name,options=(),master=None,secondary=None,valid=False,error_token=None):
    w=root/name;w.mkdir(exist_ok=True)
    texts=dict(baseline)
    if master is not None:texts['master.PRM']=master
    if secondary is not None:texts['aligned.PRM']=secondary
    for path,text in texts.items():(w/path).write_text(text)
    for path in ('master.SLC','aligned.SLC'):
        link=w/path
        if not link.exists():link.symlink_to((a.fixtures/path).resolve())
    out=w/'freq_xcorr.dat';sentinel=b'previous complete output\n';out.write_bytes(sentinel)
    env=dict(os.environ)
    if not valid:env['CUDA_VISIBLE_DEVICES']=''
    cmd=[str(binary),'master.PRM','aligned.PRM','-nx','7','-ny','5','-xsearch','32','-ysearch','32',*options]
    done=subprocess.run(cmd,cwd=w,env=env,capture_output=True,text=True,timeout=60)
    if valid:
        rows=[[float(v) for v in line.split()] for line in out.read_text().splitlines()]
        passed=done.returncode==0 and len(rows)==35 and all(len(row)==6 and all(math.isfinite(x) for x in row) and row[1]==3 and row[3]==2 for row in rows)
    else:
        passed=bool(done.returncode>0 and out.read_bytes()==sentinel and error_token and error_token in done.stderr and 'cudaGetDeviceCount' not in done.stderr and 'No CUDA' not in done.stderr)
    records.append(dict(name=name,pass_=passed,expected_error=error_token,command=cmd,returncode=done.returncode,stderr=done.stderr))
run('unknown', ['-typo'],error_token='-typo')
run('missing_value',['-ny'],error_token='-ny')
run('wrapped_integer',['-nx','4294967297'],error_token='-nx requires an integer')
run('empty_key',master='= value\n',error_token='empty key')
run('missing_path',master='num_rng_bins = 1536\n',error_token='SLC_file')
run('range_shift_bounds',secondary=baseline['aligned.PRM']+'rshift = 2147483647\n',error_token='range sampling grid')
run('azimuth_shift_bounds',secondary=baseline['aligned.PRM']+'ashift = -2147483648\n',error_token='azimuth sampling grid')
run('stretch_bounds',secondary=baseline['aligned.PRM']+'PRF = 1e100\n',error_token='azimuth sampling grid')
run('invalid_stretch',master=baseline['master.PRM']+'PRF = 1e-100\n',secondary=baseline['aligned.PRM']+'PRF = 1e300\n',error_token='invalid azimuth stretch')
run('cancelled_stretch_overflow',['-ny','1'],master=baseline['master.PRM']+'PRF = 544\n',secondary=baseline['aligned.PRM']+'PRF = 2147484448\nashift = -2147483648\n',error_token='PRF stretch exceeds the azimuth index range')
run('normal',['-nointerp'],valid=True)
# Internal spaces and '=' must survive parsing through actual binary I/O.
w=root/'paths_with_spaces';w.mkdir(exist_ok=True)
for src,dst in [('master.SLC','master = image.SLC'),('aligned.SLC','aligned image.SLC')]:
    if not (w/dst).exists():(w/dst).symlink_to((a.fixtures/src).resolve())
run('paths_with_spaces',['-nointerp'],master=baseline['master.PRM'].replace('master.SLC','master = image.SLC'),secondary=baseline['aligned.PRM'].replace('aligned.SLC','aligned image.SLC'),valid=True)
report={'pass':all(r['pass_'] for r in records),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'cases':records}
(root/'results.json').write_text(json.dumps(report,indent=2))
print(f"{sum(r['pass_'] for r in records)}/{len(records)} production input cases passed")
sys.exit(0 if report['pass'] else 1)
