#!/usr/bin/env python3
"""Mode and failure regressions. Requires NumPy, GMT and built executables.

python3 test/test_modes.py --original /path/to/xcorr --candidate ./xcorr_mt \
    --work /tmp/xcorr-modes
All generated inputs/outputs stay in --work. Run separately from benchmarks.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--original',required=True,type=Path)
p.add_argument('--candidate',required=True,type=Path)
p.add_argument('--work',required=True,type=Path)
args=p.parse_args()
orig=args.original.resolve();mt=args.candidate.resolve();root=args.work.resolve()
root.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OPENBLAS_NUM_THREADS='1');env.pop('OMP_NUM_THREADS',None)
rng=np.random.default_rng(93851)
nx,ny=1536,1024
m=rng.integers(-2000,2000,(ny,nx,2),dtype=np.int16)
s=np.zeros_like(m);s[2:,3:]=m[:-2,:-3]
fixtures=root/'fixtures';fixtures.mkdir(exist_ok=True)
for name,a in [('master',m),('aligned',s)]:
    a.astype('<i2').tofile(fixtures/(name+'.SLC'))
    amp=np.sqrt((a.astype(float)**2).sum(axis=2)).astype('<f4')
    amp.tofile(fixtures/(name+'.flt'))
    subprocess.run(['gmt','xyz2grd',name+'.flt','-ZTLf',f'-R0/{nx-1}/0/{ny-1}','-I1','-G'+name+'.grd'],cwd=fixtures,check=True,env=env)
    for suffix,slc in [('.PRM','.SLC'),('F.PRM','.flt')]:
        (fixtures/(name+suffix)).write_text(f'SLC_file = {name}{slc}\nnum_rng_bins = {nx}\nnum_patches = 1\nnum_valid_az = {ny}\nrshift = 0\nashift = 0\nPRF = 2000\n')
base=['-nx','7','-ny','5','-xsearch','32','-ysearch','32']
results=[]
def execute(tag,exe,options,inputs=('master.PRM','aligned.PRM'),extra_env=None,limit=None,sentinel=False,directory_output=False):
    work=root/tag;work.mkdir(exist_ok=True)
    for f in fixtures.iterdir():
        dst=work/f.name
        if not dst.exists():dst.symlink_to(f)
    outfile='time_xcorr_Gatelli.dat' if '-time4' in options else 'time_xcorr.dat' if '-time' in options else 'freq_xcorr.dat'
    out=work/outfile
    if sentinel:out.write_bytes(b'previous complete result\n')
    if directory_output:out.mkdir(exist_ok=True)
    def restrict():resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
    cmd=[str(exe),*inputs,*base,*options]
    proc=subprocess.run(cmd,cwd=work,env={**env,**(extra_env or {})},capture_output=True,preexec_fn=restrict if limit else None)
    (work/'stdout.log').write_bytes(proc.stdout);(work/'stderr.log').write_bytes(proc.stderr)
    result={'tag':tag,'command':cmd,'returncode':proc.returncode}
    if proc.returncode==0:
        a=np.loadtxt(out)
        expected={'-nx':7,'-ny':5}
        for i,opt in enumerate(options[:-1]):
            if opt in expected:expected[opt]=int(options[i+1])
        a=np.atleast_2d(a)
        assert a.shape==(expected['-nx']*expected['-ny'],5),tag
        if exe==orig and '-nointerp' in options:
            # Undefined stock fractions can be finite garbage or NaN/Inf.
            # Only coordinates/correlation are a usable reference in this mode.
            assert np.isfinite(a[:,[0,2,4]]).all(),tag
            result['stock_nonfinite_offset_rows']=int((~np.isfinite(a[:,[1,3]])).any(axis=1).sum())
        else:
            assert np.isfinite(a).all(),tag
        result['sha256']=hashlib.sha256(out.read_bytes()).hexdigest()
    assert not list(work.glob('.xcorr_mt-*')),(tag,'temporary files left behind')
    if sentinel:assert out.read_bytes()==b'previous complete result\n',tag
    results.append(result)
    (root/'results.json').write_text(json.dumps(results,indent=2))
    return result

for mode,options,inputs in [
    ('frequency',[],('master.PRM','aligned.PRM')),
    ('time',['-time'],('master.PRM','aligned.PRM')),
    ('time4',['-time4'],('master.PRM','aligned.PRM')),
    ('real',['-real'],('masterF.PRM','alignedF.PRM')),
    ('grid',[],('master.grd','aligned.grd')),
    ('nointerp',['-nointerp'],('master.PRM','aligned.PRM')),
    ('norange',['-norange'],('master.PRM','aligned.PRM')),
    ('interp8',['-interp','8'],('master.PRM','aligned.PRM')),
    ('range4',['-range_interp','4'],('master.PRM','aligned.PRM')),
    ('noshift',['-noshift'],('master.PRM','aligned.PRM')),
]:
    ref=execute(mode+'_original',orig,options,inputs)
    assert ref['returncode']==0,ref
    for workers in [1,3,'auto']:
        opts=options+([] if workers=='auto' else ['-nproc',str(workers)])
        actual=execute(mode+'_mt'+str(workers),mt,opts,inputs)
        assert actual['returncode']==0,actual
        if mode=='nointerp':
            # Stock print_results reads uninitialized xfrac/yfrac when
            # interpolation is disabled. Preserve that output as evidence;
            # the known integer translation is the independent oracle here.
            data=np.loadtxt(root/actual['tag']/'freq_xcorr.dat')
            reference=np.loadtxt(root/ref['tag']/'freq_xcorr.dat')
            assert np.array_equal(data[:,[0,2,4]],reference[:,[0,2,4]])
            assert np.array_equal(data[:,1],np.full(35,3.0))
            assert np.array_equal(data[:,3],np.full(35,2.0))
            actual['oracle']='known translation dx=3, dy=2; stock nointerp has uninitialized fractions'
        else:
            assert actual['sha256']==ref['sha256'],(ref,actual)
        print('PASS',actual['tag'],flush=True)

for name,options,extra in [('cap',['-nproc','5000'],{}),('zero',['-nproc','0'],{}),('text',['-nproc','abc'],{}),('overflow',['-nproc','999999999999999999999999'],{}),('env',[],{'OMP_NUM_THREADS':'3'}),('precedence',['-nproc','2'],{'OMP_NUM_THREADS':'3'})]:
    result=execute(name,mt,options,extra_env=extra)
    assert result['returncode']==0 and result['sha256']==results[0]['sha256'],result
    if name in ['env','precedence']:
        wanted=3 if name=='env' else 2
        assert f'using {wanted} worker processes' in (root/name/'stderr.log').read_text()
for n in [1,3]:
    result=execute('io_failure_'+str(n),mt,['-nproc',str(n)],limit=512,sentinel=True)
    assert result['returncode']!=0,result
result=execute('worker_failure',mt,['-nproc','3'],limit=64,sentinel=True)
assert result['returncode']!=0,result
assert 'worker failed' in (root/'worker_failure/stderr.log').read_text()
for name,opts in [('missing',['-nproc']),('bad_grid',['-nx','0']),('unknown',['-unknown'])]:
    result=execute(name,mt,opts,sentinel=True)
    assert result['returncode']!=0,result
result=execute('publish_failure',mt,['-nproc','3'],directory_output=True)
assert result['returncode']!=0,result
small_reference=execute('one_point_original',orig,['-nx','1','-ny','1'])
small=execute('cap_at_points',mt,['-nx','1','-ny','1','-nproc','5000'])
assert small['returncode']==0 and small['sha256']==small_reference['sha256']
assert 'using 1 worker processes' in (root/'cap_at_points/stderr.log').read_text()
(root/'PASS.json').write_text(json.dumps({'status':'PASS','cases':len(results),'original_sha256':hashlib.sha256(orig.read_bytes()).hexdigest(),'candidate_sha256':hashlib.sha256(mt.read_bytes()).hexdigest(),'nointerp_oracle':'known integer translation; stock uninitialized-fraction discrepancy retained'}))
print('ALL PASS',len(results),flush=True)
