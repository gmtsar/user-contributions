#!/usr/bin/env python3
"""Real GMT/GMTSAR synthetic geocoding through the installed xcorr3 launcher."""
import argparse,hashlib,json,os,pathlib,shlex,struct,subprocess
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--prefix',type=pathlib.Path,required=True)
p.add_argument('--fixtures',type=pathlib.Path,required=True)
p.add_argument('--work',type=pathlib.Path,required=True)
a=p.parse_args();root=a.work.resolve();root.mkdir(parents=True,exist_ok=True)
entry=a.prefix.resolve()/'bin/xcorr3';cc=a.prefix.resolve()/'bin/xcorr_cc'
env=dict(os.environ,PATH=str(a.prefix.resolve()/'bin')+':/usr/local/GMTSAR/bin:/usr/local/GMTSAR/gmtsar/csh:'+os.environ['PATH'],PYTHONDONTWRITEBYTECODE='1')
products=['freq_xcorr.dat','azi_offset.grd','rng_offset.grd','azi_offset_ll.grd','rng_offset_ll.grd']
reference=None;records=[]
for mode in ['direct','cli','config','workflow','empty_filter','failed_projection','missing_product']:
    work=root/mode;work.mkdir()
    for name in ['master','aligned']:
        (work/(name+'.SLC')).symlink_to((a.fixtures/(name+'.SLC')).resolve())
        text=(a.fixtures/(name+'.PRM')).read_text()
        (work/(name+'.PRM')).write_text(text+'SC_vel = 7500\nearth_radius = 6371000\nSC_height = 700000\nrng_samp_rate = 64345238\n')
    with (work/'trans.dat').open('wb') as f:
        for y in range(0,1025,16):
            for x in range(0,1537,16):f.write(struct.pack('5d',x,y,0,87+x*.00002,28+y*.00002))
    opts=['master.PRM','aligned.PRM','-nx','7','-ny','5','-xsearch','32','-ysearch','32','-geocode','-snr','100' if mode=='empty_filter' else '0','-psnr','0']
    failure=mode in ['empty_filter','failed_projection','missing_product']
    if failure:
        for name in products:(work/name).write_text('OLD '+name+'\n')
    testenv=env.copy()
    if mode in ['failed_projection','missing_product']:
        tools=work/'tools';tools.mkdir()
        tool=tools/'proj_ra2ll.csh';tool.write_text('#!/bin/sh\nexit '+('7' if mode=='failed_projection' else '0')+'\n');tool.chmod(0o755)
        testenv['PATH']=str(tools)+':'+testenv['PATH']
    if mode=='direct':cmd=[str(cc),*opts]
    elif mode=='cli':cmd=[str(entry),'--backend','cc',*opts]
    else:
        config=work/'config.txt';config.write_text('xcorr_backend = cc\n')
        cmd=[str(entry),'--config',str(config)]
        if mode=='config':cmd+=opts
        else:cmd+=['--run','csh','-f','-c','xcorr '+shlex.join(opts)+'; /bin/true']
    done=subprocess.run(cmd,cwd=work,env=testenv,capture_output=True,timeout=180)
    (work/'stderr.log').write_bytes(done.stderr)
    record=dict(name=mode,returncode=done.returncode,command=cmd)
    if failure:
        record['pass']=done.returncode>0 and all((work/name).read_text()=='OLD '+name+'\n' for name in products)
    else:
        record['pass']=done.returncode==0 and all((work/name).stat().st_size>0 for name in products)
        if record['pass']:
            arrays={}
            for name in products[1:]:
                raw=subprocess.check_output(['gmt','grd2xyz',str(work/name),'-ZTLf'],env=env)
                arrays[name]=np.frombuffer(raw,dtype=np.float32)
                record['pass'] &= bool(np.isfinite(arrays[name]).any())
            if reference is None:reference=(work/'freq_xcorr.dat').read_bytes(),arrays
            record['equal_to_direct_cc']=(work/'freq_xcorr.dat').read_bytes()==reference[0] and all(np.array_equal(arrays[n],reference[1][n],equal_nan=True) for n in arrays)
            record['pass'] &= record['equal_to_direct_cc']
    record['pass']=bool(record['pass']) and not list(work.glob('.xcorr_cc-*'))
    records.append(record)
    (root/'results.json').write_text(json.dumps(dict(passed=all(r['pass'] for r in records),cases=records),indent=2))
    print(('PASS ' if record['pass'] else 'FAIL ')+mode,flush=True)
assert all(r['pass'] for r in records)
