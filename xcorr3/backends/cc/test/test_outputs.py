#!/usr/bin/env python3
"""Host fault injection for actual output/postprocess code; no GPU required.

Mock commands test failure propagation, not scientific geocoding accuracy.
Use --real-tools for an additional synthetic mapping with installed GMT/GMTSAR.
"""
import argparse,json,os,pathlib,subprocess,shlex,resource,struct,hashlib,sys
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work',type=pathlib.Path,required=True)
p.add_argument('--real-tools',action='store_true')
a=p.parse_args();root=a.work.resolve();root.mkdir(parents=True,exist_ok=True)
source=pathlib.Path(__file__).resolve().parents[1]
flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','glib-2.0'],text=True))
subprocess.run(['gcc','-std=gnu99','-g','-fsanitize=address,undefined','-c',str(source/'prm_helper.c'),*flags,'-o',str(root/'prm.o')],check=True)
binary=root/'outputs'
subprocess.run(['g++','-std=c++14','-Wall','-Wextra','-Werror','-g','-no-pie','-fsanitize=address,undefined','-I'+str(source),str(source/'test/output_harness.cpp'),str(source/'xcorr_postprocess.cpp'),str(root/'prm.o'),*flags,'-o',str(binary)],check=True)
fake=root/'fake';fake.mkdir(exist_ok=True)
script='''#!/usr/bin/python3
import sys,os,pathlib
name=pathlib.Path(sys.argv[0]).name
stage=sys.argv[1] if name=='gmt' else 'projection'
if stage==os.environ.get('FAIL_STAGE'):sys.exit(7)
if stage=='blockmedian':print(pathlib.Path(sys.argv[2]).read_text(),end='')
elif stage=='grdinfo':
 print('0 1 0 1 nan nan 1 1 2 2' if os.environ.get('NAN_GRID') else '0 1 0 1 1 2 1 1 2 2')
else:
 dest=next((x[2:] for x in sys.argv if x.startswith('-G')),None) if name=='gmt' else sys.argv[3]
 if dest and os.environ.get('MISSING_GRID')!=stage:pathlib.Path(dest).write_text('NEW GRID\\n')
'''
for name in ('gmt','proj_ra2ll.csh'):
    (fake/name).write_text(script);(fake/name).chmod(0o755)
# LD_PRELOAD fault only the specified production operation, once. Rollback can run.
fault=root/'fault.c'
fault.write_text(r'''#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
static int used;
int rename(const char *a,const char *b){
 int (*real)(const char*,const char*)=dlsym(RTLD_NEXT,"rename");
 const char *mode=getenv("RENAME_FAIL");
 if(!used&&mode&&strstr(a,".xcorr_cc-")&&!strstr(a,"backup-")&&strstr(b,mode)){used=1;errno=EIO;return -1;}
 return real(a,b);
}
int fclose(FILE *f){
 int (*real)(FILE*)=dlsym(RTLD_NEXT,"fclose"); int fd=fileno(f);char path[4096],link[64];
 snprintf(link,sizeof(link),"/proc/self/fd/%d",fd);ssize_t n=readlink(link,path,sizeof(path)-1);if(n>=0)path[n]=0;else path[0]=0;
 int result=real(f);
 if(!used&&getenv("CLOSE_FAIL")&&strstr(path,".xcorr_cc-")&&strstr(path,"freq_xcorr.dat")){used=1;errno=EIO;return EOF;}
 return result;
}
int fsync(int fd){
 int (*real)(int)=dlsym(RTLD_NEXT,"fsync");
 if(!used&&getenv("SYNC_FAIL")){used=1;errno=EIO;return -1;}
 return real(fd);
}
''')
subprocess.run(['gcc','-shared','-fPIC',str(fault),'-ldl','-o',str(root/'fault.so')],check=True)
asan=subprocess.check_output(['gcc','-print-file-name=libasan.so'],text=True).strip()
names=['freq_xcorr.dat','azi_offset.grd','rng_offset.grd','azi_offset_ll.grd','rng_offset_ll.grd']
records=[]
def case(name,geo=True,extra=None,expect=None,limit=None,missing=False,badprm=False,real=False,target=None,existing=True):
    w=root/name;w.mkdir(exist_ok=True)
    if existing:
        for item in names:(w/item).write_bytes(('OLD '+item+'\n').encode())
    sentinel={x:((w/x).read_bytes() if existing else None) for x in names}
    if target:
        (w/'freq_xcorr.dat').unlink()
        if target=='directory':(w/'freq_xcorr.dat').mkdir()
        else:(w/'freq_xcorr.dat').symlink_to('/dev/full')
    (w/'master.PRM').write_text('PRF = 2000\nSC_vel = 7500\nearth_radius = 6371000\nSC_height = 700000\nrng_samp_rate = 64345238\n' if not badprm else 'PRF = bad\n')
    if not missing:
        with (w/'trans.dat').open('wb') as f:
            for y in range(0,705,16):
                for x in range(0,705,16):f.write(struct.pack('5d',x,y,0,87+x*.00002,28+y*.00002))
    env=dict(os.environ,PATH=str(fake)+':'+os.environ['PATH'],ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1')
    env.update(extra or {})
    if real:env['PATH']='/usr/local/GMTSAR/bin:/usr/local/GMTSAR/gmtsar/csh:'+os.environ['PATH']
    if extra and any(key in extra for key in ('CLOSE_FAIL','RENAME_FAIL','SYNC_FAIL')):env['LD_PRELOAD']=asan+':'+str(root/'fault.so')
    def restrict():resource.setrlimit(resource.RLIMIT_FSIZE,(limit,limit))
    run=subprocess.run([str(binary),'geo' if geo else 'plain'],cwd=w,env=env,capture_output=True,text=True,timeout=90,preexec_fn=restrict if limit else None)
    sanitizer=any(x in run.stderr for x in ('AddressSanitizer','runtime error:','LeakSanitizer'))
    passed=not sanitizer and not list(w.glob('.xcorr_cc-*'))
    if expect:
        passed=passed and run.returncode>0 and expect in run.stderr
        if target:passed=passed and ((w/'freq_xcorr.dat').is_dir() if target=='directory' else (w/'freq_xcorr.dat').is_symlink())
        else:passed=passed and all(((w/x).read_bytes()==sentinel[x]) if existing else not (w/x).exists() for x in names)
    else:
        passed=passed and run.returncode==0 and (w/'freq_xcorr.dat').read_bytes()!=sentinel['freq_xcorr.dat']
        if geo:passed=passed and all((w/x).read_bytes()!=sentinel[x] for x in names)
        else:passed=passed and all((w/x).read_bytes()==sentinel[x] for x in names[1:])
    records.append(dict(name=name,pass_=bool(passed),returncode=run.returncode,stderr=run.stderr))
case('plain',geo=False)
case('geo_success')
case('raw_success',extra={'RAW_GRID':'1'})
case('missing_trans',missing=True,expect='trans.dat')
case('bad_prm',badprm=True,expect='PRF')
case('empty_filter',extra={'EMPTY_FILTER':'1'},expect='no points passed')
for stage in ('blockmedian','xyz2grd','grdfilter','grdinfo','projection'):
    case('failure_'+stage,extra={'FAIL_STAGE':stage},expect='stage failed')
case('zero_exit_missing_grid',extra={'MISSING_GRID':'projection'},expect='missing output')
case('all_nan_grid',extra={'NAN_GRID':'1'},expect='invalid or all-NaN')
case('write_limit',geo=False,limit=128,expect='flush')
case('immediate_write_failure',geo=False,limit=128,extra={'UNBUFFERED':'1'},expect='write test correlation')
case('close_failure',geo=False,extra={'CLOSE_FAIL':'1'},expect='close')
case('sync_failure',geo=False,extra={'SYNC_FAIL':'1'},expect='sync')
case('publish_failure',extra={'RENAME_FAIL':'rng_offset_ll.grd'},expect='publish')
case('publish_without_old',extra={'RENAME_FAIL':'rng_offset_ll.grd'},expect='publish',existing=False)
case('missing_tool',extra={'PATH':str(root/'absent-tools')},expect='stage failed: gmt')
case('directory_output',geo=False,target='directory',expect='not a regular file')
case('device_symlink',geo=False,target='symlink',expect='not a regular file')
if a.real_tools:case('real_gmt_gmtsar',real=True)
report=dict(pass_=all(x['pass_'] for x in records),cases=records,source_sha256={n:hashlib.sha256((source/n).read_bytes()).hexdigest() for n in ('xcorr_io.h','xcorr_postprocess.cpp','xcorr_postprocess.h','prm_helper.c')})
(root/'results.json').write_text(json.dumps(report,indent=2))
print(f"{sum(x['pass_'] for x in records)}/{len(records)} output cases passed")
for x in records:
    if not x['pass_']:print(x['name'],x['stderr'][-2500:])
sys.exit(0 if report['pass_'] else 1)
