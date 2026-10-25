#!/usr/bin/env python3
"""One sequential command shared by README, CI and clean-clone preflight."""
import os,sys,time,resource,json,subprocess
if not __debug__:raise SystemExit('Optimized Python is unsupported')
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
resource.setrlimit(resource.RLIMIT_CPU,(180,180))
from verify import ROOT,run as verify,require
from regenerate import generate
from negative_controls import run as negative
from check_packet import run as packet
start=time.monotonic()
r=verify();tables=generate()
for name,text in tables.items():require((ROOT/name).read_text()==text,'generated output mismatch '+name)
n=negative();p=packet()
opt=subprocess.run([sys.executable,'-O','-B',str(ROOT/'verify.py')],capture_output=True,text=True,timeout=10)
require(opt.returncode!=0 and 'Optimized Python' in opt.stderr+opt.stdout,'-O must fail explicitly')
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024);require(rss<2_000_000_000,'RSS limit');require(time.monotonic()-start<180,'wall limit')
print(json.dumps({'status':'PASS','certificates':r['certificates'],'same_zero_rows':r['same_zero_rows'],'tables':7,'negative_controls':n['rejected'],'optimized_python_rejected':True,'packet':p,'seconds':time.monotonic()-start,'peak_rss_bytes':rss,'ODE_integrations':0,'full_trace_replays':0,'new_independent_reviews':0},indent=2))
