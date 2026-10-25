#!/usr/bin/env python3
"""Explicit opt-in bounded reintegration. Never imported by default checks."""
from pathlib import Path
import argparse,json,subprocess,sys,os,time,signal,shlex,importlib.util
ROOT=Path(__file__).resolve().parents[1]
if not __debug__:raise SystemExit('Run without -O')
p=argparse.ArgumentParser();p.add_argument('--run-heavy',action='store_true');p.add_argument('--family',choices=['A','B','C'],required=True);p.add_argument('--row',type=int,required=True);p.add_argument('--seconds',type=int,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--phase',choices=['c0','c1'],default='c1');p.add_argument('--binary',type=Path);a=p.parse_args()
if not a.run_heavy:raise SystemExit('No integration run. Supply --run-heavy explicitly after reading REPRODUCING.md.')
if not 1<=a.seconds<=21600:raise SystemExit('seconds must be in [1,21600]')
entry=next((x for x in json.loads((ROOT/'data/index.json').read_text()) if x['family']==a.family and x['row']==a.row),None)
if entry is None:raise SystemExit('No adopted certificate for family/row')
out=a.output.resolve()
if out.exists():raise SystemExit('Output exists; use a new directory')
# Keep arbitrary writes away from distributed input directories.
if ROOT in out.parents and ROOT/'work' not in out.parents:raise SystemExit('Inside repository, output must be below work/')
d=json.loads((ROOT/entry['path']).read_text())
if a.family=='A':
 if importlib.util.find_spec('flint') is None:raise SystemExit('Optional prerequisite missing: python-flint 0.9.0. No automatic installation.')
 runner=ROOT/'optional/arb/run_certificate.py'
 # The frozen driver has one input-path constant. Adapt only that constant in a new work copy.
 text=runner.read_text();old="data=Path('studies/R011/prestudy/50-files/data/euler_100_4860.txt')"
 if text.count(old)!=1:raise SystemExit('Unknown frozen input adapter')
 out.mkdir(parents=True);work=out/'driver.py';work.write_text(text.replace(old,'data=Path('+repr(str(ROOT/'data/catalogue/euler-first20.txt'))+')'))
 cmd=[sys.executable,'-B',str(work),'--orbit',str(a.row),'--bits',str(d['precision_bits']),'--order',str(d['order']),'--box-exp',str(d['box_radius_power_two']),'--tolbits',str(-d['local_remainder_tolerance_power_two']),'--maxstep',d['maxstep'],'--minstep-exp',str(d['minimum_step_power_two']),'--output',str(out/'result')]
 if d.get('picard_norm')=='weighted_hessian':cmd.append('--weighted-hessian')
 elif d.get('picard_norm')=='weighted_diagonal':cmd.append('--weighted')
else:
 if not a.binary or not a.binary.resolve().is_file():raise SystemExit('Supply a locally compiled --binary; no automatic build or fetch.')
 out.mkdir(parents=True)
 if a.family=='B':cmd=[str(a.binary.resolve()),str(ROOT/'data/catalogue/euler-first20.txt'),str(a.row),str(d['precision_bits']),str(d['order']),str(d['box_radius_power_two']),str(d['local_remainder_tolerance_power_two']),'1000000',str(a.seconds),str(out/'result')]
 else:
  cfg=json.loads((ROOT/f'data/settings/C{a.row:02}-{a.phase}.json').read_text());rowfile=out/'row.txt';rowfile.write_text((ROOT/'data/catalogue/euler-first20.txt').read_text().splitlines()[a.row-1]+'\n')
  cmd=[str(a.binary.resolve()),str(rowfile),a.phase,str(cfg['precision']),str(cfg['order']),str(cfg['tol_exp']),'1000000',str(a.seconds)]
  if a.row!=13:cmd.append(str(d['parameter_radius_exponent']))
start=time.monotonic();peak=0;reason=None
with (out/'stdout.txt').open('w') as f,(out/'stderr.txt').open('w') as e:
 child=subprocess.Popen(cmd,cwd=out,stdout=f,stderr=e,start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',PYTHONPATH=str(ROOT/'optional/arb')))
 while child.poll() is None:
  info=subprocess.run(['ps','-axo','pgid=,rss='],capture_output=True,text=True)
  rss=sum(int(s.split()[1])*1024 for s in info.stdout.splitlines() if len(s.split())==2 and int(s.split()[0])==child.pid);peak=max(peak,rss)
  if time.monotonic()-start>a.seconds or rss>2_000_000_000:reason='time/RSS limit';os.killpg(child.pid,signal.SIGTERM);break
  time.sleep(.25)
 try:rc=child.wait(timeout=3)
 except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);rc=child.wait()
(out/'execution.json').write_text(json.dumps({'command':cmd,'exit_code':rc,'stop':reason,'wall_seconds':time.monotonic()-start,'peak_process_group_rss_bytes':peak,'memory_cap_bytes':2_000_000_000,'input_excerpt_note':'20-row excerpt; whole-file source hash differs from historical catalogue hash','not_a_new_independent_implementation':True},indent=2)+'\n')
raise SystemExit(1 if rc or reason else 0)
