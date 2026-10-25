#!/usr/bin/env python3
import argparse,hashlib,json,sys,time,traceback,datetime
from pathlib import Path
from fractions import Fraction
from flint import arb,ctx
from validated_flow import integrate,finish,encode,encode_mat,norm


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--orbit',type=int,default=1)
    p.add_argument('--bits',type=int,default=320)
    p.add_argument('--order',type=int,default=24)
    p.add_argument('--box-exp',type=int,default=-80)
    p.add_argument('--tolbits',type=int,default=170)
    p.add_argument('--maxstep',default='1/16')
    p.add_argument('--stop-time')
    norms=p.add_mutually_exclusive_group()
    norms.add_argument('--weighted', action='store_true')
    norms.add_argument('--weighted-hessian', action='store_true')
    p.add_argument('--minstep-exp', type=int, default=-28)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    ctx.prec=a.bits;ctx.cap=a.order+4;ctx.threads=1
    data=Path('studies/R011/prestudy/50-files/data/euler_100_4860.txt')
    raw=data.read_bytes();row=raw.decode().splitlines()[a.orbit-1].split()
    vx,vy,T=map(arb,row[:3]);r=arb(2)**a.box_exp
    start=time.monotonic();last=[start];latest=[None]
    mode='hessian' if a.weighted_hessian else a.weighted
    meta=dict(schema='R015.validated-quarter.v3' if a.weighted_hessian else
                    'R015.validated-quarter.v2' if a.weighted else 'R015.validated-quarter.v1',orbit=a.orbit,
              picard_norm='weighted_hessian' if a.weighted_hessian else 'weighted_diagonal' if a.weighted else 'unweighted',
              source_sha256=hashlib.sha256(raw).hexdigest(),center=row[:3],Tstar=row[3],
              precision_bits=a.bits,order=a.order,box_radius_power_two=a.box_exp,
              local_remainder_tolerance_power_two=-a.tolbits,maxstep=a.maxstep,
              minimum_step_power_two=a.minstep_exp,
              started_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              status='RUNNING_NOT_CERTIFIED')
    (a.output/'run.json').write_text(json.dumps(meta,indent=2)+'\n')
    def progress(d):
        latest[0]=d
        if time.monotonic()-last[0]>10:
            d=dict(d,elapsed_seconds=time.monotonic()-start)
            print(json.dumps(d),flush=True)
            (a.output/'progress.json').write_text(json.dumps(d,indent=2)+'\n')
            last[0]=time.monotonic()
    def checkpoint(d):
        tmp=a.output/'checkpoint.tmp';tmp.write_text(json.dumps(d)+'\n');tmp.replace(a.output/'checkpoint.json')
    try:
        f=integrate(vx,vy,T,r,a.order,a.tolbits,Fraction(a.maxstep),progress,
                    Fraction(a.stop_time) if a.stop_time else None,checkpoint,weighted=mode,
                    minstep_exp=a.minstep_exp)
        out=dict(meta,elapsed_seconds=time.monotonic()-start,steps=f['steps'],rejected=f['rejected'],
                 rejection_reasons=f['rejection_reasons'],
                 integrated_time=[str(f['time'].numerator),str(f['time'].denominator)],
                 endpoint={k:encode_mat(f[k]) for k in ['c','Q','rp','rb','R']},
                 distance_squared_lower=encode(f['min_distance_squared_lower']))
        if a.stop_time:
            out['status']='PARTIAL_FLOW_ENCLOSED_NOT_PERIODIC_CERTIFICATE'
        else:
            k=finish(f,T,r,weighted=mode)
            out.update(status='CERTIFIED' if k['passed'] else 'KRAWCZYK_NOT_PROVED',
                       krawczyk={x:encode_mat(k[x]) for x in ['F','J','C','K','Z','endpoint']},
                       krawczyk_ratios=[encode(x) for x in k['ratios']],
                       contraction=encode(k['contraction']),nondegenerate=k['nondegenerate'],
                       nondegeneracy_interval=encode(k['g']),distance_squared_lower=encode(k['distance_squared_lower']))
            print('Krawczyk',k['passed'],'ratios',[float(x) for x in k['ratios']],
                  'contraction',float(k['contraction']),'g',str(k['g']),flush=True)
        (a.output/'result.json').write_text(json.dumps(out,indent=2)+'\n')
        print(out['status'],'seconds',out['elapsed_seconds'],'steps',len(f['steps']),flush=True)
        return 0
    except BaseException as e:
        out=dict(meta,status='INCONCLUSIVE_ERROR',elapsed_seconds=time.monotonic()-start,error=repr(e),last_progress=latest[0])
        (a.output/'result.json').write_text(json.dumps(out,indent=2)+'\n')
        traceback.print_exc();return 1


if __name__=='__main__':
    sys.exit(main())
