#!/usr/bin/env python3
"""Finite certificate algebra, conditional on saved flow enclosures. No ODE replay."""
import json,sys,resource,time
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
if not __debug__:raise SystemExit('Optimized Python is unsupported; run without -O or PYTHONOPTIMIZE.')
A_ROWS=[1,2,3,4,5,6,10,11,12,13,14,15,16]
B_ROWS=[1,2,3,12];C_ROWS=[4,5,6,10,13,14];TWO=sorted(B_ROWS+C_ROWS)
def require(ok,msg):
 if not ok:raise ValueError(msg)
def load(p):
 b=(ROOT/p).read_bytes();require(len(b)<32*1024**2,'input size limit');return json.loads(b)
def rat(v):
 require(isinstance(v,list) and len(v)==2,'rational pair');return F(int(v[0]),int(v[1]))
def dyad(v):
 if isinstance(v,dict):n,e=int(v['n']),int(v['e'])
 else:n,e=int(v[0]),int(v[1])
 require(abs(e)<100000 and len(str(abs(n)))<10000,'dyadic size limit')
 return F(n)*F(2)**e

def interval(v):
 if isinstance(v,dict):
  m,r=dyad(v['mid']),dyad(v['rad']);require(r>=0,'negative radius');a,b=m-r,m+r
 else:a,b=dyad(v[0]),dyad(v[1])
 require(a<=b,'reversed interval');return a,b

def plus(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
 x=[u*v for u in a for v in b];return min(x),max(x)
def ab(a):return max(abs(x) for x in a)
def isum(xs):
 v=(F(0),F(0))
 for x in xs:v=plus(v,x)
 return v

def check(entry,d,expected_center):
 family,row=entry['family'],entry['row'];is_c=family=='C'
 t=d['terminal']['mpfr'] if is_c else d['krawczyk']
 cm=[[interval(v) for v in row_] for row_ in t['C']];jm=[[interval(v) for v in row_] for row_ in t['J']]
 fv=[interval(v[0] if family=='A' else v) for v in t['F0' if is_c else 'F']]
 require(len(cm)==len(jm)==len(fv)==3 and all(len(x)==3 for x in cm+jm),'3x3 dimensions')
 require(all(a==b for rr in cm for a,b in rr),'C must be a point matrix')
 m=[[v[0] for v in rr] for rr in cm]
 det=sum((1 if (j,k)==(1,2) else -1)*m[0][i]*m[1][(i+j)%3]*m[2][(i+k)%3] for i in range(3) for j,k in [(1,2),(2,1)])
 require(det!=0,'singular C')
 z=[[plus((F(i==j),F(i==j)),neg(isum(mul(cm[i][k],jm[k][j]) for k in range(3)))) for j in range(3)] for i in range(3)]
 q=max(sum(ab(v) for v in rr) for rr in z);require(q<1,'contraction q >= 1')
 exp=d['parameter_radius_exponent' if is_c else 'box_radius_power_two'];require(isinstance(exp,int) and -10000<exp<0,'radius exponent')
 require(exp==entry['radius_exponent'],'radius index mismatch');radius=F(2)**exp
 kval=[plus(neg(isum(mul(cm[i][k],fv[k]) for k in range(3))),isum(mul(z[i][j],(-radius,radius)) for j in range(3))) for i in range(3)]
 ratio=max(ab(v) for v in kval)/radius;require(ratio<1,'Krawczyk strict inclusion failed')
 nd=interval(t['nondegeneracy'] if is_c else d['nondegeneracy_interval']);require(not nd[0]<=0<=nd[1],'nondegeneracy contains zero')
 distance=rat(d['all_time_distance2_lower']) if is_c else interval(d['distance_squared_lower'])[0];require(distance>0,'nonpositive distance')
 center=list(map(F,entry['center']));require(center==list(map(F,expected_center)),'catalogue center mismatch')
 if not is_c:require(center==list(map(F,d['center'])),'certificate center mismatch')
 elif 'source_row_text' in d:require(center==list(map(F,d['source_row_text'].split()[:3])),'CAPD source center mismatch')
 elif row==13:require(center==list(map(F,(ROOT/'data/catalogue/C13-row.txt').read_text().split()[:3])),'CAPD13 source center mismatch')
 else:raise ValueError('Missing CAPD center binding')
 require(center[2]>radius and any(abs(x)>radius for x in center[:2]),'parameter domain')
 require(d.get('orbit',d.get('row',row))==row,'row identity')
 if not is_c:
  times=[rat(v) for v in d['steps']];require(times and all(x>0 for x in times),'nonpositive step')
  require(sum(times)==rat(d['integrated_time']),'step sum mismatch')
 else:
  require(d['clock_target']==128 and d['variable_T_included'],'CAPD clock/period domain')
  require(d['validated_c0_slabs']>0 and d['validated_c1_slabs']>0,'CAPD slab counts')
 require(d['status']==('PASS' if is_c else 'CERTIFIED'),'saved status mismatch')
 return {'family':family,'row':row,'det_C_nonzero':True,'q':str(q),'K_radius_ratio':str(ratio),'nondegenerate':True,'distance_squared_lower':str(distance),'radius_exponent':exp,'center':entry['center']}

def run():
 entries=load('data/index.json');require(len(entries)==23,'23 certificates')
 require(len({(x['family'],x['row']) for x in entries})==23,'duplicate certificate')
 for fam,rows in [('A',A_ROWS),('B',B_ROWS),('C',C_ROWS)]:require(sorted(x['row'] for x in entries if x['family']==fam)==rows,'family coverage')
 raw=[s.split()[:3] for s in (ROOT/'data/catalogue/euler-first20.txt').read_text().splitlines()];require(len(raw)==20,'20 input rows')
 results=[check(e,load(e['path']),raw[e['row']-1]) for e in entries]
 a={x['row']:x for x in results if x['family']=='A'};paired=[]
 for y in results:
  if y['family']=='A':continue
  x=a[y['row']]
  require(list(map(F,x['center']))==list(map(F,y['center'])),'same-center comparison')
  require(x['radius_exponent']<=y['radius_exponent'],'nested boxes comparison')
  paired.append(y['row'])
 require(sorted(paired)==TWO,'ten comparisons')
 return {'status':'PASS','certificates':len(results),'same_zero_rows':sorted(paired),'Arb_rows':A_ROWS,'MPFR_rows':B_ROWS,'CAPD_rows':C_ROWS,'results':results,'scope':'Finite exact-rational terminal implications, conditional on saved flow/J enclosures. No ODE, Taylor, source audit or full-trace replay.'}
if __name__=='__main__':
 resource.setrlimit(resource.RLIMIT_CPU,(120,120));t=time.monotonic()
 try:
  result=run();rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024);require(rss<2_000_000_000,'RSS limit');result.update(seconds=time.monotonic()-t,peak_rss_bytes=rss)
  print(json.dumps(result,indent=2))
 except Exception as e:print('FAIL:',str(e),file=sys.stderr);raise SystemExit(1)
