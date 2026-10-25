#!/usr/bin/env python3
"""Input mutations must be rejected by arithmetic checks, not only hash checks."""
import copy,json
from verify import load,check,require,ROOT

def zero(v):
 if isinstance(v,dict) and 'mid' in v:return {'mid':[0,0],'rad':[0,0]}
 return [{'n':'0','e':0},{'n':'0','e':0}]
def huge(v):
 if isinstance(v,dict) and 'mid' in v:return {'mid':[1,1000],'rad':[0,0]}
 return [{'n':'1','e':1000},{'n':'1','e':1000}]
def run():
 raw=[s.split()[:3] for s in (ROOT/'data/catalogue/euler-first20.txt').read_text().splitlines()];controls=[]
 for family in ('A','B','C'):
  e=next(x for x in load('data/index.json') if x['family']==family);d=load(e['path'])
  check(e,d,raw[e['row']-1])
  for kind in ['singular-C','zero-J','large-F','zero-ND','negative-distance','wrong-center','wrong-radius','invalid-status','reversed-interval']:
   z=copy.deepcopy(d);ee=copy.deepcopy(e);t=z['terminal']['mpfr'] if family=='C' else z['krawczyk']
   if kind=='singular-C':t['C'][0]=[zero(v) for v in t['C'][0]]
   if kind=='zero-J':t['J']=[[zero(v) for v in rr] for rr in t['J']]
   if kind=='large-F':
    key='F0' if family=='C' else 'F'
    if family=='A':t[key][0][0]=huge(t[key][0][0])
    else:t[key][0]=huge(t[key][0])
   if kind=='zero-ND':
    if family=='C':t['nondegeneracy']=zero(t['nondegeneracy'])
    else:z['nondegeneracy_interval']=zero(z['nondegeneracy_interval'])
   if kind=='negative-distance':
    if family=='C':z['all_time_distance2_lower']=['-1','1']
    else:z['distance_squared_lower']={'mid':[-1,0],'rad':[0,0]}
   if kind=='wrong-center':ee['center'][0]='0'
   if kind=='wrong-radius':ee['radius_exponent']-=1
   if kind=='invalid-status':z['status']='NOT_CERTIFIED'
   if kind=='reversed-interval':t['J'][0][0]=[{'n':'2','e':0},{'n':'1','e':0}]
   try:check(ee,z,raw[e['row']-1])
   except (ValueError,KeyError,TypeError) as exc:controls.append({'family':family,'mutation':kind,'rejected':True,'reason':str(exc)})
   else:raise ValueError(f'Accepted mutation {family} {kind}')
 require(len(controls)==27,'control count');return {'status':'PASS','rejected':len(controls),'controls':controls}
if __name__=='__main__':print(json.dumps(run(),indent=2))
