#!/usr/bin/env python3
"""Deterministic table fragments from saved records, with explicit manual text."""
import argparse,json,re,math
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
from verify import ROOT,load,F,rat,dyad,interval,require,A_ROWS,B_ROWS,C_ROWS,TWO

def dec(v):return Decimal(v.numerator)/Decimal(v.denominator)
def fixed(v,n):
 with localcontext() as c:c.prec=100;return format(dec(v).quantize(Decimal(1).scaleb(-n),rounding=ROUND_HALF_EVEN),f'.{n}f')
def sci(v,n):
 with localcontext() as c:
  c.prec=100;d=dec(v);e=d.adjusted();s=format(d.scaleb(-e).quantize(Decimal(1).scaleb(-n),rounding=ROUND_HALF_EVEN),f'.{n}f');return '$'+s+r'\times10^{'+str(e)+'}$'
def square_root(v,n,down=False):
 scale=10**n;num=(v.numerator*scale*scale)//v.denominator;k=math.isqrt(num)
 if not down and 4*v*scale*scale>(2*k+1)**2:k+=1
 return str(k//scale)+'.'+str(k%scale).zfill(n)
def radius(e):return '$2^{'+str(e)+'}$'
def line(c):return ' & '.join(map(str,c))+r'\\'
def generate():
 index=load('data/index.json');ds={f:{e['row']:load(e['path']) for e in index if e['family']==f} for f in 'ABC'};A,B,C=(ds[f] for f in 'ABC');generated={}
 templates={int(p.name.split('-')[0]):p for p in (ROOT/'data/table-templates').glob('*.tex')}
 # Table1 is an explicitly manual account of proof layers; counts checked against certificate index.
 generated[1]=templates[1].read_text()
 runs=[]
 for label,records in [('A, rows 2, 3',[A[r] for r in (2,3)]),('A, other successes',[A[r] for r in A_ROWS if r not in (2,3)]),('B, all successes',list(B.values()))]:
  settings={(d['precision_bits'],d['order'],d['box_radius_power_two'],d['local_remainder_tolerance_power_two']) for d in records};require(len(settings)==1,'uniform family settings');p,o,r,t=settings.pop();runs.append(line([label,p,o,radius(r),radius(t)]))
 for rows in ([4,13,14],[5,6,10]):
  for phase in (0,1):
   settings={(d['precision'],d['order'],d['tol_exp']) for d in [load(f'data/settings/C{r:02}-c{phase}.json') for r in rows]};require(len(settings)==1,'C settings');p,o,t=settings.pop()
   rs=sorted({C[r]['parameter_radius_exponent'] for r in rows},reverse=True)
   text='Point enclosure' if phase==0 else ' or '.join(map(radius,rs))
   runs.append(line(['C, rows '+', '.join(map(str,rows))+r': $C^'+str(phase)+'$',p,o,text,radius(t)]))
 generated[2]='\n'.join(runs)
 bm={x['orbit']:x for x in load('data/status/B.json')};asum=load('data/status/A.json');rows=[]
 for r in range(1,21):
  a='Cert' if r in A else 'Uncert';require(r in A or r in asum['uncertified_rows'],'A status basis')
  b='Cert' if r in B else 'Partial' if bm[r]['attempted'] else 'Notrun'
  if r in C:c='Cert'
  elif r in (11,16):require(load(f'data/status/C{r}.json')['proof_status'].startswith('PARTIAL'),'hold basis');c='Hold'
  elif r in (7,15,17,20):
   require(not load(f'data/status/C{r}.json')['complete'],'incomplete basis');c=r'Uncert/\Partial' if r==15 else 'Partial'
   if r==15:require(not load('data/status/C15-terminal-failure.json')['certified'],'failed terminal basis')
  elif r in (8,9,18,19):require(load(f'data/status/C{r}.json')['pilot'],'pilot basis');c='Notrun^{*}'
  else:c='Notrun'
  rows.append(line([r,*['$\\'+x+'$' for x in (a,b,c)],'Yes' if r in TWO else '---']))
 rows.extend([r'\midrule',line(['Certified rows',len(A),len(B),len(C),len(TWO)])]);generated[3]='\n'.join(rows)
 generated[4]='\n'.join(line([r,f"{len(A[r]['steps']):,}",radius(A[r]['box_radius_power_two']),square_root(interval(A[r]['distance_squared_lower'])[0],6,True)]) for r in A_ROWS)
 generated[5]='\n'.join(line([r,f"{len(B[r]['steps']):,}",radius(B[r]['box_radius_power_two']),square_root(interval(B[r]['distance_squared_lower'])[0],9),sci(interval(B[r]['contraction_upper'])[1],2)]) for r in B_ROWS)
 generated[6]='\n'.join(line([r,f"{C[r]['validated_c0_slabs']:,}",f"{C[r]['validated_c1_slabs']:,}",radius(C[r]['parameter_radius_exponent']),sci(rat(C[r]['all_time_distance2_lower']),9),fixed(dyad(C[r]['terminal']['mpfr']['q_upper']),10)]) for r in C_ROWS)
 x1=load('data/catalogue/crosswalk-first.json');x2=load('data/catalogue/crosswalk-added.json');matched={r['row']:r for r in load('data/catalogue/matched-records.json')['rows']};xs={r['row']:(r,False) for r in x1['rows']};xs.update({r['row']:(r,True) for r in x2['rows']})
 free=[s.split() for s in (ROOT/'data/catalogue/free-first20.txt').read_text().splitlines()];euler=[s.split() for s in (ROOT/'data/catalogue/euler-first20.txt').read_text().splitlines()];lines=[];positive=[];negative=[]
 for r in TWO:
  d,added=xs[r];require(list(map(F,d['free']))==list(map(F,free[r-1])),'free tokens');require(list(map(F,d['euler']))==list(map(F,euler[r-1])),'Euler tokens');cells=[r]
  for i,key in enumerate(['H2024-full','H2024-selfdual'] if added else ['H2024-free24582','H2024-selfdual236']):
   ms=d['matches'][key]
   if not ms:cells.append('Not matched' if i==0 else 'Earlier label unconfirmed');continue
   require(len(ms)==1,'unique saved match');nr=ms[0]['old_row'];tokens=d['matched_records'][key][0]['tokens'] if added else next(v['original_tokens'] for v in matched[r]['matching_records'] if v['catalogue']==key and v['old_row']==nr)
   require(len(tokens)==len(free[r-1]),'match dimensions');require(all(abs(F(a)-F(b))<F('5e-79') for a,b in zip(tokens,free[r-1],strict=True)),'positive crosswalk decimal agreement');cells.append(nr if i==0 else r'self-dual \#'+str(nr))
  (negative if cells[1]=='Not matched' else positive).append(r);lines.append(line(cells))
 require(len(positive)==7 and negative==[4,10,12],'crosswalk scope');generated[7]='\n'.join(lines)
 output={}
 for i,p in templates.items():output['generated/'+p.name]=generated[i] if i==1 else p.read_text().replace('%% GENERATED_ROWS',generated[i])
 tex=(ROOT/'paper/article.tex').read_text();blocks=[m.group()+'\n' for m in re.finditer(r'\\begin\{(table|longtable)\}.*?\\end\{\1\}',tex,re.S)];require(len(blocks)==7,'seven paper tables')
 for i,p in templates.items():require(output['generated/'+p.name]==blocks[i-1],f'paper table {i} differs from regenerated fragment')
 claims=load('claims.json');require(claims['Arb_rows']==A_ROWS and claims['MPFR_rows']==B_ROWS and claims['CAPD_rows']==C_ROWS and claims['same_zero_rows']==TWO,'claim rows')
 return output
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();output=generate()
 for name,text in output.items():
  path=ROOT/name
  if args.check:require(path.read_text()==text,'generated file mismatch '+name)
  else:path.parent.mkdir(exist_ok=True);path.write_text(text)
 print(json.dumps({'status':'PASS','tables':[1,2,3,4,5,6,7],'table1':'manual qualitative proof layers, with certificate row counts bound by claims','tables2to6':'numeric cells regenerated; retained incomplete statuses are saved-record evidence','table7':'7 positive token comparisons recomputed; 3 nonmatches are saved records, not an exhaustive new catalogue search','write':not args.check}))
