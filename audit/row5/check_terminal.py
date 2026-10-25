"""Row5 terminal binding, using exact functions from own prior R015 audit."""
from pathlib import Path
from fractions import Fraction as Q
import json,resource,time
from dyadic import D,matrix,product
ROOT=Path(__file__).resolve().parent

def pair(q):return [str(q.numerator),str(q.denominator)]

def unpair(x):return Q(int(x[0]),int(x[1]))

def contains_point(interval,q):
    a,b=interval.endpoints();return a<=q<=b

def mv(a,x):return [r[0] for r in product(a,[[v] for v in x])]

def terminal(certificate,c0,c1,d):
    t=certificate['terminal']['mpfr'];f=[D.read(x) for x in t['F0']];j=matrix(t['J']);c=matrix(t['C'])
    assert f==[D.read(c0['final']['endpoint'][i]) for i in (4,5,6)]
    assert j==[d[i] for i in (4,5,6)]
    nd=D.read(t['nondegeneracy']);assert nd==D.read(c1['final']['nondegeneracy'])
    endpoint=[D.read(x) for x in c1['final']['endpoint']];assert nd.contains(endpoint[0]+D(2)*endpoint[2])
    radius=Q(2)**certificate['parameter_radius_exponent'];assert D.read(t['radius']).endpoints()==(radius,radius)
    assert all(x.a==x.b for row in c for x in row)
    cc=[[x.endpoints()[0] for x in row] for row in c]
    det=(cc[0][0]*cc[1][1]*cc[2][2]+cc[0][1]*cc[1][2]*cc[2][0]+cc[0][2]*cc[1][0]*cc[2][1]
         -cc[0][2]*cc[1][1]*cc[2][0]-cc[0][1]*cc[1][0]*cc[2][2]-cc[0][0]*cc[1][2]*cc[2][1])
    cj=product(c,j);w=[[D(i==k)-cj[i][k] for k in range(3)] for i in range(3)]
    cf=mv(c,f);box=D(-1,1,certificate['parameter_radius_exponent'])
    delta=[-v+sum((x*box for x in row),D(0)) for v,row in zip(cf,w)]
    q=max(sum(max(abs(a),abs(b)) for a,b in (x.endpoints() for x in row)) for row in w)
    assert all(D.read(x).contains(y) for a,b in zip(t['W'],w) for x,y in zip(a,b))
    assert all(D.read(x).contains(y) for x,y in zip(t['K_delta'],delta))
    assert contains_point(D.read(t['det_C']),det)
    assert D.read([t['q_upper']]*2).endpoints()[0]>=q
    nlo,nhi=nd.endpoints();ratio=max(max(abs(a),abs(b))/radius for a,b in (x.endpoints() for x in delta))
    passed=det!=0 and q<1 and ratio<1 and (nlo>0 or nhi<0)
    saved=certificate['terminal']['fraction'];assert unpair(saved['q_exact_upper'])==q and unpair(saved['det_C_exact'])==det
    assert [[pair(a),pair(b)] for a,b in (x.endpoints() for x in delta)]==saved['K_delta']
    assert passed==saved['pass']==t['pass']
    return {'terminal_inequalities_pass':passed,'F0_J_ND_bound_to_saved_flow':True,'q_exact_upper':pair(q),
            'q_float':float(q),'strict_inclusion_ratio':float(ratio),'det_C_nonzero':det!=0,
            'ND_excludes_zero':nlo>0 or nhi<0,'radius_exponent':certificate['parameter_radius_exponent']}

if __name__=='__main__':
    start=time.monotonic()
    full=json.loads((ROOT/'full-check.json').read_text());assert full['status']=='PASS' and full['limit']==0
    cases={r['phase']:r for r in full['cases']}
    d=[[D.pairs(v) for v in r] for r in cases['c1']['final_D']]
    metadata=json.loads((ROOT/'cache/full-row5-c1-r140.json').read_text())
    assert d==[[D.pairs(v) for v in r] for r in metadata['verification']['last_parameter_derivative']]
    certificate=json.loads((ROOT/'cache/certificate-row5-r140.json').read_text())
    result=terminal(certificate,cases['c0'],cases['c1'],d)
    assert result['terminal_inequalities_pass']
    result.update(status='PASS',old_final_D_equals_independent_reduction=True,elapsed_seconds=time.monotonic()-start,max_rss_bytes_macos=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    destination=ROOT/'terminal-check.json';assert not destination.exists();destination.write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='q_exact_upper'})
