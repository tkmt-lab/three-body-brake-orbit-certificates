"""Read-only audit of saved row4 evidence; no ODE integration or producer import.

Uses integer shifts for every enclosure comparison, interval product and outward
grid operation. Only time choices and the finite decimal center use Fraction.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import argparse,gzip,hashlib,json,resource,time
from dyadic import D,matrix,product
ROOT=Path(__file__).resolve().parent
RAW=ROOT/'external-evidence'

def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def fraction(x):
    n=int(x['n']);e=int(x['e'])
    return Q(0) if not n else (Q(n<<e) if e>=0 else Q(n,1<<(-e)))
def same(a,b):
    if a==b:return True
    if isinstance(a,dict) and isinstance(b,dict):
        if a.keys()!=b.keys():return False
        if set(a)=={'n','e'}:return fraction(a)==fraction(b)
        return all(same(a[k],b[k]) for k in a)
    if isinstance(a,list) and isinstance(b,list):return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return False
def canonical_hash(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def records(path,hashes,name):
    digest=hashlib.sha256();size=0;lines=0
    with gzip.open(path,'rb') as f:
        for raw in f:
            digest.update(raw);size+=len(raw);lines+=1
            assert raw.endswith(b'\n') and raw.strip(),(name,'partial or blank line')
            yield json.loads(raw)
    hashes[name]={'sha256':digest.hexdigest(),'bytes':size,'lines':lines}
def math_records(stream):
    for row in stream:
        assert row['kind'] in ('smoke','setup','initial','step','progress','final')
        if row['kind']!='progress':yield row
def mathematical(row):
    excluded=('elapsed_seconds','cpu_seconds','max_rss_kib_linux') if row['kind']=='final' else ()
    return {k:v for k,v in row.items() if k not in excluded}

def check_call(call,phase):
    assert call['branch']=='normal','Only observed normal branch is certified here'
    assert call['factor_before_shrink'] is None
    assert same(call['trial_H'],call['returned_H'])
    H=D.read(call['trial_H']);assert H.a>0
    assert len(call['input_box'])==9
    fields=['polynomial_range','predicted_margin','enclosure','evaluated_remainder_at_trial_H']
    for field in fields:
        obj=call[field];assert len(obj['state'])==9
        assert obj['local_c1'] is None if phase=='c0' else (len(obj['local_c1'])==9 and all(len(r)==9 for r in obj['local_c1']))
    counts=Counter();minimum=None
    for kind in (['state'] if phase=='c0' else ['state','local_c1']):
        arrays=[call[key][kind] for key in fields]
        indexed=[[(i,None,vals[i]) for i in range(9)] if kind=='state' else [(i,j,vals[i][j]) for i in range(9) for j in range(9)] for vals in arrays]
        for entries in zip(*indexed):
            i,j=entries[0][:2];pp,margin,E,rem=[D.read(a[2]) for a in entries]
            exponent=min(v.e for v in [pp,margin,E,rem])
            pl,pu=pp.aligned(exponent);ml,mu=margin.aligned(exponent);el,eu=E.aligned(exponent);rl,ru=rem.aligned(exponent)
            assert ml<0<mu and el<=pl+ml and pu+mu<=eu,(call['step'],kind,i,j,'P+margin not in E')
            assert ml<rl<=ru<mu,(call['step'],kind,i,j,'strict remainder containment')
            initial=D.read(call['input_box'][i]) if kind=='state' else D(i==j)
            assert pp.contains(initial),(call['step'],kind,i,j,'polynomial range excludes time zero')
            assert E.contains(initial),(call['step'],kind,i,j,'enclosure excludes time zero')
            slack=D(min(rl-ml,mu-ru),e=exponent)
            if minimum is None or (slack-minimum).b<0:minimum=slack
            counts[kind]+=1
    return dict(counts),{'n':str(minimum.a),'e':minimum.e}

def check_distances(raw_tube,raw_saved,previous=None):
    tube=[D.read(v) for v in raw_tube];assert len(tube)==9
    distances=[dx.square()+dy.square() for dx,dy in [(tube[2]-tube[0],tube[3]-tube[1]),(-D(2)*tube[0]-tube[2],-D(2)*tube[1]-tube[3]),(-tube[0]-D(2)*tube[2],-tube[1]-D(2)*tube[3])]]
    lower=min(v.endpoints()[0] for v in distances);saved=fraction(raw_saved)
    assert lower>0,'pairwise distance not positive'
    assert 0<saved<=lower,'invalid saved distance lower bound'
    assert previous is None or saved<=previous,'running distance lower bound increased'
    return saved

def final_duration(prefix,H,native):
    selected=Q(128)-prefix;lo,hi=H.endpoints()
    assert selected>0,'nonpositive selected duration'
    assert 0<lo<=selected<=hi,'selected duration outside H'
    assert native.endpoints()[0]<=128<=native.endpoints()[1],'target outside native clock'
    return {'prefix':[str(prefix.numerator),str(prefix.denominator)],'selected_last_duration':[str(selected.numerator),str(selected.denominator)],'selected_equals_left':selected==lo,'selected_duration_positive':True,'final_H':H.pairs_out(),'target':128}

def check_anchor(anchor,call,row,previous_endpoint):
    assert anchor['final_enclosure_call']==call['call'],'stale final call anchor'
    assert same(anchor['accepted_H'],row['h']) and same(row['h'],call['returned_H'])
    assert same(call['input_box'],anchor['input_box']) and same(anchor['input_box'],row['input_box']) and same(row['input_box'],previous_endpoint)
    assert same(call['enclosure']['state'],anchor['tube']) and same(anchor['tube'],row['all_time_tube'])

def run_case(cfg,limit,write):
    phase=cfg['phase'];tag='full-'+phase;start=time.monotonic();rawhash={}
    oldpath=Path(cfg['legacy_local_path']);runmeta=json.loads((ROOT/'cache'/tag/'audit-run.json').read_text())
    cm=json.loads((ROOT/'cache/cloud-manifest.json').read_text());bycloud={v['path']:v for v in cm['files']}
    files={'stdout':RAW/tag/'audit-stdout.jsonl.gz','stderr':RAW/tag/'audit-stderr.jsonl.gz','old':oldpath}
    compressed={}
    for name,path in files.items():
        expected=bycloud['inputs/'+phase+'.jsonl.gz' if name=='old' else tag+'/audit-'+name+'.jsonl.gz']
        digest=sha(path);assert digest==expected['sha256'] and path.stat().st_size==expected['bytes']
        compressed[name]={'path':str(path),'sha256':digest,'bytes':path.stat().st_size}
    old=iter(math_records(records(oldpath,rawhash,'old')));audit=iter(records(files['stderr'],rawhash,'stderr'))
    headers={};previous_endpoint=None;previous_time=D(0);clock=Q(0);before_clock=Q(0);count=0;last_call=0;gaps=[];nonsingular=[];components=Counter();minimums={};calls_total=0;branches=Counter();final=None;running_min=None
    derivative=[[D(int((i in (4,6) and j==0) or (i in (5,7) and j==1) or (i==8 and j==2))) for j in range(3)] for i in range(9)]
    for row in math_records(records(files['stdout'],rawhash,'stdout')):
        assert final is None,'record after final';kind=row['kind'];legacy=next(old)
        assert same(mathematical(row),mathematical(legacy)),('old record differs',phase,kind,row.get('index'))
        if kind in ('smoke','setup','initial'):
            assert count==0 and kind not in headers
            headers[kind]=row
            if kind=='setup':
                assert set(headers)=={'smoke','setup'}
                assert all(row[k]==cfg[k] for k in ['phase','precision','order','tol_exp'])
                assert D.read(row['parameter_radius']).endpoints()==(Q(2)**cfg['radius'],)*2
            if kind=='initial':
                assert set(headers)=={'smoke','setup','initial'}
                a,b,t=map(Q,cfg['row_text'].split()[:3]);center=[Q(-1),Q(0),Q(1),Q(0),a,b,a,b,t]
                assert len(row['box'])==9
                for i,bound in enumerate(row['box']):
                    lo,hi=D.read(bound).endpoints();rad=Q(2)**cfg['radius'] if phase=='c1' and i in (4,5,6,7,8) else Q(0)
                    assert lo<=center[i]-rad<=center[i]+rad<=hi
                assert phase!='c0' or same(row['box'],headers['setup']['center'])
                assert t-(Q(2)**cfg['radius'] if phase=='c1' else 0)>0 and (a!=0 or b!=0)
                previous_endpoint=row['box']
            continue
        if kind=='step':
            count+=1;assert row['index']==count and row['phase']==phase and previous_endpoint is not None
            returns=[]
            while True:
                event=next(audit);assert event['step']==count and event['logger_environment_unchanged'] is True
                if event['audit_kind']=='accepted_step':anchor=event;break
                assert event['audit_kind']=='enclosure_return' and event['call']>last_call
                if event['call']!=last_call+1:gaps.append([last_call,event['call']])
                last_call=event['call'];returns.append(event)
            assert returns;call=returns[-1];calls_total+=len(returns);branches[call['branch']]+=1
            assert call['order']==cfg['order'] and call['mpfr_precision']==cfg['precision']
            check_anchor(anchor,call,row,previous_endpoint)
            checked,minimum=check_call(call,phase);components.update(checked)
            H=D.read(row['h']);native=D.read(row['time']);before=D.read(call['time'])
            assert before==previous_time and native.contains(before+H)
            before_clock=clock;clock+=H.endpoints()[0];assert native.endpoints()[0]<=clock<=native.endpoints()[1]
            if count<cfg['required_steps']:assert H.a==H.b,'nonterminal common duration must be singleton'
            if H.a!=H.b:nonsingular.append(count)
            running_min=check_distances(row['all_time_tube'],row['distance2_running_lower'],running_min)
            if phase=='c1':
                L=matrix(row['local_jacobian']);G=matrix(row['global_parameter_derivative']);assert len(L)==9 and all(len(v)==9 for v in L) and len(G)==9 and all(len(v)==3 for v in G)
                propagated=product(L,derivative,512);derivative=[[x.meet(y) for x,y in zip(a,b)] for a,b in zip(propagated,G)]
            write({'phase':phase,'step':count,'final_call':call['call'],'successful_return_count':len(returns),'branch':call['branch'],'components':checked,'minimum_slack':minimum,
                   'old_math_sha256':canonical_hash(mathematical(legacy)),'state_tube_and_input_match':True,'uniform_H_obligation_pass':True,
                   'margin_zero_components':sum(checked.values()),'time0_state_components':9,'time0_local_C1_components':81 if phase=='c1' else 0,'positive_pairwise_distances':3})
            previous_endpoint=row['endpoint'];previous_time=native;last_H=H
            if limit and count==limit:break
            if count%512==0:print(phase,count,round(time.monotonic()-start,3),flush=True)
        elif kind=='final':
            final=row;assert row['steps']==count==cfg['required_steps'] and row['complete'] and row['positive']
            assert same(row['endpoint'],previous_endpoint) and D.read(row['time'])==previous_time
            assert fraction(row['distance2_lower'])==running_min
        else:raise AssertionError(kind)
    connection=None
    if not limit:
        assert final is not None and next(old,None) is None and next(audit,None) is None
        connection=final_duration(before_clock,last_H,previous_time)
        for name in ['stdout','stderr']:
            assert rawhash[name]['sha256']==runmeta['outputs'][name]['raw_sha256'] and rawhash[name]['bytes']==runmeta['outputs'][name]['raw_bytes']
        oldmeta=json.loads((ROOT/'cache'/('old-'+phase+'.json')).read_text())
        # Locate the original native runtime receipt without assuming a top-level layout.
        def find(value):
            if isinstance(value,dict):
                if 'stdout_sha256' in value and 'binary_sha256' in value:yield value
                for v in value.values():yield from find(v)
            elif isinstance(value,list):
                for v in value:yield from find(v)
        matching=[v for v in find(oldmeta) if v['stdout_sha256']==rawhash['old']['sha256']]
        assert matching and all(v['binary_sha256']==cfg['legacy_binary_sha256'] for v in matching)
    return {'phase':phase,'status':'PREFIX_PASS' if limit else 'PASS','steps':count,'components':dict(components),'successful_returns':calls_total,'final_branch_counts':dict(branches),
            'call_id_gaps':gaps,'nondegenerate_H_steps':nonsingular,'connection':connection,'compressed_files':compressed,'raw_streams':rawhash,
            'final_D':[[v.pairs_out() for v in r] for r in derivative] if phase=='c1' else None,'final':final,'elapsed_seconds':time.monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=0);p.add_argument('--output',required=True);args=p.parse_args()
    assert __debug__, 'Run without -O';start=time.monotonic();assert __import__('shutil').disk_usage(ROOT).free >= 30000000000;resource.setrlimit(resource.RLIMIT_CPU,(1500,1500));dest=ROOT/args.output;assert not dest.exists()
    cfg=json.loads((ROOT/'cache/config.json').read_text());results=[]
    with dest.with_suffix('.steps.jsonl').open('w') as out:
        for case in cfg['cases']:results.append(run_case(case,args.limit,lambda row:out.write(json.dumps(row,separators=(',',':'))+'\n')))
    final={'status':'PASS','kind':'saved-evidence inspection only; no ODE recomputation','limit':args.limit,'cases':results,'elapsed_seconds':time.monotonic()-start,'max_rss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dest.write_text(json.dumps(final,indent=2)+'\n');print({k:v for k,v in final.items() if k!='cases'})
