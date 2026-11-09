"""Exact bounded binary factor elimination. No third-party dependencies."""
from fractions import Fraction as F
from itertools import permutations, combinations
from math import prod
import json,re

def integer(x): return type(x) is int

def rational(x):
    if not isinstance(x,str) or len(x)>12 or not re.fullmatch(r'[0-9]+(?:/[1-9][0-9]*)?',x): raise ValueError('values must be nonnegative integer/fraction strings, <=12 characters')
    q=F(x)
    if q.numerator>10000 or q.denominator>10000: raise ValueError('reduced numerator/denominator must be <=10000')
    return q

def indices(a,n):
    return isinstance(a,list) and all(integer(i) and 0<=i<n for i in a) and len(a)==len(set(a))

def validate(d):
    if not isinstance(d,dict) or set(d)!={'n','factors','evidence','map_vars','order'}: raise ValueError('expected n, factors, evidence, map_vars, order only')
    n=d['n']
    if not integer(n) or not 1<=n<=5:raise ValueError('n must be integer 1..5')
    fs=d['factors']
    if not isinstance(fs,list) or not 1<=len(fs)<=10:raise ValueError('1..10 factors required')
    factors=[]
    for i,f in enumerate(fs):
        if not isinstance(f,dict) or set(f)!={'scope','values'} or not indices(f['scope'],n) or f['scope']!=sorted(f['scope']):raise ValueError('factor scope must be sorted distinct indices')
        if not isinstance(f['values'],list) or len(f['values'])!=2**len(f['scope']):raise ValueError('factor values length must be 2^scope length')
        factors.append({'id':f'f{i}','scope':f['scope'][:],'values':[rational(x) for x in f['values']]})
    e=d['evidence']
    if not isinstance(e,dict) or any(k not in [str(i) for i in range(n)] or not integer(v) or v not in (0,1) for k,v in e.items()):raise ValueError('evidence uses canonical index keys and integer bits')
    active=[i for i in range(n) if str(i) not in e]
    if not indices(d['map_vars'],n) or d['map_vars']!=sorted(d['map_vars']) or not set(d['map_vars'])<=set(active):raise ValueError('map_vars must be sorted unobserved indices')
    if not indices(d['order'],n) or sorted(d['order'])!=active:raise ValueError('order must contain each unobserved index once')
    return n,factors,{int(k):v for k,v in e.items()},active

def global_mask(scope,mask):return sum(((mask>>j)&1)<<v for j,v in enumerate(scope))
def local_mask(scope,mask):return sum(((mask>>v)&1)<<j for j,v in enumerate(scope))
def value(f,mask):return f['values'][local_mask(f['scope'],mask)]
def strict_loads(text):
    if not isinstance(text,str) or len(text.encode('utf-8'))>32768:raise ValueError('input exceeds 32KiB')
    def pairs(xs):
        r={}
        for k,v in xs:
            if k in r:raise ValueError('duplicate JSON key')
            r[k]=v
        return r
    def bad(x):raise ValueError('nonfinite JSON')
    return json.loads(text,object_pairs_hook=pairs,parse_constant=bad)

def safe_order(order,maps):
    seen=False
    for v in order:
        if v in maps:seen=True
        elif seen:return False
    return True

def graph(scopes):return {tuple(x) for s in scopes for x in combinations(s,2)}

def structure(scopes,order,maps):
    scopes=[s[:] for s in scopes];edges=graph(scopes);steps=[]
    for v in order:
        bucket=[s for s in scopes if v in s];scope=sorted({v}|{i for s in bucket for i in s});out=[i for i in scope if i!=v]
        fill=sorted(set(combinations(out,2))-edges);edges=(edges|set(fill));edges={e for e in edges if v not in e}
        scopes=[s for s in scopes if v not in s]+[out]
        steps.append({'variable':v,'neighbors':out,'fill_edges':[list(x) for x in fill],'edges_after':[list(x) for x in sorted(edges)],'join_entries':2**len(scope),'output_entries':2**len(out),'multiplications':2**len(scope)*max(len(bucket)-1,0),'additions':2**len(out) if v not in maps else 0,'comparisons':2**len(out) if v in maps else 0})
    return {'order':order[:],'sum_before_max':safe_order(order,maps),'induced_width':max((len(s['neighbors']) for s in steps),default=0),'join_entries':sum(s['join_entries'] for s in steps),'peak_join_entries':max((s['join_entries'] for s in steps),default=0),'peak_output_entries':max((s['output_entries'] for s in steps),default=0),'steps':steps}

def audit(d):
    n,fs,e,active=validate(d);maps=d['map_vars'];order=d['order'];em=sum(bit<<v for v,bit in e.items())
    sliced=[]
    for f in fs:
        s=[v for v in f['scope'] if v not in e]
        sliced.append({'id':f['id'],'scope':s,'values':[value(f,global_mask(s,k)|em) for k in range(2**len(s))]})
    # Explicit unit factors preserve counting measure for absent variables.
    absent=[v for v in active if not any(v in f['scope'] for f in sliced)]
    sliced += [{'id':f'u{v}','scope':[v],'values':[F(1),F(1)]} for v in absent]
    profiles=[structure([f['scope'] for f in sliced],list(p),maps) for p in permutations(active)]
    stats=next(p for p in profiles if p['order']==order);current=[{'id':f['id'],'scope':f['scope'][:],'values':f['values'][:]} for f in sliced];trace=[]
    for t,v in enumerate(order):
        bucket=[f for f in current if v in f['scope']];scope=sorted({v}|{i for f in bucket for i in f['scope']});out=[i for i in scope if i!=v];joined=[];jvalues=[]
        for k in range(2**len(scope)):
            m=global_mask(scope,k);terms=[{'factor':f['id'],'local_mask':local_mask(f['scope'],m),'value':value(f,m)} for f in bucket];z=prod((x['value'] for x in terms),start=F(1));jvalues.append(z);joined.append({'mask':k,'bits':[(k>>j)&1 for j in range(len(scope))],'terms':terms,'product':z})
        rows=[]
        for k in range(2**len(out)):
            m=global_mask(out,k);inds=[local_mask(scope,m|(b<<v)) for b in (0,1)];vals=[jvalues[i] for i in inds];z=max(vals) if v in maps else sum(vals,F(0));rows.append({'mask':k,'bits':[(k>>j)&1 for j in range(len(out))],'join_masks':inds,'branches':vals,'value':z,'maximizers':[i for i,x in enumerate(vals) if x==z] if v in maps else []})
        generated={'id':f'g{t}','scope':out,'values':[r['value'] for r in rows]};current=[f for f in current if v not in f['scope']]+[generated]
        trace.append({'step':t+1,'variable':v,'operator':'max' if v in maps else 'sum','bucket':[f['id'] for f in bucket],'join_scope':scope,'output_scope':out,'join_rows':joined,'output_rows':rows,'generated':generated,**stats['steps'][t]})
    scalar=prod((f['values'][0] for f in current),start=F(1));joint=[];marg=[F(0)]*(2**len(maps))
    for k in range(2**len(active)):
        m=global_mask(active,k)|em;terms=[value(f,m) for f in fs];w=prod(terms,start=F(1));idx=local_mask(maps,m);marg[idx]+=w;joint.append({'mask':m,'bits':[(m>>v)&1 for v in range(n)],'factor_values':terms,'weight':w,'map_mask':idx})
    z=sum((r['weight'] for r in joint),F(0));truth=max(marg);bestjoint=max(r['weight'] for r in joint)
    for r in joint:r['probability']=r['weight']/z if z else None
    margrows=[{'mask':k,'bits':[(k>>j)&1 for j in range(len(maps))],'mass':w,'probability':w/z if z else None,'optimal':w==truth} for k,w in enumerate(marg)]
    # Invalid interleavings yield an upper bound, not a feasible MAP assignment.
    report={'input':d,'active':active,'unit_variables':absent,'sliced_factors':sliced,'initial_edges':[list(x) for x in sorted(graph([f['scope'] for f in sliced]))],'trace':trace,'final_factors':current,'final_multiplications':max(len(current)-1,0),'elimination_value':scalar,'partition':z,'probability_defined':bool(z),'exact_marginal_map_value':truth,'marginal_map_probability':truth/z if z else None,'joint_map_value':bestjoint,'joint_map_assignments':[r['bits'] for r in joint if r['weight']==bestjoint],'joint_rows':joint,'map_rows':margrows,'sum_before_max':stats['sum_before_max'],'matches_exact_value':scalar==truth,'upper_bound_gap':scalar-truth,'normalized_elimination_value':scalar/z if z else None,'order_stats':stats,'order_census':[{k:v for k,v in p.items() if k!='steps'} for p in profiles],'minimum_induced_width':min(p['induced_width'] for p in profiles),'minimum_safe_join_entries':min(p['join_entries'] for p in profiles if p['sum_before_max'])}
    def encode(x):
        if isinstance(x,F):return str(x)
        if isinstance(x,list):return [encode(a) for a in x]
        if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
        return x
    return encode(report)
