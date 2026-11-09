import unittest,copy,json,random,itertools,subprocess,sys,tempfile
from fractions import Fraction as F
from pathlib import Path
from math import prod
from factorlab import audit,strict_loads
from tools.build_audits import cases

def direct(d):
    rows={}
    for bits in itertools.product([0,1],repeat=d['n']):
        if any(bits[int(k)]!=v for k,v in d['evidence'].items()):continue
        w=prod((F(f['values'][sum(bits[v]<<j for j,v in enumerate(f['scope']))]) for f in d['factors']),start=F(1));rows[bits]=w
    masses={}
    for bits,w in rows.items():
        k=tuple(bits[v] for v in d['map_vars']);masses[k]=masses.get(k,F(0))+w
    active=[i for i in range(d['n']) if str(i) not in d['evidence']]
    dense={tuple(bits[v] for v in active):w for bits,w in rows.items()};stages=[]
    for v in d['order']:
        j=active.index(v);groups={}
        for bits,w in dense.items():
            k=bits[:j]+bits[j+1:];groups.setdefault(k,[]).append(w)
        dense={k:(max(vs) if v in d['map_vars'] else sum(vs,F(0))) for k,vs in groups.items()};active.remove(v);stages.append((active[:],dense.copy()))
    return rows,masses,next(iter(dense.values())),stages

def random_input(rng,n):
    fs=[]
    for _ in range(rng.randrange(1,7)):
        s=sorted(rng.sample(range(n),rng.randrange(min(n,3)+1)));fs.append({'scope':s,'values':[f'{rng.randrange(6)}/{rng.randrange(1,6)}' for _ in range(2**len(s))]})
    e={str(v):rng.randrange(2) for v in range(n) if rng.random()<.18};a=[v for v in range(n) if str(v) not in e];maps=[v for v in a if rng.random()<.5];order=a[:];rng.shuffle(order)
    return dict(n=n,factors=fs,evidence=e,map_vars=maps,order=order)

class Tests(unittest.TestCase):
    def test_counterexample(self):
        a,b,c,p=[audit(x['input']) for x in cases()[:4]]
        self.assertEqual([a['elimination_value'],b['elimination_value'],c['elimination_value'],p['elimination_value']],['10','11','6','16']);self.assertEqual(a['marginal_map_probability'],'5/8');self.assertEqual(b['upper_bound_gap'],'1');self.assertEqual(a['joint_map_assignments'],[[0,0]]);self.assertEqual([r['bits'] for r in a['map_rows'] if r['optimal']],[[1]])
    def test_star_cost_and_fill(self):
        a,b=[audit(x['input']) for x in cases()[4:6]]
        self.assertEqual(a['partition'],b['partition']);self.assertEqual(a['order_stats']['join_entries'],62);self.assertEqual(b['order_stats']['join_entries'],18);self.assertEqual(a['order_stats']['induced_width'],4);self.assertEqual(b['order_stats']['induced_width'],1);self.assertEqual(len(a['trace'][0]['fill_edges']),6);self.assertEqual(a['minimum_induced_width'],1)
    def test_160_dense_oracles(self):
        rng=random.Random(947)
        for _ in range(160):
            d=random_input(rng,rng.randrange(1,6));r=audit(d);rows,masses,scalar,_=direct(d)
            self.assertEqual(F(r['partition']),sum(rows.values(),F(0)));self.assertEqual(F(r['exact_marginal_map_value']),max(masses.values()));self.assertEqual(F(r['elimination_value']),scalar);self.assertGreaterEqual(scalar,max(masses.values()))
            if r['sum_before_max']:self.assertEqual(scalar,max(masses.values()))
    def test_100_full_partial_factor_oracles(self):
        rng=random.Random(431)
        for _ in range(100):
            d=random_input(rng,rng.randrange(1,6));r=audit(d);_,_,_,stages=direct(d);current=copy.deepcopy(r['sliced_factors'])
            for t,(active,dense) in zip(r['trace'],stages):
                current=[f for f in current if t['variable'] not in f['scope']]+[t['generated']]
                for bits,truth in dense.items():
                    assignment=dict(zip(active,bits));got=prod((F(f['values'][sum(assignment[v]<<j for j,v in enumerate(f['scope']))]) for f in current),start=F(1));self.assertEqual(got,truth)
    def test_all_join_and_reduction_records(self):
        for x in cases():
            r=audit(x['input']);fs={f['id']:f for f in r['sliced_factors']}
            for t in r['trace']:
                for row in t['join_rows']:
                    assignment=dict(zip(t['join_scope'],row['bits']));vals=[]
                    for term in row['terms']:
                        f=fs[term['factor']];idx=sum(assignment[v]<<j for j,v in enumerate(f['scope']));self.assertEqual(idx,term['local_mask']);self.assertEqual(f['values'][idx],term['value']);vals.append(F(term['value']))
                    self.assertEqual(F(row['product']),prod(vals,start=F(1)))
                for row in t['output_rows']:
                    bs=[F(t['join_rows'][i]['product']) for i in row['join_masks']];self.assertEqual(bs,list(map(F,row['branches'])));v=max(bs) if t['operator']=='max' else sum(bs,F(0));self.assertEqual(F(row['value']),v);self.assertEqual(row['maximizers'],[i for i,w in enumerate(bs) if w==v] if t['operator']=='max' else [])
                fs[t['generated']['id']]=t['generated']
    def test_all_order_widths_independent_graph(self):
        for x in cases():
            r=audit(x['input']);a=r['active']
            for p in r['order_census']:
                adj={v:set() for v in a}
                for f in r['sliced_factors']:
                    for u,v in itertools.combinations(f['scope'],2):adj[u].add(v);adj[v].add(u)
                widths=[];entries=0
                for v in p['order']:
                    ns=adj.pop(v);widths.append(len(ns));entries+=2**(len(ns)+1)
                    for u in ns:adj[u].discard(v);adj[u].update(ns-{u})
                self.assertEqual(max(widths,default=0),p['induced_width']);self.assertEqual(entries,p['join_entries'])
    def test_all_orders_on_20_random_graphs(self):
        rng=random.Random(13)
        for _ in range(20):
            d=random_input(rng,4);d['evidence']={};d['order']=list(range(4));d['map_vars']=[]
            for maps in [[],list(range(4)),[0,2]]:
                d['map_vars']=maps;_,masses,_,_=direct(d)
                for order in itertools.permutations(range(4)):
                    d['order']=list(order);r=audit(d);self.assertGreaterEqual(F(r['elimination_value']),max(masses.values()))
                    if r['sum_before_max']:self.assertEqual(F(r['elimination_value']),max(masses.values()))
    def test_40_positive_scale_controls(self):
        rng=random.Random(29)
        for _ in range(40):
            d=random_input(rng,4);a=audit(d);d['factors'].append({'scope':[],'values':['3/2']});b=audit(d)
            for key in ['partition','elimination_value','exact_marginal_map_value','upper_bound_gap']:self.assertEqual(F(b[key]),F(a[key])*F(3,2))
            self.assertEqual(a['marginal_map_probability'],b['marginal_map_probability'])
    def test_disconnected_invalid_can_be_exact(self):
        d=dict(n=2,factors=[dict(scope=[0],values=['1','2']),dict(scope=[1],values=['3','4'])],evidence={},map_vars=[0],order=[0,1]);r=audit(d);self.assertFalse(r['sum_before_max']);self.assertTrue(r['matches_exact_value']);self.assertEqual(r['elimination_value'],'14')
    def test_zero_mass_refuses_probabilities(self):
        r=audit(cases()[8]['input']);self.assertEqual(r['partition'],'0');self.assertFalse(r['probability_defined']);self.assertIsNone(r['marginal_map_probability']);self.assertTrue(all(x['probability'] is None for x in r['joint_rows']))
    def test_unmentioned_variable(self):
        r=audit(cases()[10]['input']);self.assertEqual(r['unit_variables'],[1,2]);self.assertEqual(r['partition'],'12')
    def test_all_evidence(self):
        r=audit(cases()[12]['input']);self.assertEqual(r['trace'],[]);self.assertEqual(r['partition'],'6');self.assertEqual(r['marginal_map_probability'],'1')
    def test_structural_counts_no_zero_sparsity(self):
        d=cases()[4]['input'];a=audit(d)
        for f in d['factors']:f['values']=['0']*len(f['values'])
        b=audit(d);self.assertEqual(a['order_stats'],b['order_stats'])
    def test_no_input_mutation(self):
        for c in cases():d=c['input'];s=copy.deepcopy(d);audit(d);self.assertEqual(d,s)
    def test_strict_keys(self):
        d=cases()[0]['input'];d['extra']=0
        with self.assertRaises(ValueError):audit(d)
    def test_strict_integer_types(self):
        for key,val in [('n',True),('n',2.0),('order',[0,True]),('map_vars',[False]),('evidence',{'1':True})]:
            d=cases()[0]['input'];d[key]=val
            with self.assertRaises(ValueError):audit(d)
    def test_bad_scopes(self):
        for scope in [[1,0],[0,0],[0,9],[True]]:
            d=cases()[0]['input'];d['factors'][0]['scope']=scope
            with self.assertRaises(ValueError):audit(d)
    def test_bad_rationals(self):
        for value in ['-1','1/0','1.5','NaN','1e2','10001','1/10001',1,True,'9'*13]:
            d=cases()[0]['input'];d['factors'][0]['values'][0]=value
            with self.assertRaises(ValueError):audit(d)
    def test_evidence_and_map_scope(self):
        for key,val in [('evidence',{'01':1}),('map_vars',[1,0]),('evidence',{'2':0}),('order',[0])]:
            d=cases()[0]['input'];d[key]=val
            with self.assertRaises(ValueError):audit(d)
    def test_invalid_table_length(self):
        d=cases()[0]['input'];d['factors'][0]['values']=['1']
        with self.assertRaises(ValueError):audit(d)
    def test_strict_json(self):
        for s in ['{"n":1,"n":2}','{"x":NaN}','{"x":Infinity}',' '*32769]:
            with self.assertRaises(ValueError):strict_loads(s)
    def test_cli_invalid_preserves_output(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'bad.json').write_text('{}');(p/'out.json').write_text('keep');r=subprocess.run([sys.executable,'-m','factorlab',str(p/'bad.json'),'--output',str(p/'out.json')],capture_output=True);self.assertEqual(r.returncode,2);self.assertEqual((p/'out.json').read_text(),'keep')
    def test_cli_reproduces(self):
        r=subprocess.run([sys.executable,'-m','factorlab','examples/demo.json'],capture_output=True,text=True);self.assertEqual(r.returncode,0);self.assertEqual(json.loads(r.stdout),audit(strict_loads(Path('examples/demo.json').read_text())))
    def test_fixture_reproduction(self):
        rs=json.loads(Path('docs/audits.json').read_text());self.assertEqual(len(rs),24)
        for r in rs:self.assertEqual(audit(r['input']),r['report'])
if __name__=='__main__':unittest.main()
