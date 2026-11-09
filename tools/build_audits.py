from pathlib import Path
import json,sys,random
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from factorlab import audit

def cases():
    def item(name,n,fs,maps,order,e=None):return {'name':name,'input':{'n':n,'factors':[{'scope':s,'values':[str(x) for x in vs]} for s,vs in fs],'evidence':e or {},'map_vars':maps,'order':order}}
    base=[([0,1],[6,5,0,5])]
    r=[item('sum before max: true marginal MAP',2,base,[0],[1,0]),item('max before sum: an upper bound',2,base,[0],[0,1]),item('joint MAP is a different question',2,base,[0,1],[0,1]),item('partition function is a sum',2,base,[],[0,1])]
    star=[([0,i],[2,1,1,2]) for i in range(1,5)]
    r += [item('star: hub first',5,star,[],list(range(5))),item('star: leaves first',5,star,[],[1,2,3,4,0]),item('star: observed hub',5,star,[],[1,2,3,4],{'0':1}),item('star: all leaves are MAP',5,star,[1,2,3,4],[0,1,2,3,4])]
    r += [item('zero evidence mass',2,base,[0],[0],{'1':1}),item('all zero factor',2,[([0,1],[0]*4)],[0],[1,0]),item('unmentioned variable still counts',3,[([0],[1,2])],[],[0,1,2]),item('scalar and disconnected ties',2,[([],["1/3"]),([0],[1,1]),([1],[2,2])],[0],[0,1]),item('all observed: empty elimination',2,base,[],[],{'0':0,'1':0}),item('one variable: exact ties',1,[([0],["2/3","2/3"])],[0],[0])]
    # Separate explicit impossible evidence factor.
    r[8]=item('zero evidence mass',2,[([0,1],[1,2,0,0])],[0],[0],{'1':1})
    rng=random.Random(719)
    for i in range(10):
        n=3+i%3;fs=[([j,j+1],[f'{rng.randrange(5)}/{rng.randrange(1,5)}' for _ in range(4)]) for j in range(n-1)]+[([0,n-1],[1,2,3,4])];maps=list(range(0,n,2));order=list(range(n));rng.shuffle(order);r.append(item('seeded cyclic '+str(i+1),n,fs,maps,order))
    return r

def main():
    root=Path(__file__).resolve().parents[1];records=[{**x,'report':audit(x['input'])} for x in cases()];(root/'docs/audits.json').write_text(json.dumps(records,separators=(',',':'))+'\n');(root/'examples/demo.json').write_text(json.dumps(records[1]['input'],indent=2)+'\n');print('Built',len(records),'reports')
if __name__=='__main__':main()
