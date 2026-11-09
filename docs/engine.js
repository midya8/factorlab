// Independent BigInt rational finite factor elimination, not a Python bridge.
export class Q {
 constructor(n,d=1n){n=BigInt(n);d=BigInt(d);if(!d)throw Error('zero denominator');if(d<0n){n=-n;d=-d}let a=n<0n?-n:n,b=d;while(b)[a,b]=[b,a%b];this.n=n/a;this.d=d/a}
 add(b){return new Q(this.n*b.d+b.n*this.d,this.d*b.d)}
 sub(b){return new Q(this.n*b.d-b.n*this.d,this.d*b.d)}
 mul(b){return new Q(this.n*b.n,this.d*b.d)}
 div(b){return new Q(this.n*b.d,this.d*b.n)}
 cmp(b){const v=this.n*b.d-b.n*this.d;return v<0n?-1:v>0n?1:0}
 toString(){return this.d===1n?''+this.n:this.n+'/'+this.d}
}
const zero=()=>new Q(0),one=()=>new Q(1),sum=a=>a.reduce((x,y)=>x.add(y),zero()),product=a=>a.reduce((x,y)=>x.mul(y),one()),maximum=a=>a.reduce((x,y)=>x.cmp(y)>=0?x:y),range=n=>Array.from({length:n},(_,i)=>i),sorted=a=>[...a].sort((x,y)=>x-y),unique=a=>[...new Set(a)],bits=(k,n)=>range(n).map(i=>(k>>i)&1),glob=(s,k)=>s.reduce((m,v,j)=>m|(((k>>j)&1)<<v),0),loc=(s,m)=>s.reduce((k,v,j)=>k|(((m>>v)&1)<<j),0),at=(f,m)=>f.values[loc(f.scope,m)];
const keys=(x,ks)=>x&&typeof x==='object'&&!Array.isArray(x)&&Object.keys(x).sort().join('|')===[...ks].sort().join('|');
function parse(d){
 if(!keys(d,['n','factors','evidence','map_vars','order'])||!Number.isInteger(d.n)||d.n<1||d.n>5)throw Error('expected n=1..5 and exact input keys');
 const valid=a=>Array.isArray(a)&&a.every(i=>Number.isInteger(i)&&i>=0&&i<d.n)&&unique(a).length===a.length;
 const ascending=a=>a.every((v,i)=>i===0||v>a[i-1]);
 if(!Array.isArray(d.factors)||d.factors.length<1||d.factors.length>10)throw Error('1..10 factors required');
 const fs=d.factors.map((f,i)=>{if(!keys(f,['scope','values'])||!valid(f.scope)||!ascending(f.scope)||!Array.isArray(f.values)||f.values.length!==2**f.scope.length)throw Error('invalid factor scope or table');return {id:'f'+i,scope:[...f.scope],values:f.values.map(s=>{if(typeof s!=='string'||s.length>12||!/^\d+(?:\/[1-9]\d*)?$/.test(s))throw Error('nonnegative integer/fraction strings required');let a=s.split('/'),v=new Q(a[0],a[1]??1);if(v.n>10000n||v.d>10000n)throw Error('rational limit');return v})}});
 if(!d.evidence||typeof d.evidence!=='object'||Array.isArray(d.evidence)||Object.entries(d.evidence).some(([k,v])=>!range(d.n).map(String).includes(k)||!Number.isInteger(v)||v<0||v>1))throw Error('invalid evidence');
 const active=range(d.n).filter(v=>!(String(v) in d.evidence));
 if(!valid(d.map_vars)||!ascending(d.map_vars)||!d.map_vars.every(v=>active.includes(v)))throw Error('invalid MAP variables');
 if(!valid(d.order)||sorted(d.order).join(',')!==active.join(','))throw Error('order must permute unobserved variables');
 return {fs,active};
}
export function safe(order,maps){let seen=false;for(const v of order){if(maps.includes(v))seen=true;else if(seen)return false}return true}
function edges(scopes){const es=[];for(const s of scopes)for(let i=0;i<s.length;i++)for(let j=i+1;j<s.length;j++)es.push(s[i]+','+s[j]);return unique(es).sort((a,b)=>{let x=a.split(',').map(Number),y=b.split(',').map(Number);return x[0]-y[0]||x[1]-y[1]})}
const edgeLists=es=>es.map(s=>s.split(',').map(Number));
function profile(scopes,order,maps){scopes=scopes.map(s=>[...s]);let es=edges(scopes),steps=[];
 for(const v of order){const bucket=scopes.filter(s=>s.includes(v)),s=sorted(unique([v,...bucket.flat()])),out=s.filter(i=>i!==v),fill=edges([out]).filter(e=>!es.includes(e));es=unique([...es,...fill]).filter(e=>!e.split(',').map(Number).includes(v));es=edges(edgeLists(es));scopes=[...scopes.filter(s=>!s.includes(v)),out];steps.push({variable:v,neighbors:out,fill_edges:edgeLists(fill),edges_after:edgeLists(es),join_entries:2**s.length,output_entries:2**out.length,multiplications:2**s.length*Math.max(bucket.length-1,0),additions:maps.includes(v)?0:2**out.length,comparisons:maps.includes(v)?2**out.length:0})}
 return {order:[...order],sum_before_max:safe(order,maps),induced_width:Math.max(0,...steps.map(s=>s.neighbors.length)),join_entries:steps.reduce((a,s)=>a+s.join_entries,0),peak_join_entries:Math.max(0,...steps.map(s=>s.join_entries)),peak_output_entries:Math.max(0,...steps.map(s=>s.output_entries)),steps};
}
function perms(a){if(!a.length)return [[]];return a.flatMap(v=>perms(a.filter(x=>x!==v)).map(p=>[v,...p]))}
function encoded(x){if(x instanceof Q)return x.toString();if(Array.isArray(x))return x.map(encoded);if(x&&typeof x==='object')return Object.fromEntries(Object.entries(x).map(([k,v])=>[k,encoded(v)]));return x}
export function audit(d){
 const {fs,active}=parse(d),maps=d.map_vars,order=d.order,em=Object.entries(d.evidence).reduce((m,[v,b])=>m|(b<<Number(v)),0);
 const sliced=fs.map(f=>{let scope=f.scope.filter(v=>!(String(v) in d.evidence));return {id:f.id,scope,values:range(2**scope.length).map(k=>at(f,glob(scope,k)|em))}}),absent=active.filter(v=>!sliced.some(f=>f.scope.includes(v)));
 for(const v of absent)sliced.push({id:'u'+v,scope:[v],values:[one(),one()]});
 const profiles=perms(active).map(p=>profile(sliced.map(f=>f.scope),p,maps)),stats=profiles.find(p=>p.order.join(',')===order.join(','));let current=[...sliced],trace=[];
 for(let t=0;t<order.length;t++){const v=order[t],bucket=current.filter(f=>f.scope.includes(v)),scope=sorted(unique([v,...bucket.flatMap(f=>f.scope)])),out=scope.filter(i=>i!==v),joined=range(2**scope.length).map(k=>{const m=glob(scope,k),terms=bucket.map(f=>({factor:f.id,local_mask:loc(f.scope,m),value:at(f,m)}));return {mask:k,bits:bits(k,scope.length),terms,product:product(terms.map(t=>t.value))}}),rows=range(2**out.length).map(k=>{const m=glob(out,k),inds=[0,1].map(b=>loc(scope,m|(b<<v))),vals=inds.map(i=>joined[i].product),z=maps.includes(v)?maximum(vals):sum(vals);return {mask:k,bits:bits(k,out.length),join_masks:inds,branches:vals,value:z,maximizers:maps.includes(v)?[0,1].filter(i=>vals[i].cmp(z)===0):[]}}),generated={id:'g'+t,scope:out,values:rows.map(r=>r.value)};
 current=[...current.filter(f=>!f.scope.includes(v)),generated];trace.push({step:t+1,variable:v,operator:maps.includes(v)?'max':'sum',bucket:bucket.map(f=>f.id),join_scope:scope,output_scope:out,join_rows:joined,output_rows:rows,generated,...stats.steps[t]});}
 const scalar=product(current.map(f=>f.values[0])),marg=range(2**maps.length).map(zero),joint=range(2**active.length).map(k=>{const m=glob(active,k)|em,terms=fs.map(f=>at(f,m)),w=product(terms),idx=loc(maps,m);marg[idx]=marg[idx].add(w);return {mask:m,bits:bits(m,d.n),factor_values:terms,weight:w,map_mask:idx}}),z=sum(joint.map(r=>r.weight)),defined=z.n!==0n,truth=maximum(marg),best=maximum(joint.map(r=>r.weight));for(const r of joint)r.probability=defined?r.weight.div(z):null;
 return encoded({input:d,active,unit_variables:absent,sliced_factors:sliced,initial_edges:edgeLists(edges(sliced.map(f=>f.scope))),trace,final_factors:current,final_multiplications:Math.max(current.length-1,0),elimination_value:scalar,partition:z,probability_defined:defined,exact_marginal_map_value:truth,marginal_map_probability:defined?truth.div(z):null,joint_map_value:best,joint_map_assignments:joint.filter(r=>r.weight.cmp(best)===0).map(r=>r.bits),joint_rows:joint,map_rows:marg.map((w,k)=>({mask:k,bits:bits(k,maps.length),mass:w,probability:defined?w.div(z):null,optimal:w.cmp(truth)===0})),sum_before_max:stats.sum_before_max,matches_exact_value:scalar.cmp(truth)===0,upper_bound_gap:scalar.sub(truth),normalized_elimination_value:defined?scalar.div(z):null,order_stats:stats,order_census:profiles.map(p=>Object.fromEntries(Object.entries(p).filter(([k])=>k!=='steps'))),minimum_induced_width:Math.min(...profiles.map(p=>p.induced_width)),minimum_safe_join_entries:Math.min(...profiles.filter(p=>p.sum_before_max).map(p=>p.join_entries))});
}
