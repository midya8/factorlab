import {test} from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import {audit} from '../docs/engine.js';
const records=JSON.parse(fs.readFileSync(new URL('../docs/audits.json',import.meta.url)));let leaves=0;
function count(x){if(x&&typeof x==='object')Object.values(x).forEach(count);else leaves++}
test('independent BigInt all complete reports',()=>{for(const x of records){assert.deepEqual(audit(x.input),x.report,x.name);count(x.report)}console.log(JSON.stringify({cases:records.length,exactLeaves:leaves}))});
test('strict browser model schema and no mutation',()=>{for(const x of records){let d=structuredClone(x.input),s=JSON.stringify(d);audit(d);assert.equal(JSON.stringify(d),s)}for(const fn of [d=>d.n=true,d=>d.factors[0].values[0]='-1',d=>d.factors[0].values[0]=1,d=>d.order=[0,0],d=>d.evidence={'00':0},d=>d.map_vars=[99],d=>d.extra=1]){let d=structuredClone(records[0].input);fn(d);assert.throws(()=>audit(d))}});
