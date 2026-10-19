/* GPL-3.0. Independent regression and edge-case checks for the browser analysis. */
"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const root = path.resolve(__dirname,"..");
const app = require(path.join(root,"app.js"));
const data = JSON.parse(fs.readFileSync(path.join(root,"patient-marker-data.json"),"utf8"));
const settings = {policy:"RP_Tumor",contrast:"Club|LE",minimum:20,effect:0.25};
let checks = 0;
function close(actual,expected,label) {
  assert.ok(Number.isFinite(actual) && Math.abs(actual-expected)<=1e-11*Math.max(1,Math.abs(expected)),label);
  checks++;
}
const golden = {
  PIGR: {balanced:2.3843852610210163,pooled:2.5956945068625012,detection:0.6216915214323755,looMin:2.2630429280114046,looMax:2.580562440040905,candidate:true},
  MMP7: {balanced:2.6413943794153125,pooled:2.5770627265898844,detection:0.7292449543481621,looMin:2.364718880839915,looMax:2.7936312085458765,candidate:true},
  KRT5: {balanced:0.2496671441893921,pooled:0.2148416722072662,detection:0.1297187762954201,looMin:0.22573319134193123,looMax:0.3060440902270351,candidate:false},
  KLK3: {balanced:-4.086267150090712,pooled:-4.43211516797175,detection:-0.6444083520326717,looMin:-4.176941201732327,looMax:-3.8799999530597673,candidate:false}
};
for (const [gene,expected] of Object.entries(golden)) {
  const actual = app.summarize(data,settings,gene);
  assert.equal(actual.n,4);
  assert.equal(actual.nA,268);assert.equal(actual.nB,1607);
  assert.deepEqual(actual.included.map(x=>x.donor),["PR5196","PR5249","PR5251","PR5269"]);
  assert.equal(actual.candidate,expected.candidate);
  for(const key of ["balanced","pooled","detection","looMin","looMax"]) close(actual[key],expected[key],gene+" "+key);
}
assert.deepEqual(app.ranked(data,settings).filter(x=>x.candidate).map(x=>x.gene),["MMP7","PIGR","CP","CD74","LTF","KRT15","SERPINB1","IL8"]);
function fixture(deltas, targetSizes, referenceSizes) {
  return {meta:{sourceSha256:"synthetic-unit-test"},genes:["TEST"],groups:deltas.flatMap((delta,i)=>[
    {donor:"D"+i,tissue:"Paired_T",state:"Club",sourceSamples:["S"+i],n:targetSizes?.[i]??20,mean:[1+delta],detect:[0.5]},
    {donor:"D"+i,tissue:"Paired_T",state:"LE",sourceSamples:["S"+i],n:referenceSizes?.[i]??20,mean:[1],detect:[0.2]}
  ])};
}
const none=app.summarize(fixture([]),settings,"TEST");
assert.equal(none.n,0);assert.equal(none.balanced,null);assert.equal(none.pooled,null);assert.equal(none.looMin,null);assert.equal(none.candidate,false);
const single=app.summarize(fixture([1]),settings,"TEST");
assert.equal(single.n,1);assert.equal(single.loo.length,0);assert.equal(single.looMin,null);assert.equal(single.candidate,false);
assert.equal(app.summarize(fixture([1,1,1]),settings,"TEST").candidate,false);
assert.equal(app.summarize(fixture([0.249999,0.249999,0.249999,0.249999]),settings,"TEST").candidate,false,"Do not threshold rounded display values");
assert.equal(app.summarize(fixture([0.25,0.25,0.25,0.25]),settings,"TEST").candidate,true);
assert.equal(app.summarize(fixture([1,1,1,-0.1]),settings,"TEST").candidate,true,"Exactly75% positive");
assert.equal(app.summarize(fixture([4,0.01,-0.01,-0.01]),settings,"TEST").candidate,false,"Insufficient positive donors");
assert.equal(app.summarize(fixture([10,0.1,0.1,-0.3]),settings,"TEST").candidate,false,"LOO sign reversal");
const unavailable=fixture([1,1,1,1]);unavailable.groups[0].mean[0]=null;
const missing=app.summarize(unavailable,settings,"TEST");
assert.equal(missing.n,3);assert.equal(missing.omitted[0].reason,"unavailable gene");assert.equal(missing.candidate,false);
const weighted=app.summarize(fixture([1,3],[20,60],[60,20]),settings,"TEST");
close(weighted.balanced,2,"Equal-donor weighting");close(weighted.pooled,2.5,"Pooling the same included cells");
assert.throws(()=>app.summarize(data,settings,"NOT_A_GENE"));
assert.throws(()=>app.summarize(data,{...settings,policy:"UNKNOWN"},"PIGR"));
function parseCSV(text) {
  const rows=[];let row=[],value="",quoted=false;
  for(let i=0;i<text.length;i++) {
    const c=text[i];
    if(c==='"') {if(quoted&&text[i+1]==='"'){value+='"';i++;}else quoted=!quoted;}
    else if(c===','&&!quoted){row.push(value);value="";}
    else if((c==='\r'||c==='\n')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(value);rows.push(row);row=[];value="";}
    else value+=c;
  }
  assert.equal(quoted,false,"Unterminated CSV quote");
  if(value||row.length){row.push(value);rows.push(row);}
  return rows;
}
const rows=parseCSV(app.csv(data,settings));
assert.equal(rows.length,21);
for(const row of rows.slice(1))assert.equal(row.length,rows[0].length);
const records=rows.slice(1).map(row=>Object.fromEntries(rows[0].map((key,i)=>[key,row[i]])));
const exported=records.find(r=>r.gene==='PIGR');
close(Number(exported.balanced_mean_log_difference),golden.PIGR.balanced,"Unrounded CSV estimate");
assert.equal(exported.minimum_cells_per_state,"20");assert.equal(exported.minimum_balanced_difference,"0.25");
assert.equal(exported.source_object_sha256,data.meta.sourceSha256);
assert.equal(JSON.parse(exported.included_donor_evidence).length,4);
assert.equal(JSON.parse(exported.omitted_donor_evidence).length,4);
const sandbox={window:{}};
vm.runInNewContext(fs.readFileSync(path.join(root,"data.js"),"utf8"),sandbox,{timeout:1000});
assert.ok(JSON.stringify(sandbox.window.HUANG_DATA)===JSON.stringify(data),"Browser bundle must equal JSON (including provenance metadata)");
console.log(JSON.stringify({status:"PASS",baselineGenes:4,numericalGoldenChecks:checks,edgeCases:12,csvRows:20,bundleEqualsJSON:true},null,2));
