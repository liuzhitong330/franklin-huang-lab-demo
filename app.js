(function () {
  /* GPL-3.0. New exploratory implementation; source study credited in README. */
  "use strict";
  const POLICIES = {
    RP_Tumor: { label: "Tumor RP specimens", tissues: ["Paired_T", "Unpaired_T"], note: "Paired and unpaired tumor-source radical prostatectomy (RP) specimens; no biopsies or grossly normal samples." },
    Biopsy: { label: "Biopsy specimens", tissues: ["Biopsy"], note: "Biopsies only; this collection group is too small for the four-patient shortlist rule." },
    Grossly_normal: { label: "Grossly normal RP specimens", tissues: ["Paired_N"], note: "Grossly normal tissue from paired prostatectomies, not healthy controls; the source paper notes possible tumor contamination." }
  };
  const average = xs => xs.length ? xs.reduce((a,b) => a+b, 0) / xs.length : null;
  const finite = x => typeof x === "number" && Number.isFinite(x);
  function summarize(data, settings, gene) {
    const index = data.genes.indexOf(gene);
    if (index < 0 || !POLICIES[settings.policy]) throw new Error("Unknown gene or specimen policy");
    const [target, reference] = settings.contrast.split("|");
    const selected = data.groups.filter(g => POLICIES[settings.policy].tissues.includes(g.tissue));
    const donors = [...new Set(selected.map(g => g.donor))].sort();
    function group(donor, state) {
      const gs = selected.filter(g => g.donor === donor && g.state === state);
      const n = gs.reduce((s,g) => s+g.n, 0);
      const valid = n > 0 && gs.every(g => finite(g.mean[index]) && finite(g.detect[index]));
      return {
        n, mean: valid ? gs.reduce((s,g) => s+g.n*g.mean[index],0)/n : null,
        detect: valid ? gs.reduce((s,g) => s+g.n*g.detect[index],0)/n : null,
        sources: [...new Set(gs.flatMap(g => g.sourceSamples))], valid
      };
    }
    const included = [], omitted = [];
    for (const donor of donors) {
      const a = group(donor, target), b = group(donor, reference);
      const row = { donor, a, b };
      if (a.n < settings.minimum || b.n < settings.minimum || !a.valid || !b.valid) {
        row.reason = !a.n || !b.n ? "missing state" : !a.valid || !b.valid ? "unavailable gene" : "below cell threshold";
        omitted.push(row);
      } else {
        row.delta = a.mean - b.mean;
        row.detectionDelta = a.detect - b.detect;
        included.push(row);
      }
    }
    const n = included.length, balanced = average(included.map(r => r.delta));
    const nA = included.reduce((s,r)=>s+r.a.n,0), nB = included.reduce((s,r)=>s+r.b.n,0);
    const pooled = n ? included.reduce((s,r)=>s+r.a.n*r.a.mean,0)/nA - included.reduce((s,r)=>s+r.b.n*r.b.mean,0)/nB : null;
    const detection = average(included.map(r => r.detectionDelta));
    const loo = n > 1 ? included.map((_,i) => average(included.filter((_,j) => i !== j).map(r => r.delta))) : [];
    const positive = included.filter(r => r.delta > 0).length;
    const looMin = loo.length ? Math.min(...loo) : null, looMax = loo.length ? Math.max(...loo) : null;
    const candidate = n >= 4 && balanced >= settings.effect && positive/n >= 0.75 && looMin > 0;
    const reason = n < 4 ? "Low coverage" : candidate ? "Candidate" : balanced <= 0 ? "Opposite / zero" : "Review";
    return { gene, target, reference, included, omitted, n, nA, nB, balanced, pooled, detection, loo, looMin, looMax, positive, candidate, reason };
  }
  function ranked(data, settings) {
    return data.genes.map(g => summarize(data,settings,g)).sort((a,b) =>
      Number(b.candidate)-Number(a.candidate) || (b.balanced ?? -Infinity)-(a.balanced ?? -Infinity) || a.gene.localeCompare(b.gene));
  }
  function csv(data, settings) {
    const fields = ["gene","specimen_policy","contrast","minimum_cells_per_state","minimum_balanced_difference","eligible_patients","target_cells","reference_cells","balanced_mean_log_difference","pooled_mean_log_difference","balanced_detection_difference","positive_patients","loo_min","loo_max","review","included_donor_evidence","omitted_donor_evidence","source_object_sha256","units"];
    const rows = ranked(data,settings).map(r => [
      r.gene,settings.policy,settings.contrast,settings.minimum,settings.effect,r.n,r.nA,r.nB,
      r.balanced,r.pooled,r.detection,r.positive,r.looMin,r.looMax,r.reason,
      JSON.stringify(r.included.map(d => ({donor:d.donor,n_target:d.a.n,n_reference:d.b.n,delta:d.delta,detection_delta:d.detectionDelta,sources:[...d.a.sources,...d.b.sources]}))),
      JSON.stringify(r.omitted.map(d => ({donor:d.donor,n_target:d.a.n,n_reference:d.b.n,reason:d.reason,sources:[...d.a.sources,...d.b.sources]}))),
      data.meta.sourceSha256,"difference of mean log1p(CP10K), not log fold change"
    ]);
    const quote = x => '"' + String(x ?? "").replaceAll('"','""') + '"';
    return [fields,...rows].map(row => row.map(quote).join(",")).join("\r\n") + "\r\n";
  }
  if (typeof module !== "undefined" && module.exports) { module.exports = {summarize,ranked,csv,POLICIES}; return; }
  const data = window.HUANG_DATA;
  if (!data || !data.groups || !data.meta) { document.getElementById("cohort-note").textContent = "The cohort summary could not load. Please reload the page."; return; }
  const $ = id => document.getElementById(id);
  const format = (v,places=2) => finite(v) ? (v > 0 ? "+" : "") + v.toFixed(places) : "—";
  const integer = n => n.toLocaleString("en-US");
  const escape = v => String(v).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  let gene = "PIGR";
  const settings = () => ({policy:$("tissue").value,contrast:$("contrast").value,minimum:Number($("minimum").value),effect:Number($("effect").value)});
  for (const [value,p] of Object.entries(POLICIES)) {
    const o = document.createElement("option"); o.value=value; o.textContent=p.label; $("tissue").append(o);
  }
  for (const g of data.genes) { const o=document.createElement("option");o.value=g;$("target-options").append(o); }
  for (const g of ["PIGR","MMP7","KLK3","KRT5","IL8"]) {
    const b=document.createElement("button");b.type="button";b.textContent=g;b.addEventListener("click",()=>select(g));$("suggestions").append(b);
  }
  $("cohort-note").textContent = integer(data.meta.cells)+" deposited epithelial cells, "+data.meta.donors+" patients and "+data.meta.sourceSamples+" source sample groups. Biopsy replicates are already collapsed in the deposited identifiers; these groups are not a count of independent physical specimens. The fixed 20-marker panel compares only Club, luminal (LE) and basal (BE) states. Author-labelled tumor states are retained in the source inventory but excluded from these contrasts.";
  function select(g) {
    const found=data.genes.find(x=>x.toUpperCase()===g.trim().toUpperCase());
    if (!found) { $("input-error").textContent="This fixed panel does not include that gene. Choose one of the 20 listed markers; no transcriptome-wide search is implied."; return; }
    gene=found;$("target-select").value=found;$("input-error").textContent="";render();
  }
  $("target-form").addEventListener("submit",e=>{e.preventDefault();select($("target-select").value);});
  for (const id of ["tissue","contrast","minimum","effect"]) $(id).addEventListener("change",()=>{ $("export-status").textContent="";render(); });
  function draw(r) {
    const svg=$("hero-viz"), w=Math.min(600, Math.max(300,svg.parentElement.clientWidth-36)), left=w<440?82:108, right=58;
    const rows = r.included.map(d=>({label:d.donor,value:d.delta,count:d.a.n+"/"+d.b.n}));
    if(r.n) rows.push({label:"Cell-pooled",value:r.pooled,special:"pooled"},{label:"Equal-patient",value:r.balanced,special:"balanced"});
    const h=Math.max(190,rows.length*31+105);
    svg.setAttribute("viewBox","0 0 "+w+" "+h);
    const vals=[0,...rows.map(x=>x.value),...r.loo];
    let low=Math.min(...vals), high=Math.max(...vals), pad=Math.max(.15,(high-low)*.15);
    low-=pad;high+=pad;const x=v=>left+(v-low)/(high-low)*(w-left-right);
    const bottom=55+rows.length*31;
    let s='<title>'+escape(gene+": "+r.target+" minus "+r.reference+", "+r.n+" eligible patients")+'</title>';
    s+='<text x="0" y="18" font-size="11" fill="#666">Patient</text><text x="'+(w-2)+'" y="18" text-anchor="end" font-size="10" fill="#666">n '+escape(r.target)+' / '+escape(r.reference)+'</text>';
    s+='<line x1="'+x(0)+'" x2="'+x(0)+'" y1="28" y2="'+(bottom-20)+'" stroke="#b8b8b8" stroke-dasharray="3 3"/>';
    rows.forEach((row,i)=>{
      const y=42+i*31;
      s+='<text x="0" y="'+(y+4)+'" font-size="11" fill="#444">'+escape(row.label)+'</text>';
      if(row.special==="balanced"&&r.loo.length)s+='<line x1="'+x(r.looMin)+'" x2="'+x(r.looMax)+'" y1="'+y+'" y2="'+y+'" stroke="#333" stroke-width="2"/>';
      if(row.special==="pooled")s+='<rect x="'+(x(row.value)-4)+'" y="'+(y-4)+'" width="8" height="8" fill="#777"/>';
      else s+='<circle cx="'+x(row.value)+'" cy="'+y+'" r="'+(row.special?5:3.5)+'" fill="'+(row.special?"#1f7a8c":"#333")+'"><title>'+escape(row.label+": "+format(row.value,4))+'</title></circle>';
      s+='<text x="'+(w-2)+'" y="'+(y+4)+'" text-anchor="end" font-size="10" fill="#666">'+escape(row.count||format(row.value))+'</text>';
    });
    if(!r.n)s+='<text x="'+(w/2)+'" y="83" text-anchor="middle" font-size="12" fill="#555">No patients meet both coverage requirements.</text>';
    for(let i=0;i<=4;i++){const v=low+(high-low)*i/4;s+='<text x="'+x(v)+'" y="'+bottom+'" text-anchor="middle" font-size="10" fill="#666">'+v.toFixed(1)+'</text>';}
    s+='<text x="'+((left+w-right)/2)+'" y="'+(bottom+24)+'" text-anchor="middle" font-size="11" fill="#555">Mean log-expression difference</text>';
    svg.innerHTML=s;
  }
  function render() {
    const st=settings(), r=summarize(data,st,gene), results=ranked(data,st), candidates=results.filter(x=>x.candidate);
    draw(r);
    $("metric-one").textContent=r.n;$("metric-two").textContent=integer(r.nA+r.nB);$("metric-three").textContent=format(r.balanced);
    const coverage=r.n<4?"Fewer than four eligible patients: descriptive only; no candidate designation.":r.candidate?"Passes the current exploratory follow-up rules.":"Does not pass all current positive-contrast follow-up rules.";
    $("readout").innerHTML='<p><strong>'+escape(gene)+'</strong> · '+escape(POLICIES[st.policy].label)+'. '+escape(coverage)+'</p>'+
      (r.n?'<p>Equal-patient Δ '+format(r.balanced)+'; cell-pooled Δ '+format(r.pooled)+'. Positive in '+r.positive+'/'+r.n+' patients. Mean detection difference '+format(r.detection*100,1)+' percentage points. '+(r.loo.length?'Leave-one-out range '+format(r.looMin)+' to '+format(r.looMax)+'.':'Leave-one-out is unavailable with one patient.')+'</p>':'')+
      '<p><strong>What this lets the lab do:</strong> '+(r.candidate?'Take '+escape(gene)+' forward as a measured expression contrast to check in independent specimens.':'Inspect patient coverage and discordance before using '+escape(gene)+' as a shared-state readout.')+'</p>'+
      '<p><strong>Where Cathy adds value:</strong> pairing patient-level interpretation with transparent coverage and sensitivity checks, building on my PDAC single-cell analysis.</p>';
    $("coverage").innerHTML='<p>'+escape(POLICIES[st.policy].note)+'</p><ul>'+[...r.included.map(d=>'<li>'+escape(d.donor)+': included; '+escape(r.target)+' '+d.a.n+', '+escape(r.reference)+' '+d.b.n+'; Δ '+format(d.delta,3)+'. Source: '+escape([...new Set([...d.a.sources,...d.b.sources])].join(", "))+'.</li>'),...r.omitted.map(d=>'<li>'+escape(d.donor)+': omitted ('+escape(d.reason)+'); '+escape(r.target)+' '+d.a.n+', '+escape(r.reference)+' '+d.b.n+'.</li>')].join("")+'</ul>';
    $("review-summary").textContent=candidates.length+" of "+data.genes.length+" author-selected markers pass the current follow-up rules. "+(r.n<4?"Coverage is insufficient for a shortlist, regardless of effect size.":candidates.length?"Candidates: "+candidates.map(x=>x.gene).join(", ")+". These are measurement priorities, not novel or validated biomarkers.":"No marker passes all criteria; relax a preference only for sensitivity analysis, not to manufacture evidence.");
    $("review-body").replaceChildren();
    for(const row of results) {
      const tr=document.createElement("tr");if(row.gene===gene)tr.className="selected";
      const td=document.createElement("td"),b=document.createElement("button");b.type="button";b.className="gene-link";b.textContent=row.gene;b.setAttribute("aria-label","Inspect "+row.gene);b.addEventListener("click",()=>{select(row.gene);$("target-select").focus();$("target-form").scrollIntoView({behavior:window.matchMedia("(prefers-reduced-motion: reduce)").matches?"auto":"smooth",block:"start"});});td.append(b);tr.append(td);
      for(const v of [format(row.balanced,4),row.n?row.positive+"/"+row.n:"—",format(row.looMin,4),row.reason]){const c=document.createElement("td");c.textContent=v;if(v==="Candidate")c.className="candidate";tr.append(c);}
      $("review-body").append(tr);
    }
  }
  $("download").addEventListener("click",()=>{
    const st=settings(),url=URL.createObjectURL(new Blob([csv(data,st)],{type:"text/csv;charset=utf-8"}));
    const a=document.createElement("a");a.href=url;a.download="huang-marker-review-"+st.policy+"-"+st.contrast.replace("|","-")+"-n"+st.minimum+".csv";document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),10000);
    $("export-status").textContent="Exported the current 20-gene review, with unrounded estimates, included and omitted patient evidence, source hash and active thresholds.";
  });
  window.addEventListener("resize",()=>draw(summarize(data,settings(),gene)));
  render();
})();
