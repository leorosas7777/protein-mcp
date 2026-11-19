"""The HTML template of build_console_page, kept apart so console.py stays readable.

A Python string rather than a data file on purpose: the server is run from the source tree with
`uv run protein-mcp`, and a template that lives in a .py file cannot go missing from a package
install. `/*__DATA__*/null` is replaced with the embedded dossier.
"""

TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Target console</title>
<style>
:root{
  --bg:#04070a; --panel:#070d10; --grid:#0e2b22; --line:#15382c;
  --p:#35e39a; --p-dim:#1b7a55; --a:#ffab3d; --a-dim:#8a5c20; --r:#ff5f5a;
  --ink:#bff3da; --muted:#5f8f7c; --faint:#38604f;
  --mono:ui-monospace,"Cascadia Mono",Consolas,"DejaVu Sans Mono",monospace;
}
*{box-sizing:border-box}
html,body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--mono);font-size:12px;
  line-height:1.45;letter-spacing:.02em}
body::after{content:"";position:fixed;inset:0;pointer-events:none;z-index:99;
  background:repeating-linear-gradient(to bottom,rgba(53,227,154,.045) 0 1px,transparent 1px 3px)}
.u{text-transform:uppercase;letter-spacing:.1em}
.g{color:var(--p)} .am{color:var(--a)} .rd{color:var(--r)} .mu{color:var(--muted)} .fa{color:var(--faint)}
a{color:var(--p);text-decoration:none;border-bottom:1px solid #16402f}
a:hover{border-bottom-color:var(--p)}
.pill{display:inline-block;border:1px solid var(--a);color:var(--a);border-radius:11px;
  padding:1px 10px;font-size:10px;text-transform:uppercase;letter-spacing:.13em;white-space:nowrap}
.pill.p{border-color:var(--p);color:var(--p)}
.pill.d{border-color:var(--faint);color:var(--muted)}
.frame{border:1px solid var(--line);background:var(--panel);position:relative}
.frame>header{display:flex;gap:10px;align-items:center;padding:6px 10px;border-bottom:1px solid var(--line);
  flex-wrap:wrap}
.frame>header .t{font-size:10px;letter-spacing:.14em;text-transform:uppercase;color:var(--p)}
.frame>header .r{margin-left:auto;font-size:10px;color:var(--faint);letter-spacing:.1em}
.body{padding:9px 10px}
.band{border-bottom:1px solid var(--line);background:linear-gradient(180deg,#071015,#04070a);
  padding:8px clamp(10px,2vw,20px);position:sticky;top:0;z-index:40}
.band .row1{display:flex;gap:12px;align-items:center;flex-wrap:wrap}
.logo{font-size:15px;letter-spacing:.18em;color:var(--p);text-transform:uppercase}
.sub{color:var(--muted);font-size:10.5px;letter-spacing:.09em}
.readouts{display:flex;flex-wrap:wrap;margin-top:8px}
.ro{padding:2px 14px 2px 12px;border-left:1px solid var(--grid)}
.ro:first-child{border-left:0;padding-left:0}
.ro b{display:block;font-size:15px;line-height:1.15;font-weight:400}
.ro span{font-size:9.5px;letter-spacing:.11em;text-transform:uppercase;color:var(--faint)}
.wrap{display:grid;grid-template-columns:280px minmax(0,1fr) 300px;gap:10px;
  padding:10px clamp(10px,2vw,20px) 14px;align-items:start}
.lower{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px;
  padding:0 clamp(10px,2vw,20px) 36px}
@media (max-width:1150px){.wrap{grid-template-columns:1fr}.lower{grid-template-columns:1fr}}
.cmp{display:grid;grid-template-columns:1fr auto;gap:0 8px;padding:5px 7px;
  border-left:2px solid transparent;border-bottom:1px solid #0b1f18;cursor:pointer;
  font:inherit;color:inherit;background:none;text-align:left;width:100%}
.cmp:hover{background:#08161180}
.cmp[aria-pressed="true"]{border-left-color:var(--p);background:#0a221a}
.cmp.warn[aria-pressed="true"]{border-left-color:var(--a)}
.cmp .n{letter-spacing:.04em}
.cmp .s{text-align:right;font-size:13px}
.cmp .m{grid-column:1/-1;font-size:9.5px;color:var(--faint);letter-spacing:.07em;text-transform:uppercase}
.feed div{padding:2px 0 2px 11px;position:relative;color:var(--muted);font-size:10.5px}
.feed div::before{content:"\00b7";position:absolute;left:0;color:var(--p)}
.feed .hd{color:var(--p);padding-left:0;margin-top:7px;font-size:9.5px;letter-spacing:.12em}
.feed .hd::before{content:""}
table{width:100%;border-collapse:collapse;font-size:10.5px}
td,th{padding:2.5px 4px;border-bottom:1px solid #0b1f18;text-align:left;vertical-align:top}
th{color:var(--faint);font-weight:400;font-size:9px;letter-spacing:.1em;text-transform:uppercase}
td.n,th.n{text-align:right}
.bar{height:4px;background:#0b1f18;border:1px solid #11261d;position:relative;min-width:40px}
.bar i{position:absolute;inset:0 auto 0 0;background:var(--p);display:block}
.bar i.am{background:var(--a)}
.ph{border:1px dashed var(--line);color:var(--faint);text-align:center;padding:14px 8px;font-size:10px;
  letter-spacing:.1em;text-transform:uppercase;
  background:repeating-linear-gradient(45deg,#060e0b 0 7px,#04090700 7px 14px)}
.btn{display:inline-block;border:1px solid var(--p-dim);color:var(--p);padding:3px 10px;border-radius:3px;
  font-size:10px;letter-spacing:.11em;text-transform:uppercase;cursor:pointer;background:#07140f;
  font-family:inherit}
.btn:hover{border-color:var(--p)}
.btn.a{border-color:var(--a-dim);color:var(--a)}
.btn.doc{border-color:var(--p);background:#0a2019}
.btn[disabled],.btn.off{border-color:#1a2a24;color:#2e4a3e;cursor:not-allowed}
.note{border-left:2px solid var(--a);padding:4px 8px;color:var(--a);font-size:10px;margin:6px 0}
.note.q{border-left-color:var(--p-dim);color:var(--muted)}
.kv{display:grid;grid-template-columns:auto 1fr;gap:1px 8px;font-size:10.5px}
.kv dt{color:var(--faint)} .kv dd{margin:0}
.fold{display:flex;align-items:center;gap:7px;width:100%;margin-top:7px;background:none;border:0;
  border-top:1px solid var(--line);padding:6px 0 0;cursor:pointer;font:inherit;color:var(--p-dim);
  font-size:9.5px;letter-spacing:.12em;text-transform:uppercase}
.fold:hover{color:var(--p)}
.fold .chev{display:inline-block;width:9px;transition:transform .12s}
.fold[aria-expanded="true"] .chev{transform:rotate(90deg)}
.fold .what{color:var(--faint)} .fold .cnt{margin-left:auto;color:var(--faint)}
.more[hidden]{display:none}
.more{border-left:1px solid var(--grid);margin:7px 0 0;padding:2px 0 2px 9px}
.more .u{color:var(--faint);font-size:9px;margin:8px 0 2px}
.scroll{max-height:320px;overflow:auto}
.plots{display:grid;grid-template-columns:1fr 1fr;gap:7px;margin-top:6px}
.plot{border:1px solid var(--line);background:#050c09;padding:4px 5px 2px}
.plot .ttl{font-size:9px;letter-spacing:.07em;color:var(--muted);display:flex;gap:6px}
.plot .ttl b{margin-left:auto;font-weight:400;color:var(--faint)}
.plot svg{display:block;width:100%;height:auto}
input[type=search]{font:inherit;font-size:11px;padding:3px 7px;border:1px solid var(--line);
  background:#050c09;color:var(--ink);border-radius:3px;width:100%;margin:4px 0}
.docbar{display:flex;gap:7px;flex-wrap:wrap;align-items:center;margin-top:8px}
#gl{position:relative;width:100%;height:330px;border:1px solid var(--line);background:#02050a;
  margin-top:6px;overflow:hidden}
#gl canvas{display:block}
.glmsg{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
  text-align:center;color:var(--faint);font-size:10.5px;letter-spacing:.08em;padding:12px}
.glbar{display:flex;gap:6px;flex-wrap:wrap;align-items:center;margin-top:5px}
.swatch{display:inline-block;width:8px;height:8px;border-radius:2px;margin-right:5px;
  vertical-align:middle}
.cmp .sw{display:inline-block;width:7px;height:7px;border-radius:2px;margin-right:5px}
.copied{color:var(--p);font-size:9.5px;letter-spacing:.1em}
</style>
</head>
<body>
<header class="band">
  <div class="row1" id="brow"></div>
  <div class="readouts" id="readouts"></div>
</header>
<div class="wrap">
  <section class="frame"><header><span class="t">compound set</span><span class="r" id="setr"></span></header>
    <div class="body" style="padding:0" id="set"></div></section>
  <section class="frame"><header><span class="t" id="pockt">binding site</span><span class="r" id="pockr"></span></header>
    <div class="body" id="pocket"></div></section>
  <section class="frame"><header><span class="t">source / data feed</span><span class="r">append only</span></header>
    <div class="body feed" id="feed"></div></section>
</div>
<div class="lower">
  <section class="frame"><header><span class="t">context &middot; phase 1</span><span class="r" id="ctxr"></span></header>
    <div class="body" id="ctx"></div></section>
  <section class="frame"><header><span class="t">experiment &middot; phase 3</span><span class="r" id="expr"></span></header>
    <div class="body" id="exp"></div></section>
  <section class="frame"><header><span class="t">dynamics &middot; phase 3</span><span class="r" id="mdr"></span></header>
    <div class="body" id="md"></div></section>
</div>
<script>/*__VIEWER_LIB__*/</script>
<script>
const DATA = /*__DATA__*/null;
const C = DATA.context, X = DATA.experiments;
const $ = s => document.querySelector(s);
const el = (t,c,x) => { const n=document.createElement(t); if(c) n.className=c;
  if(x!==undefined&&x!==null) n.textContent=x; return n; };
const num = (v,d=2) => (v===null||v===undefined||v==="") ? null : Number(v).toFixed(d);
const NA = "not applicable", NR = "no record";
function kv(pairs){ const d=el("dl","kv");
  for(const [k,v] of pairs){ if(v===null||v===undefined||v==="") continue;
    d.append(el("dt",null,k)); const dd=el("dd"); if(v instanceof Node) dd.append(v); else dd.textContent=v;
    d.append(dd); } return d; }
function link(href,text){ const a=el("a",null,text); a.href=href; a.target="_blank"; a.rel="noopener"; return a; }
function fold(host, label, what, count, build){
  const b=el("button","fold"); b.type="button"; b.setAttribute("aria-expanded","false");
  b.append(el("span","chev","▶"), el("span",null,label), el("span","what",what||""),
           el("span","cnt",count===null||count===undefined?"":String(count)));
  const more=el("div","more"); more.hidden=true;
  b.addEventListener("click",()=>{ const open=b.getAttribute("aria-expanded")==="true";
    b.setAttribute("aria-expanded",open?"false":"true"); more.hidden=open;
    if(!open && !more.dataset.built){ build(more); more.dataset.built="1"; } });
  host.append(b,more);
}
function table(head, rows, opts={}){
  const t=el("table"); if(head){ const h=el("thead"), tr=el("tr");
    head.forEach(x=>tr.append(el("th",(opts.num||[]).includes(x)?"n":null,x))); h.append(tr); t.append(h); }
  const b=el("tbody");
  rows.forEach(r=>{ const tr=el("tr");
    r.forEach(c=>{ const td=el("td"); if(c instanceof Node) td.append(c);
      else if(c&&c.n!==undefined){ td.className="n"; td.textContent=c.n; }
      else td.textContent=(c===null||c===undefined)?"":c; tr.append(td); });
    b.append(tr); });
  t.append(b); return t;
}

/* ------------------------------------------------------------------ band */
const compounds = ((X&&X.runs)||[]).filter(r=>r.kind==="compound");
const comparable = compounds.filter(r=>!r.box_limited && r.score!==null && r.score!==undefined);
const scores = comparable.map(r=>r.score);
const anchorRun = compounds.find(r=>r.anchor);
const mdStages = ((X&&X.md)||[]).reduce((n,p)=>n+(p.stages||[]).length,0);

const brow=$("#brow");
brow.append(el("span","logo","Protein-MCP · Target console"));
if(C.identity.gene) brow.append(Object.assign(el("span","pill"),{textContent:
  (C.identity.gene||"")+(X&&X.mutation?" "+X.mutation:"")}));
brow.append(Object.assign(el("span","pill p"),{textContent:DATA.accession}));
if(X&&X.choice&&X.choice.structure) brow.append(Object.assign(el("span","pill d"),
  {textContent:X.choice.structure+(X.pockets&&X.pockets[0]?" · pocket "+X.pockets[0].key:"")}));
const sub=el("span","sub","built "+DATA.built_at+" · no network · quoted from the dossier · "
  +(DATA.viewer_credit||"no 3D viewer embedded"));
sub.style.marginLeft="auto"; brow.append(sub);

function readout(value, label, cls, title){
  const d=el("div","ro"); if(title) d.title=title;
  d.append(Object.assign(el("b",cls||null),{textContent:value}), el("span",null,label)); return d;
}
const ro=$("#readouts");
if(X){
  ro.append(readout(compounds.length+" / "+((X.compound_choice.compounds||[]).length||compounds.length),
    "docked of selected","g"));
  if(scores.length){
    ro.append(readout(num(Math.min(...scores))+" … "+num(Math.max(...scores)),
      "score range kcal/mol",null,"the recorded top-pose scores of the comparable runs"));
    const span=Math.max(...scores)-Math.min(...scores);
    /* the span is arithmetic on the recorded scores and nothing else. There is NO threshold
       here: this project has measured no interpretable difference, and printing one would be
       a claim with no source behind it. */
    ro.append(readout(num(span)+" kcal/mol","span of the set",null,
      "the largest recorded score minus the smallest, over the comparable runs. No threshold "
      +"of an interpretable difference is pre-registered, so this span orders nothing"));
  }
  if(anchorRun){ const a=anchorRun.anchor;
    ro.append(readout(num(a.rmsd,4)+" ±"+num(a.bar,3),
      a.decides===false ? "anchor Å · does not decide" : "anchor Å",
      a.decides===false?"am":"g"));
  } else ro.append(readout("none","anchor","fa","no T1 compound in this set"));
  const limited=compounds.filter(r=>r.box_limited).length;
  if(limited) ro.append(readout(String(limited),"excluded from ranking","am"));
  const admetOfSet = (X.admet||[]).filter(a=>compounds.some(
    r=>r.id && (a.pocket||"").toUpperCase().endsWith(r.id.toUpperCase()))).length;
  ro.append(readout(admetOfSet+" / "+compounds.length,"admet of the set",
    admetOfSet>=compounds.length?"g":"am",
    (X.admet||[]).length+" ADMET results in all, including crystal ligands of other pockets"));
  ro.append(readout(String((X.md||[]).length),"md packages",(X.md||[]).length?"am":"fa"));
  ro.append(readout(String(mdStages),"trajectories measured",mdStages?"g":"fa"));
}
ro.append(readout((C.retrieved.profile||"").slice(0,10)||NR,"sources retrieved","fa"));

/* ------------------------------------------------------------------ compound set */
let selected = anchorRun || compounds[0] || null;
function drawSet(){
  const host=$("#set"); host.replaceChildren();
  $("#setr").textContent = compounds.length+" docked · "+(C.ligands.rows||[]).length+" retrieved";
  if(!compounds.length){ host.append(Object.assign(el("div","ph"),
    {textContent:"no compound docked — run select_ligands and dock_ligand"})); return; }
  compounds.forEach(r=>{
    const b=el("button","cmp"+(r.box_limited?" warn":"")); b.type="button";
    b.setAttribute("aria-pressed", selected && r.run===selected.run ? "true":"false");
    const sw = el("span","sw"); sw.style.background = poseColour(r.run);
    const nm = el("span","n"+(r.anchor?" g":"")); nm.append(sw, document.createTextNode(r.name||r.id||r.run));
    b.append(nm,
             Object.assign(el("span","s"+(r.box_limited?" am":"")),{textContent:num(r.score,3)??NR}));
    const bits=[(r.id||"").toLowerCase(), r.level||"", "le "+(num(r.le,4)??NR)];
    if(r.anchor) bits.push("anchor "+num(r.anchor.rmsd,4));
    if(r.box_limited) bits.push("box limited "+num(r.clearance)+" Å · not ranked");
    b.append(Object.assign(el("span","m"),{textContent:bits.filter(Boolean).join(" · ")}));
    b.addEventListener("click",()=>{ selected=r; drawSet(); drawLower(); });
    host.append(b);
  });
  const foot=el("div"); foot.style.padding="7px 8px";
  if(scores.length>1){
    const span=Math.max(...scores)-Math.min(...scores);
    foot.append(Object.assign(el("div","note"),{textContent:
      "These runs span "+num(span)+" kcal/mol. No threshold for an interpretable difference "+
      "is pre-registered in this project, because none has been measured here, so this page "+
      "says nothing about which of them is better. The order is the order the manifest holds, "+
      "nothing more, and the set has no inactive compound in it, so it cannot show that the "+
      "score separates actives from inactives either."}));
  }
  host.append(foot);
}

/* ------------------------------------------------------------------ pocket */
function pymolLine(what){
  const s = X&&X.choice&&X.choice.structure ? X.choice.structure : "<structure>";
  return 'open_in_pymol("'+DATA.accession+'", structure="'+s+'")  '+(what?"# "+what:"");
}
/* the calls these two buttons copy name a run that IS in the record: the first measured
   stage, taken from the manifest, never a structure or a stage typed in here */
function measuredStage(){
  for(const m of (X&&X.md)||[]) for(const s of m.stages||[])
    if(s.label) return {pdb:m.pdb, pocket:m.pocket, stage:s.label};
  return null;
}
function trajectoryLine(){
  const r=measuredStage(); if(!r) return "";
  return 'open_md_in_pymol("'+DATA.accession+'", structure="'+r.pdb+'", pocket="'+r.pocket+
         '", stage="'+r.stage+'")  # plays it in a PyMOL window on this machine';
}
function movieLine(){
  const r=measuredStage(); if(!r) return "";
  return 'render_md_movie("'+DATA.accession+'", structure="'+r.pdb+'", pocket="'+r.pocket+
         '", stage="'+r.stage+'")  # writes an .mp4 beside the reports';
}
function copyBtn(label, text, cls){
  const b=el("button","btn"+(cls?" "+cls:"")); b.type="button"; b.textContent=label;
  b.addEventListener("click", async () => {
    try { await navigator.clipboard.writeText(text); } catch (e) { /* file:// may refuse */ }
    const old=b.textContent; b.textContent="copied ✓";
    setTimeout(()=>{ b.textContent=old; }, 1400);
  });
  b.title=text; return b;
}
function drawPocket(){
  const host=$("#pocket"); host.replaceChildren();
  /* the pocket of the SELECTED run, matched on structure as well as key: two structures can
     hold the same pocket key, and drawing the box of the wrong one would be a believable lie */
  const pk=(X&&X.pockets||[]).find(p=>selected && p.key===selected.pocket && p.pdb===selected.pdb)
           || (selected ? null : (X&&X.pockets||[])[0]);
  $("#pockt").textContent = pk ? "binding site · "+pk.pdb+" pocket "+pk.key : "binding site";
  $("#pockr").textContent = selected ? (selected.name||selected.run)+" selected" : "";
  if(pk){
    const head=el("div"); head.style.display="flex"; head.style.gap="7px"; head.style.flexWrap="wrap";
    head.append(Object.assign(el("span","pill d"),{textContent:"box "+
      (pk.box||[]).map(x=>Number(x).toFixed(2)).join(" × ")+" Å"}));
    head.append(Object.assign(el("span","pill d"),{textContent:(pk.residues||[]).length+" residues ≤ 5 Å"}));
    if(pk.large_box) head.append(Object.assign(el("span","pill"),{textContent:"large box"}));
    if(pk.spans) head.append(Object.assign(el("span","pill"),{textContent:"spans chains"}));
    host.append(head);
    host.append(viewerPanel(pk));
    const lst=el("div","fa"); lst.style.fontSize="9.5px"; lst.style.marginTop="4px";
    lst.textContent = "closest residues, "+Math.min(10,(pk.residues||[]).length)+" of "+
      (pk.residues||[]).length+": " + (pk.residues||[]).slice(0,10).map((r,i)=>
      r+" "+Number((pk.distances||[])[i]).toFixed(2)+" Å").join(" · ");
    host.append(lst);
    fold(host,"flat schematic","the box to scale, no rotation",null,
         m=>m.append(pocketSvg(pk, selected)));
  } else host.append(Object.assign(el("div","ph"),{textContent: selected
    ? "no pocket of "+selected.pdb+" with the key "+selected.pocket+" is in the manifest"
    : "no pocket in the manifest"}));

  const bar=el("div","docbar");
  DATA.documents.forEach(d=>{
    if(d.exists){ const a=el("a","btn doc","▤ "+d.label); a.href=d.file; a.target="_blank"; bar.append(a); }
    else { const b=el("span","btn off"); b.textContent="▤ "+d.label+" — not built yet"; bar.append(b); }
  });
  host.append(bar);
  const bar2=el("div","docbar");
  bar2.append(copyBtn("◰ copy pymol call", pymolLine("opens the receptor and the pocket")));
  const haveTraj = mdStages>0;
  const movie=el("span","btn"+(haveTraj?" a":" off"));
  if(haveTraj){
    bar2.append(copyBtn("▶ copy trajectory call", trajectoryLine()));
    bar2.append(copyBtn("⏺ copy movie call", movieLine()));
  } else {
    const movie=el("span","btn off");
    movie.textContent="▶ trajectory — none measured yet";
    bar2.append(movie);
  }
  host.append(bar2);
  host.append(Object.assign(el("div","note q"),{textContent:
    "A page cannot start a program on your machine. The buttons above copy the call to run. "+
    "The 3D view above is a viewer embedded in this file; the trajectory itself is not in this "+
    "page (a DCD is hundreds of megabytes), so open_md_in_pymol plays it and render_md_movie "+
    "writes it as an .mp4 beside the reports."}));
}
// one colour per compound, reused by the list on the left and by the 3D view so the same
// molecule is the same colour everywhere
const POSE_COLOURS = ["#35e39a","#ffab3d","#6fb7ff","#ff8a84","#d79bff","#ffd84d","#8be36a",
                      "#9ad0c2","#ff6fae","#c9d94d"];
function poseColour(run){
  const list = ((DATA.three_d||{}).poses)||[];
  const i = list.findIndex(x=>x.run===run);
  return POSE_COLOURS[(i<0?0:i) % POSE_COLOURS.length];
}
let GLVIEW = null, GLSHOWN = "selected";
function viewerPanel(pk){
  const wrap = el("div");
  const host = el("div"); host.id = "gl";
  wrap.append(host);
  const td = DATA.three_d || {};
  const bar = el("div","glbar");
  const mk = (label, mode) => { const b=el("button","btn"+(GLSHOWN===mode?" doc":"")); b.type="button";
    b.textContent=label; b.addEventListener("click",()=>{ GLSHOWN=mode; drawPocket(); }); return b; };
  bar.append(mk("selected pose only","selected"), mk("all poses","all"), mk("pocket residues","pocket"));
  const reset=el("button","btn"); reset.type="button"; reset.textContent="recentre";
  reset.addEventListener("click",()=>{ if(GLVIEW){ GLVIEW.zoomTo({model:-1}); GLVIEW.render(); } });
  bar.append(reset);
  wrap.append(bar);
  const legend = el("div","fa"); legend.style.fontSize="9.5px"; legend.style.marginTop="4px";
  wrap.append(legend);

  // the library is embedded, but if it ever is not the panel says so instead of staying blank
  if (typeof $3Dmol === "undefined"){
    host.append(Object.assign(el("div","glmsg"),{textContent:
      "3D viewer not available: the library was not embedded in this page"}));
    return wrap;
  }
  const receptorText = (td.receptors||{})[pk && pk.pdb];
  const poses = (td.poses||[]).filter(x=>!pk || (x.pdb===pk.pdb && x.pocket===pk.key));
  if(!receptorText && !poses.length){
    host.append(Object.assign(el("div","glmsg"),{textContent:
      "nothing to show in 3D for this pocket: no receptor and no pose were embedded"}));
    return wrap;
  }
  // 3Dmol needs the element to be in the document before it builds its canvas
  setTimeout(()=>{
    try{
      const v = $3Dmol.createViewer(host, {backgroundColor:"#02050a"});
      GLVIEW = v;
      const shown = [];
      if(receptorText){
        v.addModel(receptorText, "pdb");
        v.setStyle({}, {cartoon:{color:"#3f5d57", opacity:0.9}});
        shown.push("receptor "+pk.pdb+" as cartoon");
      }
      const wanted = GLSHOWN==="all" ? poses
                   : poses.filter(x=>selected && x.run===selected.run);
      wanted.forEach(x=>{
        const col = poseColour(x.run);
        const mdl = v.addModel(x.sdf, "sdf");
        // one flat colour per compound: it matches the swatch in the list on the left, and with
        // six ligands superposed in one pocket telling them apart matters more than element colours
        mdl.setStyle({}, {stick:{radius:0.18, color: col}});
        shown.push(x.name);
      });
      if(GLSHOWN==="pocket" && receptorText && (pk.residues||[]).length){
        const nums = (pk.residues||[]).map(r=>parseInt(String(r).replace(/[^0-9]/g,""),10))
                                      .filter(n=>!isNaN(n));
        v.setStyle({resi:nums}, {stick:{radius:0.1, colorscheme:"grayCarbon"},
                                 cartoon:{color:"#3f5d57", opacity:0.9}});
        shown.push((pk.residues||[]).length+" pocket residues as sticks");
      }
      if(wanted.length) v.zoomTo({model:-1}); else v.zoomTo();
      v.render();
      legend.textContent = "showing: " + shown.join(" \u00b7 ") +
        "  \u2014  drag to rotate, scroll to zoom, right-drag or shift-drag to pan";
    }catch(err){
      host.replaceChildren(Object.assign(el("div","glmsg"),
        {textContent:"3D viewer failed to start: "+err}));
    }
  }, 0);
  return wrap;
}

function pocketSvg(pk, run){
  const NS="http://www.w3.org/2000/svg";
  const s=document.createElementNS(NS,"svg"); s.setAttribute("viewBox","0 0 640 320");
  s.style.width="100%"; s.style.height="auto"; s.style.display="block"; s.style.marginTop="6px";
  const mk=(t,at)=>{ const n=document.createElementNS(NS,t);
    for(const k in at) n.setAttribute(k,at[k]); return n; };
  const box=pk.box||[20,20,20];
  const scale=240/Math.max(...box.map(Number));
  const w=Number(box[0])*scale, h=Number(box[2])*scale, cx=320, cy=150;
  s.append(mk("rect",{x:cx-w/2,y:cy-h/2,width:w,height:h,fill:"none",stroke:"#1d6b4e","stroke-dasharray":"5 4"}));
  s.append(mk("circle",{cx:cx,cy:cy,r:Math.min(w,h)*0.42,fill:"none",stroke:"#123b2c"}));
  const n=Math.min((pk.residues||[]).length,14);
  for(let i=0;i<n;i++){
    const ang=(i/n)*2*Math.PI, rr=Math.min(w,h)*0.52;
    const t=mk("text",{x:cx+Math.cos(ang)*rr,y:cy+Math.sin(ang)*rr,fill:"#2f7a5e",
      "font-size":8.5,"font-family":"monospace","text-anchor":"middle"});
    t.textContent=(pk.residues||[])[i]; s.append(t);
  }
  const lbl=(x,y,txt,col)=>{ const t=mk("text",{x:x,y:y,fill:col,"font-size":9,"font-family":"monospace"});
    t.textContent=txt; s.append(t); };
  lbl(cx-w/2, cy-h/2-8, "SEARCH BOX · "+box.map(x=>Number(x).toFixed(2)).join(" × ")+" Å", "#38604f");

  // the pose itself, from the coordinates the run wrote. x across the box's first axis and
  // z down its third, the same two the box rectangle uses, so atom and box share one scale.
  const geo = run && run.geometry;
  if(geo && (geo.atoms||[]).length){
    const px = a => cx + a.x*scale, py = a => cy - a.z*scale;
    (geo.bonds||[]).forEach(([i,j])=>{
      const a=geo.atoms[i], b=geo.atoms[j]; if(!a||!b) return;
      s.append(mk("line",{x1:px(a),y1:py(a),x2:px(b),y2:py(b),
        stroke:"#35e39a","stroke-width":1.6,"stroke-linecap":"round"}));
    });
    geo.atoms.forEach(a=>{
      if(a.e==="C") return;                       // carbons are the skeleton, drawn as bonds only
      const col = a.e==="N" ? "#6fb7ff" : a.e==="O" ? "#ff8a84"
                : a.e==="S" ? "#ffd84d" : a.e==="F"||a.e==="Cl" ? "#9be36a" : "#bff3da";
      s.append(mk("circle",{cx:px(a),cy:py(a),r:3.1,fill:"#04070a",stroke:col,"stroke-width":1.2}));
      const t=mk("text",{x:px(a),y:py(a)+2.6,fill:col,"font-size":5.5,"font-family":"monospace",
        "text-anchor":"middle"}); t.textContent=a.e; s.append(t);
    });
    lbl(cx-w/2, cy+h/2+18, (run.name||run.run)+" · top pose · "+
      geo.atoms.length+" heavy atoms · projected on the box's 1st and 3rd axes", "#35e39a");
  } else if(run){
    lbl(cx-w/2, cy+h/2+18, (run.name||run.run)+" · pose not drawn: "+
      ((geo&&geo.error)||"no pose file recorded"), "#8a5c20");
  }
  if(run){
    if(run.anchor) lbl(cx-w/2, cy+h/2+31,
      "its reference pose in "+run.anchor.reference+" is Δ "+num(run.anchor.rmsd,4)+
      " Å away ±"+num(run.anchor.bar,3)+" (measured by cross_dock_rmsd, not drawn here)",
      "#ffab3d");
    else lbl(cx-w/2, cy+h/2+31, "no reference pose exists anywhere — RMSD "+NA, "#38604f");
  }
  return s;
}

/* ------------------------------------------------------------------ feed */
function drawFeed(){
  const f=$("#feed"); f.replaceChildren();
  const line=(t,cls)=>f.append(Object.assign(el("div",cls||null),{textContent:t}));
  const head=t=>f.append(Object.assign(el("div","hd u"),{textContent:t}));
  head("RETRIEVED");
  line("UniProt "+DATA.accession+" · "+(C.retrieved.profile||NR)+" · "+
       (C.variants||[]).length+" natural variants");
  if(C.assoc.total!==null) line("Open Targets "+(C.identity.ot_release||"")+" · "+
       C.assoc.total+" associations · "+C.assoc.examined+" examined");
  if(C.structures.total!==null) line("RCSB · "+C.structures.total+" found · "+
       (C.structures.confirmed??NR)+" confirmed · "+(C.structures.rows||[]).length+" offered");
  if(C.ligands.qualifying!==null&&C.ligands.qualifying!==undefined)
    line("ChEMBL "+(C.ligands.target||"")+" · "+(C.ligands.lower_bound?"≥":"")+
         C.ligands.qualifying+" ligands · "+((C.ligands.activities||{}).fetched??NR)+" activities read");
  if(C.literature.total!==null) line("Europe PMC · \""+(C.literature.query||"")+"\" · "+
       C.literature.total+" hits · "+(C.literature.rows||[]).length+" kept");
  if(X){
    head("PREPARED");
    if(X.choice.structure) line(X.choice.structure+" chosen by the researcher"+
      (X.choice.resolution?" · "+X.choice.resolution+" Å":""));
    if(X.choice.reason) line("reason: “"+X.choice.reason+"”");
    (X.pockets||[]).forEach(p=>line("pocket "+p.key+" · "+(p.residues||[]).length+
      " residues · box "+(p.box||[]).map(x=>Number(x).toFixed(2)).join("×")));
    head("DOCKED");
    line((X.runs||[]).length+" runs recorded in "+
         new Set((X.runs||[]).map(r=>r.pdb+" "+r.pocket)).size+" pocket(s)");
    (X.pockets||[]).forEach(p=>line("box of "+p.key+": "+
         (p.box||[]).map(x=>Number(x).toFixed(2)).join("×")+" Å"));
    (X.runs||[]).filter(r=>r.box_limited).forEach(r=>line(
      (r.name||r.run)+" box_limited · clearance "+num(r.clearance)+" Å on the "+
      num(r.limiting_axis)+" Å axis", "am"));
    if(anchorRun){ const a=anchorRun.anchor;
      head("ANCHORED");
      line("superposition "+num(a.sup_rmsd)+" Å over "+a.sup_pairs+" CA pairs");
      line("reference matched by "+(a.matched_by||NR));
      line("bar ±"+num(a.bar,3)+" Å · criterion "+num(a.criterion)+
           " Å · decides: "+String(a.decides), a.decides===false?"am":null);
    }
    head("NOT MEASURED");
    const noRmsd=compounds.filter(r=>!r.anchor).length;
    if(noRmsd) line("RMSD for "+noRmsd+" compound(s): "+NA+", no reference pose exists anywhere","fa");
    line("predicted ADMET: none — this project holds no trained model","fa");
    if(!mdStages) line("molecular dynamics: "+(X.md||[]).length+" package(s) written, 0 trajectories measured","fa");
  }
  fold(f,"full feed","every count and every refusal",null,m=>{
    m.append(Object.assign(el("div","u"),{textContent:"SELECTION, AS THE TOOLS REPORTED IT"}));
    [C.structures.selection, C.ligands.selection, C.assoc.selection]
      .filter(Boolean).forEach(t=>m.append(Object.assign(el("div","fa"),
        {textContent:t, style:"font-size:10px;margin:3px 0"})));
    if(X&&X.compound_choice.reason){
      m.append(Object.assign(el("div","u"),{textContent:"WHY THESE COMPOUNDS"}));
      m.append(Object.assign(el("div","fa"),{textContent:X.compound_choice.reason,
        style:"font-size:10px"}));
    }
    if(X){ m.append(Object.assign(el("div","u"),{textContent:"SCOPE OF EACH RUN"}));
      (X.runs||[]).forEach(r=>{ if(r.scope) m.append(Object.assign(el("div","fa"),
        {textContent:(r.name||r.run)+": "+r.scope, style:"font-size:9.5px;margin:3px 0"})); }); }
  });
}

/* ------------------------------------------------------------------ lower band */
function ligandOf(run){
  if(!run) return null;
  const want=(run.id||"").toUpperCase(), key=(run.inchikey||"").toUpperCase();
  return (C.ligands.rows||[]).find(l=>(l.id||"").toUpperCase()===want
    || (key && (l.inchikey||"").toUpperCase()===key)) || null;
}
function drawContext(){
  const host=$("#ctx"); host.replaceChildren();
  const lig=ligandOf(selected);
  $("#ctxr").textContent = selected ? (selected.name||selected.run) : "";
  if(!lig){
    host.append(Object.assign(el("div","ph"),{textContent:
      selected ? "this compound is not among the ligands retrieved for this target" : "no compound selected"}));
    host.append(Object.assign(el("div","note"),{textContent:
      "That is a fact about the retrieval, not about the compound: the ChEMBL download was "+
      "truncated, so a compound can be absent from it and still have measurements."}));
  } else {
    host.append(kv([
      ["CHEMBL", lig.link?link(lig.link,lig.id):lig.id],
      ["STATUS", lig.phase>=4 ? "approved drug"+(lig.approval?" · first approved "+lig.approval:"")
                              : (lig.phase?"clinical phase "+lig.phase:"no clinical phase recorded")],
      ["MECHANISM", lig.moa||NR],
      ["BEST pCHEMBL", lig.pchembl!==null&&lig.pchembl!==undefined?String(lig.pchembl):NR],
      ["MEASUREMENTS", lig.n!==null&&lig.n!==undefined?String(lig.n):NR],
    ]));
    if((lig.potency||[]).length){
      host.append(Object.assign(el("div","u fa"),{textContent:"MEASURED AGAINST THIS TARGET",
        style:"margin-top:8px;font-size:9px"}));
      host.append(table(["mut","type · assay","best","median","worst","n"],
        lig.potency.map(p=>[p.mutation||"any", (p.type||"")+" · "+(p.assay||""),
          {n:p.best??""},{n:p.median??""},{n:p.worst??""},{n:p.n??""}]),
        {num:["best","median","worst","n"]}));
      host.append(Object.assign(el("div","fa"),{textContent:"nM · values as ChEMBL reports them",
        style:"font-size:9.5px"}));
      host.append(Object.assign(el("div","note"),{textContent:
        "MEASURED values. The spread between best and worst is the spread of the LITERATURE, not "+
        "an error bar of the compound, and none of it is comparable with the docking score beside "+
        "it: different quantity, different units, different experiment."}));
    }
  }
  fold(host,"all retrieved context","associations · variants · structures · papers",
       (C.assoc.rows||[]).length+(C.variants||[]).length, m=>{
    m.append(Object.assign(el("div","u"),{textContent:
      "TARGET–DISEASE ASSOCIATIONS · "+(C.assoc.rows||[]).length+" of "+C.assoc.total}));
    const q=el("input"); q.type="search"; q.placeholder="filter diseases";
    const shown=el("div","note"); const box=el("div","scroll"); m.append(q,shown,box);
    const draw=()=>{ const t=q.value.trim().toLowerCase();
      const all=(C.assoc.rows||[]).filter(r=>!t||(r.name||"").toLowerCase().includes(t));
      const rows=all.slice(0,120);
      shown.textContent = rows.length<all.length
        ? "showing the first "+rows.length+" of "+all.length+" rows that match; filter to see more"
        : rows.length+" of "+(C.assoc.rows||[]).length+" rows shown";
      box.replaceChildren(table(["disease","score",""],rows.map(r=>{
        const bar=el("div","bar"); const i=el("i"); i.style.width=Math.round((r.score||0)*100)+"%";
        bar.append(i);
        return [r.name+(r.onc?" · onc":""), {n:num(r.score,3)}, bar]; }),{num:["score"]})); };
    q.addEventListener("input",draw); draw();
    m.append(Object.assign(el("div","u"),{textContent:"NATURAL VARIANTS · "+(C.variants||[]).length}));
    const vbox=el("div","scroll");
    vbox.append(table(["change","pos","uniprot description","pdb","compounds","papers"],
      (C.variants||[]).map(v=>[v.change, {n:v.pos??""}, v.description,
        {n:v.pdb??"—"},{n:v.chembl??"—"},{n:v.papers??"—"}]),
      {num:["pos","pdb","compounds","papers"]}));
    m.append(vbox);
    m.append(Object.assign(el("div","note"),{textContent:
      "The last three columns are empty on purpose: counting structures, compounds and papers per "+
      "variant needs one query per variant per source, and those were never run for any variant "+
      "except the one this dossier was filtered by. A number here would be invented."}));
    m.append(Object.assign(el("div","u"),{textContent:"STRUCTURES OFFERED · "+(C.structures.rows||[]).length+
      " of "+C.structures.total}));
    m.append(table(["pdb","resolution","ligands","title"],(C.structures.rows||[]).map(s=>[
      s.rcsb?link(s.rcsb,s.pdb):s.pdb, {n:s.res!==null?s.res+" Å":NR},
      (s.ligands||[]).join(", "), s.title]),{num:["resolution"]}));
    m.append(Object.assign(el("div","u"),{textContent:"PAPERS · "+(C.literature.rows||[]).length+
      " of "+C.literature.total}));
    (C.literature.rows||[]).forEach(p=>{
      const d=el("div"); d.style.margin="4px 0";
      d.append(Object.assign(el("div",null),{textContent:p.title}));
      const meta=el("div","fa"); meta.style.fontSize="9.5px";
      meta.textContent=[p.journal,p.year,p.oa?"open access":null].filter(Boolean).join(" · ");
      if(p.pmid){ meta.append(" "); meta.append(link("https://pubmed.ncbi.nlm.nih.gov/"+p.pmid+"/","PMID "+p.pmid)); }
      d.append(meta); m.append(d);
    });
  });
}
function drawExperiment(){
  const host=$("#exp"); host.replaceChildren();
  const r=selected;
  $("#expr").textContent = r ? r.run : "";
  if(!r){ host.append(Object.assign(el("div","ph"),{textContent:"no docking run recorded"})); return; }
  const rows=[["vina score",{n:num(r.score,3)??NR}],
              ["ligand efficiency",{n:num(r.le,4)??NR}],
              ["heavy atoms",{n:r.heavy??NR}],
              ["poses returned",{n:(r.poses||[]).length}]];
  if(r.anchor) rows.push(["anchor rmsd",{n:num(r.anchor.rmsd,4)}]);
  host.append(table(null,rows));
  if(r.anchor){
    const a=r.anchor;
    host.append(Object.assign(el("div","note"),{textContent:
      a.decides===false
        ? "±"+num(a.bar,3)+" Å puts the RMSD at "+num(a.range[0])+"–"+num(a.range[1])+
          " Å, both sides of the "+num(a.criterion)+" Å criterion. THIS RUN DOES NOT DECIDE "+
          "whether the pose was reproduced. It is not a pass."
        : "With the ±"+num(a.bar,3)+" Å bar the measurement stays on one side of the "+
          num(a.criterion)+" Å criterion, so it decides."}));
  } else {
    host.append(Object.assign(el("div","note"),{textContent:
      (r.rmsd_na||"No RMSD: no reference pose exists for this compound.")+" This is "+NA+
      ", not a missing value. The score is a hypothesis about where the compound would sit."}));
  }
  if(r.box_limited) host.append(Object.assign(el("div","note"),{textContent:
    "box_limited: "+num(r.clearance)+" Å of clearance on the "+num(r.limiting_axis)+
    " Å axis (compound extent "+num(r.extent)+" Å). Vina never explored the orientations it "+
    "explored for the others, so this score is not comparable with theirs and is not ranked."}));
  const ad=(X.admet||[]).find(a=>selected && (a.pocket||"").toUpperCase().endsWith((r.id||"\u0000").toUpperCase()));
  if(ad){
    host.append(Object.assign(el("div","u fa"),{textContent:"ADMET · CALCULATED + RETRIEVED, NEVER PREDICTED",
      style:"margin-top:8px;font-size:9px"}));
    const p=ad.props||{};
    const viol=Object.entries(ad.rules||{}).map(([n,v])=>n+" "+(v.violations??"?"));
    host.append(table(null,[
      ["MW / cLogP",{n:(num(p.molecular_weight,3)??NR)+" / "+(num(p.clogp,3)??NR)}],
      ["TPSA / rot. bonds",{n:(num(p.tpsa,2)??NR)+" / "+(p.rotatable_bonds??NR)}],
      ["rule violations",{n:viol.join(" · ")}],
      ["measured ADME rows",{n:ad.retrieved?(ad.adme_rows||[]).length+" of "+(ad.rows_read??"?")+" read":"UNKNOWN"}],
      ["predicted",{n:"none — by design"}],
    ]));
  }
  fold(host,"all poses · anchor · admet rows","nothing summarised away",
       (r.poses||[]).length, m=>{
    m.append(Object.assign(el("div","u"),{textContent:"THE POSES AS VINA RETURNED THEM"}));
    m.append(table(["mode","score",r.anchor?"rmsd to reference":"rmsd"],
      (r.poses||[]).map(p=>[{n:p.mode},{n:num(p.score,3)},{n:p.rmsd!==null&&p.rmsd!==undefined?num(p.rmsd,4):NA}]),
      {num:["mode","score","rmsd to reference","rmsd"]}));
    if(r.anchor){ const a=r.anchor;
      m.append(Object.assign(el("div","u"),{textContent:"HOW THE ANCHOR WAS MEASURED"}));
      m.append(Object.assign(el("div","fa"),{style:"font-size:10px",textContent:
        "reference pose from "+a.reference+", matched by "+a.matched_by+" · proteins superposed "+
        "with "+a.sup_pairs+" CA pairs at "+num(a.sup_rmsd)+" Å · the pocket-restricted "+
        "superposition moves that pose "+num(a.bar,3)+" Å, which is the bar"}));
    }
    if(ad && (ad.adme_rows||[]).length){
      m.append(Object.assign(el("div","u"),{textContent:"MEASURED ADME ROWS, UNCONVERTED"}));
      m.append(table(["property","value","units","assay"], ad.adme_rows.map(x=>[
        x.cls||x.type, ((x.relation||"")+" "+(x.value??"")).trim(), x.units||"—", x.assay||""])));
      m.append(Object.assign(el("div","fa"),{style:"font-size:10px",textContent:
        "Values keep their own units and their relation; nothing is averaged or converted."}));
    }
    if(r.scope){ m.append(Object.assign(el("div","u"),{textContent:"SCOPE"}));
      m.append(Object.assign(el("div","fa"),{style:"font-size:10px",textContent:r.scope})); }
  });
}
function plot(title, line, colour, dashed, extra){
  const NS="http://www.w3.org/2000/svg";
  const d=el("div","plot");
  const t=el("div","ttl"); t.append(el("span",null,title),
    Object.assign(el("b"),{textContent:num(line.min,2)+" … "+num(line.max,2)}));
  d.append(t);
  const W=250,H=46,P=3;
  const s=document.createElementNS(NS,"svg"); s.setAttribute("viewBox","0 0 "+W+" "+H);
  s.setAttribute("preserveAspectRatio","none");
  const frame=document.createElementNS(NS,"rect");
  frame.setAttribute("width",W); frame.setAttribute("height",H);
  frame.setAttribute("fill","none"); frame.setAttribute("stroke","#0e2b22"); s.append(frame);
  const add=(ln,col,dash)=>{
    const lo=Math.min(...ln.v), hi=Math.max(...ln.v), rng=(hi-lo)||1, tmax=ln.t[ln.t.length-1]||1;
    const pts=ln.v.map((v,i)=>{
      const x=P+(ln.t[i]/tmax)*(W-2*P), y=H-P-((v-lo)/rng)*(H-2*P);
      return x.toFixed(1)+","+y.toFixed(1); }).join(" ");
    const pl=document.createElementNS(NS,"polyline");
    pl.setAttribute("points",pts); pl.setAttribute("fill","none");
    pl.setAttribute("stroke",col); pl.setAttribute("stroke-width","1.1");
    if(dash) pl.setAttribute("stroke-dasharray","3 2");
    s.append(pl);
  };
  add(line,colour,dashed);
  if(extra) add(extra,colour,true);
  d.append(s);
  return d;
}
function poseRmsd(mine){
  /* three different things that must not be printed as one another: a crystal pose in THIS
     structure, a reference pose that cross_dock_rmsd measured in another one, and nothing at
     all. The anchor has the second, so "not applicable" on its own was wrong for it. */
  if(mine.pose_rmsd!==null && mine.pose_rmsd!==undefined) return num(mine.pose_rmsd,4)+" Å";
  const run=(X.runs||[]).find(r=>r.run===mine.run && r.anchor);
  if(run) return NA+" in this crystal — "+num(run.anchor.rmsd,4)+" ±"+num(run.anchor.bar,3)+
    " Å against its pose in "+run.anchor.reference+
    (run.anchor.decides===false ? ", which does not decide the criterion" : "");
  return NA+" — no crystal holds it in this pocket";
}
function drawDynamics(){
  const host=$("#md"); host.replaceChildren();
  const pockets=(X&&X.md)||[];
  const mine=pockets.find(p=>selected && p.id && (p.pocket||"").toUpperCase().endsWith(selected.id.toUpperCase()));
  const measured=pockets.filter(p=>(p.stages||[]).length);
  $("#mdr").textContent = measured.length ? measured.length+" measured" : "queued";
  if(!pockets.length){ host.append(Object.assign(el("div","ph"),
    {textContent:"no molecular dynamics record"})); return; }
  if(mine){
    host.append(kv([
      ["SYSTEM", mine.pdb+" + "+(mine.compound||mine.pocket)],
      ["STARTING POSE", mine.pose!==null&&mine.pose!==undefined
        ? "mode "+mine.pose+" of "+(mine.run||NR)+", "+(num(mine.pose_score,3)??NR)+" kcal/mol" : NR],
      ["POSE RMSD", poseRmsd(mine)],
      ["BOX / WATERS", mine.box_side? num(mine.box_side,1)+" Å · "+(mine.waters??NR) : NR],
      ["LENGTH AS WRITTEN", mine.ns_total? mine.ns_total+" ns" : NR],
      ["TRAJECTORIES", (mine.stages||[]).length ? String((mine.stages||[]).length) : "none measured yet"],
    ]));
    if(!(mine.stages||[]).length) host.append(Object.assign(el("div","note"),{textContent:
      "Nothing is measured here yet. This panel shows what was SET UP, which is not a result."}));
  } else host.append(Object.assign(el("div","ph"),
    {textContent:"no md package for the selected compound"}));

  measured.forEach(p=>{
    const order={production:0, equilibration:1};
    (p.stages||[]).slice().sort((a,b)=>(order[a.label]??9)-(order[b.label]??9)).forEach(st=>{
      const label=p.pdb+" "+p.pocket+(p.compound?" · "+p.compound:"")+" · "+st.label;
      fold(host, st.label+" trajectory", label, (Object.keys(st.series||{}).length)||0, m=>{
        if(mine!==p) m.append(Object.assign(el("div","note"),{textContent:
          "A DIFFERENT SYSTEM from the compound selected above: "+label+
          ". Its numbers say nothing about the selected compound."}));
        if(st.status) m.append(Object.assign(el("div","note"),{textContent:
          "EQUILIBRATION CONTROL: "+st.status+(st.drift!==null&&st.drift!==undefined
            ? " · "+st.segid+" backbone drifts "+num(st.drift,4)+" Å against a "+
              num(st.threshold)+" Å threshold" : "")}));
        if(st.restraints) m.append(Object.assign(el("div","fa"),
          {style:"font-size:10px",textContent:"restraints: "+st.restraints}));
        const ser=st.series||{};
        if(ser.error) m.append(Object.assign(el("div","ph"),{textContent:ser.error}));
        else {
          const g=el("div","plots");
          const segs=Object.keys(ser).filter(k=>k.startsWith("backbone")&&!k.endsWith("observed"));
          segs.forEach(k=>g.append(plot(k.toUpperCase()+" RMSD (Å)", ser[k], "#ffab3d", false,
            ser[k+" observed"])));
          ["ligand RMSD","COM distance","H-bonds"].forEach(k=>{
            if(ser[k]) g.append(plot(k.toUpperCase()+(k==="H-bonds"?" (count)":" (Å)"),
              ser[k], "#35e39a", false)); });
          m.append(g);
          m.append(Object.assign(el("div","fa"),{style:"font-size:9.5px",textContent:
            (st.frames||"")+" · "+(st.span_ns||"?")+" ns · solid = full segment (the "+
            "pre-registered control), dashed = observed residues only, a second look taken AFTER "+
            "the control failed and not a planned test"}));
        }
        if((st.hbonds||[]).length){
          m.append(Object.assign(el("div","u"),{textContent:"CONTACTS THE LIGAND KEEPS"}));
          m.append(table(["donor","acceptor","frames","share"], st.hbonds.map(h=>[
            h.donor,h.acceptor,{n:h.present+" / "+h.analysed},{n:num(h.share,2)}]),
            {num:["frames","share"]}));
          m.append(Object.assign(el("div","fa"),{style:"font-size:10px",textContent:
            "Share of STORED FRAMES, never a fraction of time: bonds form and break faster than "+
            "frames are written. With no apo simulation this cannot say the contact holds the ligand."}));
        }
        if((st.rmsf||[]).length){
          m.append(Object.assign(el("div","u"),{textContent:"HOW MUCH EACH POCKET RESIDUE MOVED (RMSF)"}));
          const box=el("div","scroll");
          box.append(table(["residue","rmsf Å"], st.rmsf.map(x=>[x.residue,{n:num(x.rmsf,4)}]),
            {num:["rmsf Å"]}));
          m.append(box);
        }
        if((st.offered||[]).length){
          m.append(Object.assign(el("div","u"),{textContent:"WHAT THIS CAN AND CANNOT BE READ AS"}));
          st.offered.forEach(o=>m.append(Object.assign(el("div","fa"),
            {style:"font-size:10px;margin:3px 0",textContent:"· "+o})));
        }
      });
    });
  });
}
function drawLower(){ drawContext(); drawExperiment(); drawDynamics(); drawPocket(); }
drawSet(); drawFeed(); drawLower();
</script>
</body>
</html>
"""
