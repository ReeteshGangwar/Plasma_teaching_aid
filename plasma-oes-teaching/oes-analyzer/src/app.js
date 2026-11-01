// ===================== UI =====================
(() => {
const $=id=>document.getElementById(id);
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const f=(v,d=2)=>Number.isFinite(v)?v.toFixed(d):'—';
const fe=(v,d=1)=>Number.isFinite(v)&&v>0? v.toExponential(d).replace('e+','×10^').replace(/\^(\d+)/,(m,p)=>'<sup>'+p+'</sup>') : '—';
const esc=s=>String(s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
let uid=1;

// ---------- state ----------
const ST = {
  conds:[], active:0, tab:'data',
  units:'energy', shiftMode:'auto', shiftManual:0, sat:65000,
  peak:{thr:8, sep:0.3, tol:0.25}, hw:0.35,
  lines: defaultLines(),
  rot:{lo:334.6, hi:337.6, fixW:false, w:0.15},
  hr:{wI:0.025, tgMode:'rot', Tg:400, P:1, wS:STARK696.wS, C:null, Cunc:15, freeG:false},
  neLTE:null, imode:'area', manual:{}, hideMissing:false
};
function defaultLines(){ return AR_BOLTZ_DEFAULT.map(d=>({...d, include: d.src!=='verify' && ![811.531,763.511,419.832,420.068].includes(d.l)})); }
const R2_696 = r2coul(AR696.Eion,AR696.Eup,AR696.lup)-r2coul(AR696.Eion,AR696.Elow,AR696.llow);
ST.hr.C = +vdwC(AR696.l,R2_696).toFixed(3);
const newCond = name => ({id:uid++, name, plasma:[], bg:[], expP:null, expB:null, hr:{plasma:[], bg:[]}, example:false, res:{}, ver:0});
const cond = () => ST.conds[ST.active];

// ---------- toasts & saving ----------
let tt; function toast(m){ const t=$('toast'); t.textContent=m; t.classList.add('show'); clearTimeout(tt); tt=setTimeout(()=>t.classList.remove('show'),2600); }
let dlP=null; function dl(){ if(!dlP) dlP = (window.claude && window.claude.use) ? window.claude.use('downloads').catch(()=>null) : Promise.resolve(null); return dlP; }
async function saveFile(name, data){
  const d = await dl();
  if(d){ try{ await d.save({filename:name, data}); toast('Saved '+name); } catch(e){ if(e && e.code==='declined') toast('Save cancelled'); else toast('Could not save: '+(e&&e.message||e)); } return; }
  try{ const blob = data instanceof Blob ? data : new Blob([data],{type:'text/plain'}); const a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download=name; document.body.appendChild(a); a.click(); setTimeout(()=>{ URL.revokeObjectURL(a.href); a.remove(); },500); toast('Saved '+name); }
  catch(e){ toast('Saving is not available here'); }
}
const csv = rows => rows.map(r=>r.map(v=> typeof v==='number' ? (Number.isFinite(v)? String(+v.toPrecision(8)) : '') : '"'+String(v??'').replace(/"/g,'""')+'"').join(',')).join('\n');
const safe = s => s.replace(/[^\w.-]+/g,'_').replace(/_+/g,'_').slice(0,60);
async function savePlot(id,name){ const gd=$(id); if(!gd || !gd.data) return; try{ const url=await Plotly.toImage(gd,{format:'png',scale:2,width:gd.clientWidth,height:gd.clientHeight}); const b=await (await fetch(url)).blob(); saveFile(name,b); }catch(e){ toast('Could not render PNG'); } }

// ---------- processing cache ----------
const cache=new Map();
function proc(c){ if(!c||!c.plasma.length) return null; const key=c.id+':'+c.ver+':'+ST.sat+':'+ST.shiftMode+':'+ST.shiftManual;
  if(cache.has(key)) return cache.get(key);
  const P=processCondition(c,{sat:ST.sat}); if(!P) return null;
  let sh={shift:0,n:0}; if(ST.shiftMode==='auto') sh=autoShift(P.grid,P.mean); else if(ST.shiftMode==='manual') sh={shift:ST.shiftManual,n:0};
  P.shift=sh; P.x=Float64Array.from(P.grid,v=>v-sh.shift);
  cache.set(key,P); return P; }
function procHr(c){ if(!c||!c.hr.plasma.length) return null; const key='hr'+c.id+':'+c.ver; if(cache.has(key)) return cache.get(key); const P=processCondition({plasma:c.hr.plasma,bg:c.hr.bg,expP:c.expP,expB:c.expB},{sat:ST.sat}); cache.set(key,P); return P; }
const touch=c=>{ c.ver++; c.res={}; };

// ---------- plotting helpers ----------
function baseLayout(extra={}){ const ink=css('--ink'), mut=css('--muted'), line=css('--line');
  return Object.assign({ paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)', font:{family:css('--f-body'),size:12,color:ink},
    margin:{l:62,r:16,t:16,b:46}, hovermode:'closest', hoverlabel:{font:{family:css('--f-mono')}},
    xaxis:{gridcolor:line, zerolinecolor:line, linecolor:mut, tickcolor:mut, ticks:'outside'},
    yaxis:{gridcolor:line, zerolinecolor:line, linecolor:mut, tickcolor:mut, ticks:'outside'},
    legend:{orientation:'h', y:1.02, yanchor:'bottom', x:0, font:{size:11.5}} }, extra); }
const PCONF={responsive:true, displaylogo:false, modeBarButtonsToRemove:['toImage','select2d','lasso2d','autoScale2d']};
const hasPlotly = () => typeof Plotly!=='undefined';
function draw(id,data,layout){ if(!hasPlotly()){ $(id).innerHTML='<p class="sub" style="padding:20px">The plotting library did not load. Use the downloadable offline version of this file.</p>'; return; } Plotly.react(id,data,layout,PCONF); }
const alpha=(hex,a)=>{ const m=hex.match(/^#?([\da-f]{2})([\da-f]{2})([\da-f]{2})$/i); return m? `rgba(${parseInt(m[1],16)},${parseInt(m[2],16)},${parseInt(m[3],16)},${a})` : hex; };
const SPC = { // species colours (categorical, both themes)
  'Ar I':'#3f7fd6','Ar II':'#7a62c9','H':'#d4453b','O I':'#2f9e6b','OH':'#c77a12','N₂ SPS':'#0e8fa3','N₂⁺':'#b04a91','Na I':'#8a8a3a' };
const spColor = sp => { for(const k of Object.keys(SPC)) if(sp.startsWith(k)) return SPC[k]; if(sp.startsWith('H I')) return SPC.H; if(sp.startsWith('N₂ ')) return SPC['N₂ SPS']; return '#888'; };

// ===================== DATA TAB =====================
function sortFile(name){ const n=name.toLowerCase(); const hr=/(696|hires|hi-res|high.?res|_hr\b|\bhr_|-hr)/.test(n); const bg=/(bg|background|dark|blank)/.test(n); return {hr, bg}; }
async function readFiles(fileList){ const out=[]; for(const file of fileList){ try{ const t=await file.text(); const s=parseSpectrum(t); out.push({name:file.name, ...s}); } catch(e){ out.push({name:file.name, error:e.message}); } } return out; }
function addParsed(c, parsed, role){ const errs=[];
  for(const p of parsed){ if(p.error){ errs.push(`${p.name}: ${p.error}`); continue; }
    let r=role; if(r==='auto'){ const s=sortFile(p.name); r = s.hr ? (s.bg?'hrbg':'hrp') : (s.bg?'bg':'plasma'); }
    ({plasma:c.plasma, bg:c.bg, hrp:c.hr.plasma, hrbg:c.hr.bg})[r].push({name:p.name, x:p.x, y:p.y}); }
  touch(c); if(errs.length) toast('Skipped: '+errs.join('; ')); }
function bindDrop(el,onFiles){ el.addEventListener('dragover',e=>{ e.preventDefault(); el.classList.add('over'); }); el.addEventListener('dragleave',()=>el.classList.remove('over'));
  el.addEventListener('drop',async e=>{ e.preventDefault(); el.classList.remove('over'); if(e.dataTransfer.files.length) onFiles(await readFiles(e.dataTransfer.files)); }); }
bindDrop($('quickDrop'), parsed=>{ const c=newCond('Condition '+(ST.conds.length+1)); addParsed(c,parsed,'auto'); ST.conds.push(c); ST.active=ST.conds.length-1; refreshAll(); });
$('quickInput').addEventListener('change',async e=>{ const parsed=await readFiles(e.target.files); e.target.value=''; const c=newCond('Condition '+(ST.conds.length+1)); addParsed(c,parsed,'auto'); ST.conds.push(c); ST.active=ST.conds.length-1; refreshAll(); });
$('addCond').onclick=()=>{ ST.conds.push(newCond('Condition '+(ST.conds.length+1))); ST.active=ST.conds.length-1; refreshAll(); };
$('loadEx').onclick=()=>{ for(const v of [0,1]){ const e=makeExample(v); const c=newCond(e.name); Object.assign(c,{plasma:e.plasma,bg:e.bg,example:true}); c.hr.plasma=e.hr.plasma; c.hr.bg=e.hr.bg; ST.conds.push(c); } ST.active=ST.conds.length-2; refreshAll(); toast('Two synthetic example conditions added'); };
$('clearAll').onclick=()=>{ if(!ST.conds.length) return; if($('clearAll').dataset.arm!=='1'){ $('clearAll').dataset.arm='1'; $('clearAll').textContent='Click again to remove all'; setTimeout(()=>{ $('clearAll').dataset.arm=''; $('clearAll').textContent='Remove all'; },3000); return; } ST.conds=[]; ST.active=0; cache.clear(); $('clearAll').dataset.arm=''; $('clearAll').textContent='Remove all'; refreshAll(); };

function fileList(c, key, label, hint){
  const arr = key.startsWith('hr') ? (key==='hrp'? c.hr.plasma : c.hr.bg) : c[key];
  const moves = {plasma:[['bg','→ background']], bg:[['plasma','→ plasma']], hrp:[['hrbg','→ hi-res bg']], hrbg:[['hrp','→ hi-res plasma']]}[key];
  return `<div class="drop" data-c="${c.id}" data-k="${key}"><div class="hd"><h3>${label} <span class="sub">(${arr.length})</span></h3><label class="btn" style="display:inline-block">Add<input type="file" multiple hidden data-c="${c.id}" data-k="${key}"></label></div>
    <div class="chips">${arr.map((fl,i)=>`<div class="chip"><span class="nm" title="${esc(fl.name)}">${esc(fl.name)}<br><span class="mt">${fl.x.length} pts · ${fl.x[0].toFixed(1)}–${fl.x[fl.x.length-1].toFixed(1)} nm</span></span><span class="acts">${moves.map(([to,t])=>`<button class="link" data-mv="${to}" data-c="${c.id}" data-k="${key}" data-i="${i}">${t}</button>`).join('')}<button class="link" data-rm="1" data-c="${c.id}" data-k="${key}" data-i="${i}" aria-label="Remove ${esc(fl.name)}">✕</button></span></div>`).join('') || `<p class="sub">${hint}</p>`}</div></div>`;
}
function renderData(){
  const L=$('condList');
  if(!ST.conds.length){ L.innerHTML='<div class="msg info">No measurements yet. Drop files in the box above, add a condition, or load the synthetic example data to explore the tool.</div>'; return; }
  L.innerHTML = ST.conds.map((c,ci)=>{ const P=proc(c);
    const st=[];
    if(c.plasma.length && c.plasma.length<3) st.push(`<div class="msg warn">${c.plasma.length} plasma spectrum(s): three replicates are needed for a standard deviation.</div>`);
    if(c.plasma.length && !c.bg.length) st.push('<div class="msg warn">No background file: the net spectrum includes dark counts and ambient light.</div>');
    if(P){ st.push(`<div class="msg good">Net spectrum: ${P.nRep} replicate(s) − mean of ${c.bg.length} background(s), ${P.grid.length} points, ${P.x[0].toFixed(1)}–${P.x[P.x.length-1].toFixed(1)} nm, step ${f(P.step,3)} nm, noise σ ≈ ${f(P.noise,1)}. Wavelength shift applied: ${f(P.shift.shift,3)} nm${ST.shiftMode==='auto'?` (from ${P.shift.n} Ar I lines)`:''}.</div>`);
      if(P.satFiles.length) st.push(`<div class="msg bad">Saturated pixels (≥ ${ST.sat}): ${P.satFiles.join('; ')}. Line intensities near saturation are unreliable.</div>`);
      const grids=new Set(c.plasma.concat(c.bg).map(fl=>fl.x.length+':'+fl.x[0].toFixed(3))); if(grids.size>1) st.push('<div class="msg info">The files use different wavelength grids; all were interpolated onto the first plasma file.</div>'); }
    return `<div class="cond ${ci===ST.active?'active':''}" data-ci="${ci}">
      <div class="row" style="justify-content:space-between"><div class="row"><input type="text" class="wide" style="min-width:240px;font-weight:600" value="${esc(c.name)}" data-name="${ci}" aria-label="Condition name">${c.example?'<span class="badge">synthetic</span>':''}</div>
        <div class="row"><button class="btn" data-act="${ci}">${ci===ST.active?'Active':'Make active'}</button><button class="btn" data-del="${ci}">Delete</button></div></div>
      <div class="files">${fileList(c,'plasma','Plasma spectra','Drop the 3 plasma files')}${fileList(c,'bg','Background spectra','Drop the 3 background files')}${fileList(c,'hrp','Hi-res Ar 696.5 nm','Optional, for nₑ')}${fileList(c,'hrbg','Hi-res background','Optional')}</div>
      <div class="row"><label class="f">Plasma exposure (ms)<input type="number" value="${c.expP??''}" data-exp="P" data-ci="${ci}" placeholder="optional"></label><label class="f">Background exposure (ms)<input type="number" value="${c.expB??''}" data-exp="B" data-ci="${ci}" placeholder="optional"></label><span class="sub">If both are given, the background is scaled by their ratio before subtraction.</span></div>
      <div class="status">${st.join('')}</div></div>`; }).join('');
  L.querySelectorAll('.drop[data-c]').forEach(el=>{ const c=ST.conds.find(x=>x.id==el.dataset.c); bindDrop(el,parsed=>{ addParsed(c,parsed,el.dataset.k); refreshAll(); }); });
  L.querySelectorAll('input[type=file][data-c]').forEach(inp=>inp.addEventListener('change',async e=>{ const c=ST.conds.find(x=>x.id==inp.dataset.c); addParsed(c,await readFiles(e.target.files),inp.dataset.k); refreshAll(); }));
  const arrOf=(c,k)=>({plasma:c.plasma,bg:c.bg,hrp:c.hr.plasma,hrbg:c.hr.bg})[k];
  L.querySelectorAll('[data-mv]').forEach(b=>b.onclick=()=>{ const c=ST.conds.find(x=>x.id==b.dataset.c); const fl=arrOf(c,b.dataset.k).splice(+b.dataset.i,1)[0]; arrOf(c,b.dataset.mv).push(fl); touch(c); refreshAll(); });
  L.querySelectorAll('[data-rm]').forEach(b=>b.onclick=()=>{ const c=ST.conds.find(x=>x.id==b.dataset.c); arrOf(c,b.dataset.k).splice(+b.dataset.i,1); touch(c); refreshAll(); });
  L.querySelectorAll('[data-name]').forEach(i=>i.onchange=()=>{ ST.conds[+i.dataset.name].name=i.value; renderCondSel(); });
  L.querySelectorAll('[data-act]').forEach(b=>b.onclick=()=>{ ST.active=+b.dataset.act; refreshAll(); });
  L.querySelectorAll('[data-del]').forEach(b=>b.onclick=()=>{ ST.conds.splice(+b.dataset.del,1); ST.active=Math.max(0,Math.min(ST.active,ST.conds.length-1)); refreshAll(); });
  L.querySelectorAll('[data-exp]').forEach(i=>i.onchange=()=>{ const c=ST.conds[+i.dataset.ci]; const v=parseFloat(i.value); c[i.dataset.exp==='P'?'expP':'expB']= v>0? v : null; touch(c); refreshAll(); });
}
function renderCondSel(){ const s=$('activeCond'); s.innerHTML = ST.conds.length ? ST.conds.map((c,i)=>`<option value="${i}" ${i===ST.active?'selected':''}>${esc(c.name)}</option>`).join('') : '<option>— none —</option>'; $('exBadge').hidden=!(cond()&&cond().example); }
$('activeCond').onchange=e=>{ ST.active=+e.target.value; refreshAll(); };
function segBind(id,cb){ $(id).addEventListener('click',e=>{ const b=e.target.closest('button'); if(!b) return; for(const x of $(id).children) x.setAttribute('aria-pressed',String(x===b)); cb(b.dataset.v); }); }
const segSet=(id,v)=>{ for(const x of $(id).children) x.setAttribute('aria-pressed',String(x.dataset.v===v)); };
segBind('unitSeg',v=>{ ST.units=v; refreshAll(); });
segBind('shiftSeg',v=>{ ST.shiftMode=v; refreshAll(); });
$('shiftVal').onchange=e=>{ ST.shiftManual=parseFloat(e.target.value)||0; if(ST.shiftMode==='manual') refreshAll(); };
$('satLvl').onchange=e=>{ ST.sat=parseFloat(e.target.value)||Infinity; refreshAll(); };

// ===================== SPECTRUM TAB =====================
const ZOOMS=[['Full',null],['OH 306–313',[304,314]],['N₂ 313–340',[312,340]],['N₂ 350–382',[349,382]],['N₂⁺ 391',[386,395]],['Ar 5p 414–432',[414,432]],['Hβ',[484,488]],['Hα',[654,659]],['Ar 4p 690–812',[690,812]],['O 777',[774,781]]];
$('zooms').innerHTML = ZOOMS.map((z,i)=>`<button class="btn" data-z="${i}">${z[0]}</button>`).join('');
$('zooms').addEventListener('click',e=>{ const b=e.target.closest('[data-z]'); if(!b||!hasPlotly()) return; const z=ZOOMS[+b.dataset.z][1]; Plotly.relayout('plotSpec', z? {'xaxis.range':z,'yaxis.autorange':true} : {'xaxis.autorange':true,'yaxis.autorange':true}); setTimeout(autoYSpec,30); });
function autoYSpec(){ const gd=$('plotSpec'), P=proc(cond()); if(!gd.layout||!P) return; const r=gd.layout.xaxis.range; if(!r) return; let mx=-Infinity, mn=Infinity; for(let i=0;i<P.x.length;i++){ if(P.x[i]<r[0]||P.x[i]>r[1]) continue; mx=Math.max(mx,P.mean[i]+P.sd[i]); mn=Math.min(mn,P.mean[i]-P.sd[i]); } if(Number.isFinite(mx) && !$('logY').checked) Plotly.relayout('plotSpec',{'yaxis.range':[Math.min(0,mn),mx*1.12]}); }
['showReps','showBg','logY','showLabels'].forEach(id=>$(id).onchange=()=>renderSpec());
['pkThr','pkSep','pkTol'].forEach(id=>$(id).onchange=()=>{ ST.peak={thr:+$('pkThr').value||8, sep:+$('pkSep').value||0.3, tol:+$('pkTol').value||0.25}; renderSpec(); });
let lastPeaks=[];
function renderSpec(){
  const c=cond(), P=proc(c);
  if(!P){ draw('plotSpec',[],baseLayout()); $('specStatus').innerHTML='<div class="msg info">Load plasma spectra in the Data tab first.</div>'; $('peakTbl').innerHTML=''; return; }
  const acc=css('--accent'), mut=css('--muted');
  const x=Array.from(P.x), up=Array.from(P.mean,(v,i)=>v+P.sd[i]), lo=Array.from(P.mean,(v,i)=>v-P.sd[i]);
  const tr=[{x,y:lo,mode:'lines',line:{width:0},hoverinfo:'skip',showlegend:false},
    {x,y:up,mode:'lines',line:{width:0},fill:'tonexty',fillcolor:alpha(acc,0.22),name:'± 1 SD (n = '+P.nRep+')',hoverinfo:'skip'},
    {x,y:Array.from(P.mean),mode:'lines',line:{color:acc,width:1.3},name:'Mean net spectrum',hovertemplate:'%{x:.3f} nm<br>%{y:.1f}<extra></extra>'}];
  if($('showReps').checked) P.nets.forEach((r,i)=>tr.push({x,y:Array.from(r),mode:'lines',line:{width:0.8,color:alpha(mut,0.75)},name:'Replicate '+(i+1),hovertemplate:'%{x:.3f} nm<br>%{y:.1f}<extra>rep '+(i+1)+'</extra>'}));
  if($('showBg').checked) tr.push({x,y:Array.from(P.bgMean),mode:'lines',line:{width:1,color:'#c77a12',dash:'dot'},name:'Mean background'});
  const pk=findPeaks(P.x,P.mean,P.noise,ST.peak); lastPeaks=pk;
  const idp=pk.filter(p=>p.id);
  const bySp={}; for(const p of idp){ const k=p.id.sp.split(' (')[0].replace(/ \(.*$/,''); (bySp[k] ||= []).push(p); }
  for(const [sp,arr] of Object.entries(bySp)) tr.push({x:arr.map(p=>p.l),y:arr.map(p=>p.h),mode:'markers',marker:{size:7,color:spColor(sp),symbol:'triangle-down'},name:sp,
    text:arr.map(p=>`${p.id.sp} ${p.id.l.toFixed(2)} nm${p.id.note?' ('+p.id.note+')':''}`),hovertemplate:'%{text}<br>measured %{x:.3f} nm<br>%{y:.0f}<extra></extra>'});
  const un=pk.filter(p=>!p.id); if(un.length) tr.push({x:un.map(p=>p.l),y:un.map(p=>p.h),mode:'markers',marker:{size:6,color:mut,symbol:'circle-open'},name:'Unassigned',hovertemplate:'%{x:.3f} nm<br>%{y:.0f}<extra>unassigned</extra>'});
  const ann=[];
  if($('showLabels').checked){ const seen=new Set(), short=sp=>sp.startsWith('H I (')? sp.slice(5,-1) : sp;
    const top=idp.slice().sort((a,b)=>b.h-a.h).filter(p=>{ const k=p.id.sp+p.id.l; if(seen.has(k)) return false; seen.add(k); return true; }).slice(0,24);
    top.forEach((p,j)=>ann.push({x:p.l,y:$('logY').checked?Math.log10(Math.max(p.h,1)):p.h,text:short(p.id.sp)+' '+p.id.l.toFixed(1),showarrow:true,arrowhead:0,ax:0,ay:-18-(j%3)*12,font:{size:10,color:spColor(p.id.sp)},arrowcolor:alpha(mut,0.6)})); }
  const gd=$('plotSpec'), keepX = gd.layout && gd.layout.xaxis && gd.layout.xaxis.range && gd.dataset.cid==c.id ? gd.layout.xaxis.range : null;
  draw('plotSpec',tr,baseLayout({xaxis:{...baseLayout().xaxis,title:{text:'Wavelength (nm, corrected)'},range:keepX||undefined,autorange:!keepX},yaxis:{...baseLayout().yaxis,title:{text:ST.units==='energy'?'Intensity (response-corrected)':'Counts'},type:$('logY').checked?'log':'linear'},annotations:ann}));
  gd.dataset.cid=c.id;
  const st=[];
  if(P.satFiles.length) st.push(`<div class="msg bad">Saturated pixels: ${P.satFiles.join('; ')}</div>`);
  st.push(`<div class="msg info">${pk.length} peaks above ${ST.peak.thr}σ (σ ≈ ${f(P.noise,1)}); ${idp.length} assigned. Band features (OH, N₂, N₂⁺) are matched at their band heads, so assignments within dense bands are approximate. Click a row to zoom.</div>`);
  $('specStatus').innerHTML=st.join('');
  $('peakTbl').innerHTML = `<thead><tr><th class="num">λ measured (nm)</th><th class="num">Height</th><th class="num">SNR</th><th>Assignment</th><th class="num">λ ref (nm)</th><th class="num">Δλ (nm)</th><th>Note</th></tr></thead><tbody>${pk.map((p,i)=>`<tr class="click" data-pk="${i}"><td class="num">${p.l.toFixed(3)}</td><td class="num">${p.h.toFixed(0)}</td><td class="num">${p.snr.toFixed(0)}</td><td>${p.id?`<span style="color:${spColor(p.id.sp)};font-weight:600">${esc(p.id.sp)}</span>`:'<span class="sub">unassigned</span>'}</td><td class="num">${p.id?p.id.l.toFixed(3):''}</td><td class="num">${p.id?(p.l-p.id.l).toFixed(3):''}</td><td>${p.id&&p.id.note?esc(p.id.note):''}</td></tr>`).join('')}</tbody>`;
  $('peakTbl').querySelectorAll('[data-pk]').forEach(r=>r.onclick=()=>{ const p=pk[+r.dataset.pk]; Plotly.relayout('plotSpec',{'xaxis.range':[p.l-3,p.l+3]}); setTimeout(autoYSpec,30); $('plotSpec').scrollIntoView({behavior:'smooth',block:'center'}); });
}
$('specPng').onclick=()=>savePlot('plotSpec',safe(cond()?.name||'spectrum')+'_spectrum.png');
$('specCsv').onclick=()=>{ const c=cond(), P=proc(c); if(!P) return; const rows=[['wavelength_nm_corrected','wavelength_nm_raw','mean_net','sd_net',...P.nets.map((_,i)=>'net_rep'+(i+1)),'background_mean']]; for(let i=0;i<P.x.length;i++) rows.push([P.x[i],P.grid[i],P.mean[i],P.sd[i],...P.nets.map(r=>r[i]),P.bgMean[i]]); saveFile(safe(c.name)+'_net_spectrum.csv',csv(rows)); };
$('peakCsv').onclick=()=>{ const c=cond(); if(!lastPeaks.length) return; const rows=[['lambda_measured_nm','height','snr','species','lambda_ref_nm','delta_nm','note']]; for(const p of lastPeaks) rows.push([p.l,p.h,p.snr,p.id?p.id.sp:'',p.id?p.id.l:NaN,p.id?p.l-p.id.l:NaN,p.id&&p.id.note||'']); saveFile(safe(c.name)+'_peaks.csv',csv(rows)); };

// ===================== BOLTZMANN TAB =====================
function boltzFor(c){ const P=proc(c); if(!P && ST.imode!=='manual') return null;
  const opts={hw:ST.hw,energy:ST.units==='energy',mode:ST.imode,manual:ST.manual,tol:ST.peak.tol};
  if(ST.imode==='peakH'||ST.imode==='peakA') opts.peaks=findPeaks(P.x,P.mean,P.noise,ST.peak);
  return P ? boltzmann(P.x,P.nets,P.noise,ST.lines,opts) : boltzmann(Float64Array.from([0,1]),[],1,ST.lines,opts); }
function neForLTE(c){ if(ST.neLTE>0) return ST.neLTE; const h=c && c.res.hr; return h && Number.isFinite(h.ne) && h.ne>0 ? h.ne : null; }
function renderBoltz(){
  const c=cond(); const B=c && boltzFor(c);
  $('unitNote').textContent = ST.units==='energy' ? 'Using response-corrected (energy) intensities: the λ factor is included.' : 'Using photon counts: the λ factor is omitted (I ∝ photon rate).';
  if(!B){ draw('plotBoltz',[],baseLayout()); $('boltzKpis').innerHTML=''; $('lteBox').innerHTML='<div class="msg info">Load plasma spectra in the Data tab first, or switch the intensity source to Manual entry.</div>'; $('boltzTbl').innerHTML=''; return; }
  c.res.boltz=B;
  const acc=css('--accent'), mut=css('--muted');
  const avail=B.rows.filter(r=>!r.missing && Number.isFinite(r.y));
  const inc=avail.filter(r=>r.include), exc=avail.filter(r=>!r.include);
  const pt=(rows,name,color,sym)=>({x:rows.map(r=>r.E),y:rows.map(r=>r.y),mode:'markers',name,marker:{size:9,color,symbol:sym,line:{width:1.2,color}},
    error_y:{type:'data',array:rows.map(r=>Number.isFinite(r.sy)?r.sy:0),visible:true,color:alpha(color,0.7),thickness:1.2,width:3},
    text:rows.map(r=>`${r.l.toFixed(2)} nm · ${r.up}→${r.low}${r.meta?' (metastable lower)':''}`),hovertemplate:'%{text}<br>E<sub>k</sub> = %{x:.3f} eV<br>y = %{y:.3f}<extra></extra>'});
  const tr=[];
  tr.push(pt(inc.filter(r=>!r.meta),'Included · resonance/other lower level',acc,'circle'));
  tr.push(pt(inc.filter(r=>r.meta),'Included · metastable lower level (1s₅, 1s₃)',acc,'diamond'));
  if(exc.length) tr.push({...pt(exc,'Excluded',mut,'circle-open'),error_y:{visible:false}});
  if(B.fit){ const Es=avail.map(r=>r.E), x0=Math.min(...Es)-0.1, x1=Math.max(...Es)+0.1; tr.push({x:[x0,x1],y:[B.fit.a+B.fit.b*x0,B.fit.a+B.fit.b*x1],mode:'lines',line:{color:acc,width:1.5,dash:'dash'},name:`Fit: T = ${f(B.T,0)} K`,hoverinfo:'skip'}); }
  const ann=avail.map(r=>({x:r.E,y:r.y,text:r.l.toFixed(1),showarrow:false,xshift:0,yshift:12,font:{size:9.5,color:r.include?css('--ink'):mut}}));
  draw('plotBoltz',tr,baseLayout({xaxis:{...baseLayout().xaxis,title:{text:'Upper-level energy E<sub>k</sub> (eV)'}},yaxis:{...baseLayout().yaxis,title:{text:ST.units==='energy'?'ln(Iλ / gA)':'ln(I / gA)'}},annotations:ann,legend:{orientation:'h',y:-0.22,x:0,font:{size:11}},margin:{l:62,r:16,t:16,b:90}}));
  const Tr=B.Treps||[];
  $('boltzKpis').innerHTML = B.fit ? [
    ['T<sub>exc</sub> (fit to mean)', `${f(B.T,0)} <small>± ${f(B.sT,0)} K</small>`,1],
    ['in eV', `${f(B.T*KB_EV,3)} <small>± ${f(B.sT*KB_EV,3)}</small>`],
    ['Replicates (mean ± SD)', Tr.length? `${f(mean(Tr),0)} <small>± ${f(sdev(Tr),0)} K (n=${Tr.length})</small>`:'—'],
    ['R²', f(B.fit.r2,3)], ['Lines used', `${B.use.length} <small>of ${avail.length} found</small>`], ['E<sub>k</sub> span', `${f(B.Espan,2)} <small>eV</small>`]
  ].map(([k,v,m])=>`<div class="kpi ${m?'main':''}"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('') : '<div class="kpi"><span class="k">Fit</span><span class="v">needs ≥ 2 lines</span></div>';
  // LTE box
  const L=lteChecks(B,{ne:neForLTE(c)}); c.res.lte=L;
  const tbl = `<table style="margin-top:4px"><tbody>
    <tr><td>T from 4p levels only</td><td class="num">${L.T4? (Number.isFinite(L.T4.T)? f(L.T4.T,0)+' K':'no slope'):'—'}</td></tr>
    <tr><td>T from higher levels only (E ≥ 14 eV)</td><td class="num">${L.T5? (Number.isFinite(L.T5.T)? f(L.T5.T,0)+' K':'no slope'):'—'}</td></tr>
    <tr><td>T from the 4p → higher-level step</td><td class="num">${Number.isFinite(L.T45)? f(L.T45,0)+' K':'—'}</td></tr>
    <tr><td>McWhirter limit, full LTE</td><td class="num">${fe(L.mcwFull)} cm⁻³</td></tr>
    <tr><td>McWhirter limit, partial LTE (4s–4p gap)</td><td class="num">${fe(L.mcwPartial)} cm⁻³</td></tr>
    <tr><td>nₑ used for the test</td><td class="num">${L.ne? fe(L.ne)+' cm⁻³':'not set'}</td></tr></tbody></table>`;
  const same = (L.sameUpper||[]).length ? `<p class="sub">Same-upper-level pairs (should agree): ${L.sameUpper.map(s=>`${s.up}: ${s.lines.map(l=>l.toFixed(1)).join('/')} Δ = ${s.dy.toFixed(2)}`).join(' · ')}</p>` : '';
  $('lteBox').innerHTML = L.verdict==='insufficient' ? '<div class="msg warn">Not enough lines for an LTE check.</div>' :
    `<div class="verdict ${L.verdict==='consistent'?'ok':'no'}"><h3>${L.verdict==='consistent'?'Consistent with a Boltzmann distribution':'Boltzmann (LTE) assumption not supported'}</h3>
     ${L.reasons.length?`<ul>${L.reasons.map(r=>`<li>${r}</li>`).join('')}</ul>`:''}
     ${L.notes.map(n=>`<p class="sub">${n}</p>`).join('')}
     ${L.verdict!=='consistent'?'<p><b>Next step:</b> treat T<sub>exc</sub> as an indicator only and estimate Tₑ with a collisional–radiative model (see the note below).</p>':'<p class="sub">Even when consistent, T<sub>exc</sub> equals Tₑ only if nₑ meets the McWhirter criterion.</p>'}
     ${tbl}${same}</div>`;
  // intensity-source note
  const nMiss=B.rows.filter(r=>r.missing).length;
  $('imodeNote').innerHTML = ({area:'Each line is integrated over ± half-width around its tabulated wavelength (after the wavelength correction), with a linear baseline, for every replicate.',
    peakA:`Only lines detected by the peak finder in tab 2 (threshold ${ST.peak.thr}σ, tolerance ${ST.peak.tol} nm) are used; each is integrated around its detected centre.`,
    peakH:`Only lines detected by the peak finder in tab 2 are used; the baseline-corrected peak height stands in for intensity (valid only if all lines have the same instrumental width).`,
    manual:'Type intensities in the "Manual I" column or paste a list below. These must already be response-corrected and background-subtracted, in the same units for every line.'})[ST.imode] + (nMiss? ` ${nMiss} line(s) unavailable.`:'');
  $('manualTools').style.display = ST.imode==='manual' ? '' : 'none';
  // table
  const shown=B.rows.map((r,i)=>[r,i]).filter(([r])=>!(ST.hideMissing && r.missing));
  $('boltzTbl').innerHTML = `<thead><tr><th>Use</th><th class="num">λ (nm)</th><th>Transition</th><th class="num">E<sub>k</sub> (eV)</th><th class="num">g<sub>k</sub></th><th class="num">A<sub>ki</sub> (s⁻¹)</th><th class="num">I used</th><th class="num">SD</th><th class="num">SNR</th><th class="num">y</th><th class="num">Manual I</th><th class="num">Manual SD</th><th>Flags</th></tr></thead><tbody>${shown.map(([r,i])=>{
    const flags=[]; if(r.src==='verify') flags.push('<span class="flag l" title="constants from memory of NIST; check before quantitative use">verify A,g,E</span>'); if(r.src==='Ohno 2006') flags.push('<span class="flag l">Ohno 2006</span>');
    if(r.meta) flags.push('<span class="flag m">metastable lower</span>'); if(r.blend&&r.blend.length) flags.push(`<span class="flag b" title="${esc(r.blend.join(', '))}">blend risk</span>`); if(r.l===427.217) flags.push('<span class="flag b">near N₂⁺ 427.8</span>'); if(!r.missing && Number.isFinite(r.snr) && r.snr<5) flags.push('<span class="flag l">low SNR</span>'); if(r.missing) flags.push(`<span class="flag l">${r.why||'unavailable'}</span>`);
    const mv=ST.manual[r.l]||{};
    return `<tr class="${r.include&&!r.missing?'':'off'}"><td><input type="checkbox" data-inc="${i}" ${r.include?'checked':''} aria-label="Include ${r.l}"></td><td class="num">${r.l.toFixed(3)}</td><td>${r.up} → ${r.low}</td>
      <td class="num"><input type="number" step="0.001" value="${r.E}" data-ed="E" data-i="${i}"></td><td class="num"><input type="number" step="1" value="${r.g}" data-ed="g" data-i="${i}" style="width:56px"></td><td class="num"><input type="number" step="any" value="${r.A.toExponential(3)}" data-ed="A" data-i="${i}" style="width:100px"></td>
      <td class="num">${r.missing?'':f(r.I,1)}</td><td class="num">${r.missing?'':f(r.sI,1)}</td><td class="num">${r.missing||!Number.isFinite(r.snr)?'':f(r.snr,0)}</td><td class="num">${r.missing?'':f(r.y,3)}</td>
      <td class="num"><input type="number" step="any" value="${mv.I??''}" data-man="I" data-l="${r.l}" style="width:90px" aria-label="Manual intensity ${r.l}"></td><td class="num"><input type="number" step="any" value="${mv.sI??''}" data-man="sI" data-l="${r.l}" style="width:76px" aria-label="Manual SD ${r.l}"></td><td>${flags.join('')}</td></tr>`; }).join('')}</tbody>`;
  $('boltzTbl').querySelectorAll('[data-inc]').forEach(cb=>cb.onchange=()=>{ ST.lines[+cb.dataset.inc].include=cb.checked; renderBoltz(); });
  $('boltzTbl').querySelectorAll('[data-ed]').forEach(inp=>inp.onchange=()=>{ const v=parseFloat(inp.value); if(Number.isFinite(v)&&v>0){ ST.lines[+inp.dataset.i][inp.dataset.ed]=v; renderBoltz(); } });
  $('boltzTbl').querySelectorAll('[data-man]').forEach(inp=>inp.onchange=()=>{ const l=+inp.dataset.l, v=parseFloat(inp.value); const m=ST.manual[l]||{}; if(Number.isFinite(v)&&v>0) m[inp.dataset.man]=v; else delete m[inp.dataset.man]; if(m.I>0||m.sI>0) ST.manual[l]=m; else delete ST.manual[l];
    if(ST.imode!=='manual' && inp.dataset.man==='I' && m.I>0){ ST.imode='manual'; segSet('imodeSeg','manual'); toast('Switched to manual intensities'); } renderBoltz(); });
}
segBind('imodeSeg',v=>{ ST.imode=v; renderBoltz(); });
$('incNone').onclick=()=>{ ST.lines.forEach(l=>l.include=false); renderBoltz(); };
$('hideMissing').onclick=()=>{ ST.hideMissing=!ST.hideMissing; $('hideMissing').setAttribute('aria-pressed',String(ST.hideMissing)); $('hideMissing').textContent=ST.hideMissing?'Show all lines':'Hide unavailable lines'; renderBoltz(); };
$('copyManual').onclick=()=>{ const c=cond(); const P=proc(c); if(!P) return toast('No spectra loaded to copy from'); const B=boltzmann(P.x,P.nets,P.noise,ST.lines,{hw:ST.hw,energy:ST.units==='energy',mode:'area'}); let n=0; for(const r of B.rows) if(!r.missing && r.I>0){ ST.manual[r.l]={I:+r.I.toPrecision(5), ...(r.sI>0?{sI:+r.sI.toPrecision(3)}:{})}; n++; } ST.imode='manual'; segSet('imodeSeg','manual'); renderBoltz(); toast(n+' integrated intensities copied; edit them as needed'); };
$('clearManual').onclick=()=>{ ST.manual={}; renderBoltz(); };
$('pasteApply').onclick=()=>{ const txt=$('pasteTxt').value; let n=0, miss=[];
  for(const raw of txt.split(/\r?\n/)){ const s=raw.trim(); if(!s) continue; let parts=s.split(/[\t;]+|\s+/).filter(Boolean); if(parts.length<2) parts=s.split(','); if(parts.length===1) continue;
    if(parts.length>=2 && s.includes(',') && !/[\t; ]/.test(s)) parts=s.split(',');
    const nums=parts.map(t=>parseFloat(t.replace(/,$/,''))); const l=nums[0], I=nums[1], sd=nums[2]; if(!(l>0&&I>0)) continue;
    const L=ST.lines.slice().sort((p,q)=>Math.abs(p.l-l)-Math.abs(q.l-l))[0]; if(!L||Math.abs(L.l-l)>0.3){ miss.push(l); continue; }
    ST.manual[L.l]={I, ...(sd>0?{sI:sd}:{})}; L.include=true; n++; }
  ST.imode='manual'; segSet('imodeSeg','manual'); renderBoltz(); toast(`${n} value(s) applied${miss.length?'; no Ar I line near '+miss.join(', '):''}`); };
$('hw').onchange=e=>{ ST.hw=parseFloat(e.target.value)||0.35; renderBoltz(); };
$('incAll').onclick=()=>{ const c=cond(); const B=c&&boltzFor(c); ST.lines.forEach(l=>{ const r=B&&B.rows.find(q=>q.l===l.l); l.include=!!(r && !r.missing); }); renderBoltz(); };
$('incDefault').onclick=()=>{ const d=defaultLines(); ST.lines.forEach((l,i)=>l.include=d[i].include); renderBoltz(); };
$('resetConst').onclick=()=>{ const inc=ST.lines.map(l=>l.include); ST.lines=defaultLines(); ST.lines.forEach((l,i)=>l.include=inc[i]); renderBoltz(); toast('Constants reset to the built-in NIST values'); };
$('neLTE').onchange=e=>{ const v=parseFloat(e.target.value); ST.neLTE=v>0?v:null; renderBoltz(); };
$('neFromHr').onclick=()=>{ ST.neLTE=null; $('neLTE').value=''; const h=cond()&&cond().res.hr; if(!h||!Number.isFinite(h.ne)) toast('Fit the 696.5 nm line in tab 5 first'); renderBoltz(); };
$('boltzPng').onclick=()=>savePlot('plotBoltz',safe(cond()?.name||'boltzmann')+'_boltzmann.png');
$('boltzCsv').onclick=()=>{ const c=cond(), B=c&&boltzFor(c); if(!B) return; const rows=[['lambda_nm','upper','lower','Ek_eV','gk','Aki_s-1','constants_source','intensity_source','included','I','I_sd','snr','y','y_rep1','y_rep2','y_rep3']]; for(const r of B.rows) if(!r.missing) rows.push([r.l,r.up,r.low,r.E,r.g,r.A,r.src||'',ST.imode,r.include?'yes':'no',r.I,r.sI,r.snr,r.y,...[0,1,2].map(k=>r.yr?r.yr[k]:NaN)]); rows.push([]); rows.push(['T_exc_K',B.T,'sigma_K',B.sT,'R2',B.fit?B.fit.r2:NaN]); saveFile(safe(c.name)+'_boltzmann.csv',csv(rows)); };

// ===================== ROTATIONAL TAB =====================
['rotLo','rotHi','rotW'].forEach(id=>$(id).onchange=()=>{ ST.rot.lo=+$('rotLo').value; ST.rot.hi=+$('rotHi').value; ST.rot.w=+$('rotW').value; });
$('rotFixW').onchange=e=>{ ST.rot.fixW=e.target.checked; };
$('measW').onclick=()=>{ const P=proc(cond()); if(!P) return; const cands=[750.387,696.543,706.722,738.398,794.818]; for(const l of cands){ if(l<P.x[0]||l>P.x[P.x.length-1]) continue; const g=fitGaussLine(P.x,P.mean,l,0.8); if(g&&g.fwhm>0.005&&g.fwhm<2){ ST.rot.w=+g.fwhm.toFixed(4); $('rotW').value=ST.rot.w; $('rotFixW').checked=true; ST.rot.fixW=true; toast(`Instrument FWHM from Ar I ${l.toFixed(1)} nm: ${g.fwhm.toFixed(3)} nm (fixed)`); return; } } toast('No suitable Ar I line found'); };
function fitRot(c){ const P=proc(c); if(!P) return null; const lo=ST.rot.lo, hi=ST.rot.hi; if(lo<P.x[0]||hi>P.x[P.x.length-1]) return {error:'The fit window is outside the measured range.'};
  const fw=ST.rot.fixW? ST.rot.w : null;
  const m=fitN2(P.x,P.mean,lo,hi,{T:400,w:ST.rot.w},fw); if(!m) return {error:'Too few points in the fit window.'};
  const reps=P.nets.map(y=>fitN2(P.x,y,lo,hi,{T:m.T,w:m.w},fw)).filter(Boolean);
  const r={mean:m, reps, T:mean(reps.map(q=>q.T)), sdT:sdev(reps.map(q=>q.T))}; c.res.rot=r; return r; }
$('rotFit').onclick=()=>{ const c=cond(); if(!c) return; $('rotStatus').innerHTML='<div class="msg info">Fitting…</div>'; setTimeout(()=>{ fitRot(c); renderRot(); },20); };
$('rotFitAll').onclick=()=>{ $('rotStatus').innerHTML='<div class="msg info">Fitting all conditions…</div>'; setTimeout(()=>{ ST.conds.forEach(fitRot); renderRot(); toast('Fitted '+ST.conds.length+' condition(s)'); },20); };
function renderRot(){
  const c=cond(); const r=c&&c.res.rot;
  if(!c||!proc(c)){ draw('plotRot',[],baseLayout()); draw('plotProf',[],baseLayout()); $('rotKpis').innerHTML=''; $('rotStatus').innerHTML='<div class="msg info">Load plasma spectra first.</div>'; return; }
  if(!r){ const P=proc(c); const i0=idxAt(P.x,ST.rot.lo-1), i1=idxAt(P.x,ST.rot.hi+0.6);
    draw('plotRot',[{x:Array.from(P.x.slice(i0,i1)),y:Array.from(P.mean.slice(i0,i1)),mode:'lines',line:{color:css('--ink'),width:1.2},name:'Measured (mean)'}],baseLayout({xaxis:{...baseLayout().xaxis,title:{text:'Wavelength (nm)'},range:[ST.rot.lo-1,ST.rot.hi+0.6]},shapes:[{type:'rect',xref:'x',yref:'paper',x0:ST.rot.lo,x1:ST.rot.hi,y0:0,y1:1,fillcolor:alpha(css('--accent'),0.08),line:{width:0}}]}));
    draw('plotProf',[],baseLayout()); $('rotKpis').innerHTML=''; $('rotStatus').innerHTML='<div class="msg info">Shaded: fit window. Press <b>Fit this condition</b>.</div>'; return; }
  if(r.error){ $('rotStatus').innerHTML=`<div class="msg bad">${r.error}</div>`; return; }
  const m=r.mean, acc=css('--accent'), ink=css('--ink'), mut=css('--muted');
  const res=m.ys.map((v,i)=>v-m.fit[i]);
  draw('plotRot',[
    {x:m.xs,y:m.ys,mode:'lines',line:{color:ink,width:1.4},name:'Measured (mean)',xaxis:'x',yaxis:'y'},
    {x:m.xs,y:m.fit,mode:'lines',line:{color:acc,width:2,dash:'dash'},name:`Simulated, T<sub>rot</sub> = ${f(m.T,0)} K`,xaxis:'x',yaxis:'y'},
    {x:m.xs,y:res,mode:'lines',line:{color:mut,width:1},name:'Residual',xaxis:'x',yaxis:'y2',showlegend:false}],
    baseLayout({grid:{rows:2,columns:1,pattern:'shared',roworder:'top to bottom'},xaxis:{...baseLayout().xaxis,title:{text:'Wavelength (nm)'}},yaxis:{...baseLayout().yaxis,domain:[0.3,1],title:{text:'Intensity'}},yaxis2:{...baseLayout().yaxis,domain:[0,0.22],title:{text:'Resid.'}}}));
  draw('plotProf',[{x:m.profile.map(p=>p.T),y:m.profile.map(p=>p.chi2-m.cmin),mode:'lines',line:{color:acc,width:1.6},name:'Δχ²'},{x:[m.T*0.4,m.T*2.2],y:[1,1],mode:'lines',line:{color:mut,dash:'dot',width:1},name:'Δχ² = 1'}],
    baseLayout({xaxis:{...baseLayout().xaxis,title:{text:'T<sub>rot</sub> (K)'}},yaxis:{...baseLayout().yaxis,title:{text:'Δχ² (scaled)'},type:'log'},margin:{l:56,r:12,t:10,b:44},showlegend:false}));
  $('rotKpis').innerHTML=[
    ['T<sub>rot</sub> replicates (mean ± SD)', r.reps.length? `${f(r.T,0)} <small>± ${f(r.sdT,0)} K (n=${r.reps.length})</small>`:'—',1],
    ['Fit to mean spectrum', `${f(m.T,0)} <small>± ${f(m.sT,0)} K</small>`],
    ['Instrument FWHM', `${f(m.w,3)} <small>nm${ST.rot.fixW?' (fixed)':''}</small>`], ['Shift', `${f(m.shift,3)} <small>nm</small>`], ['R²', f(m.r2,4)]
  ].map(([k,v,mm])=>`<div class="kpi ${mm?'main':''}"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
  const st=[]; if(m.w>0.4) st.push('<div class="msg warn">Instrument FWHM > 0.4 nm: the rotational structure is barely resolved and T<sub>rot</sub> mostly reflects the band envelope. Expect ±20% or worse.</div>');
  if(m.r2<0.97) st.push('<div class="msg warn">The residual shows structure: check the window (keep it clear of the (1,1) head at 333.9 nm) and the background subtraction.</div>');
  st.push('<div class="msg info">± from the fit is a Δχ² = 1 interval with χ² scaled by its reduced value; the replicate SD is the better estimate of repeatability.</div>');
  $('rotStatus').innerHTML=st.join('');
}
$('rotPng').onclick=()=>savePlot('plotRot',safe(cond()?.name||'rot')+'_N2_fit.png');
$('rotCsv').onclick=()=>{ const c=cond(), r=c&&c.res.rot; if(!r||!r.mean) return toast('Fit first'); const m=r.mean; const rows=[['wavelength_nm','measured_mean','fit','residual']]; m.xs.forEach((x,i)=>rows.push([x,m.ys[i],m.fit[i],m.ys[i]-m.fit[i]])); rows.push([]); rows.push(['Trot_fit_K',m.T,'sigma_K',m.sT,'fwhm_nm',m.w]); rows.push(['Trot_reps_K',...r.reps.map(q=>q.T)]); saveFile(safe(c.name)+'_N2_rotational_fit.csv',csv(rows)); };
$('moesExp').onclick=()=>{ const c=cond(), P=proc(c); if(!P) return; const lines=['# wavelength_nm\tintensity  (mean net spectrum, '+c.name+')']; for(let i=0;i<P.x.length;i++) lines.push(P.x[i].toFixed(4)+'\t'+P.mean[i].toFixed(3)); saveFile(safe(c.name)+'_massiveOES.txt',lines.join('\n')); };

// ===================== HIGH-RES 696.5 TAB =====================
$('hrC').value=ST.hr.C;
['hrWI','hrTg','hrP','hrWS','hrC','hrCunc'].forEach(id=>$(id).onchange=()=>{ ST.hr.wI=+$('hrWI').value; ST.hr.Tg=+$('hrTg').value; ST.hr.P=+$('hrP').value; ST.hr.wS=+$('hrWS').value; ST.hr.C=+$('hrC').value; ST.hr.Cunc=+$('hrCunc').value; });
$('hrFreeG').onchange=e=>{ ST.hr.freeG=e.target.checked; };
segBind('tgSeg',v=>{ ST.hr.tgMode=v; });
function tgFor(c){ if(ST.hr.tgMode==='rot' && c.res.rot && Number.isFinite(c.res.rot.T)) return {Tg:c.res.rot.T, sTg:Number.isFinite(c.res.rot.sdT)?c.res.rot.sdT:0, src:'N₂ T<sub>rot</sub>'}; return {Tg:ST.hr.Tg, sTg:0, src:ST.hr.tgMode==='rot'?'manual (no N₂ fit yet)':'manual'}; }
function fitHr(c){ const P=procHr(c); if(!P) return null; const {Tg,sTg,src}=tgFor(c);
  const wD=7.16e-7*AR696.l*Math.sqrt(Tg/AR696.M), wG=Math.hypot(ST.hr.wI,wD);
  const l0=AR696.l; if(l0<P.grid[0]||l0>P.grid[P.grid.length-1]) return {error:'696.5 nm is outside the high-resolution file range.'};
  const fit=y=>fitVoigtLine(P.grid,y,l0,Math.min(0.45,(P.grid[P.grid.length-1]-P.grid[0])/2-0.02), ST.hr.freeG? null : wG);
  const m=fit(P.mean); const reps=P.nets.map(fit).filter(Boolean);
  const wV=ST.hr.C*ST.hr.P/Math.pow(Tg,0.7);
  const neOf=wL=>Math.max(0,(wL-wV))/ST.hr.wS*1e16;
  const neR=reps.map(q=>neOf(q.wL)); const ne=neOf(m.wL);
  const sWL = reps.length>1 ? sdev(reps.map(q=>q.wL)) : 0;
  const dVc = wV*ST.hr.Cunc/100, dVt = wV*0.7*(sTg/Tg);
  const sStark=Math.sqrt(sWL*sWL+dVc*dVc+dVt*dVt);
  const wStark=m.wL-wV;
  const r={m,reps,Tg,sTg,src,wD,wG,wV,wStark,ne,neR,sne:sStark/ST.hr.wS*1e16, upper: wStark<2*sStark, limit:2*sStark/ST.hr.wS*1e16};
  c.res.hr=r; return r; }
$('hrFit').onclick=()=>{ const c=cond(); if(!c) return; fitHr(c); renderHr(); };
function renderHr(){
  const c=cond(); const P=c&&procHr(c);
  if(!P){ draw('plotHr',[],baseLayout()); $('hrKpis').innerHTML=''; $('hrTbl').innerHTML=''; $('hrStatus').innerHTML='<div class="msg info">Add high-resolution Ar I 696.5 nm files (plasma and background) to this condition in the Data tab.</div>'; return; }
  const r=c.res.hr, acc=css('--accent'), ink=css('--ink'), mut=css('--muted');
  if(!r){ draw('plotHr',[{x:Array.from(P.grid),y:Array.from(P.mean),mode:'lines',line:{color:ink,width:1.3},name:'Mean net'}],baseLayout({xaxis:{...baseLayout().xaxis,title:{text:'Wavelength (nm)'}}})); $('hrKpis').innerHTML=''; $('hrTbl').innerHTML=''; $('hrStatus').innerHTML='<div class="msg info">Set the instrument width and gas temperature, then press <b>Fit 696.5 nm line</b>. With "From N₂ fit", run tab 4 first.</div>'; return; }
  if(r.error){ $('hrStatus').innerHTML=`<div class="msg bad">${r.error}</div>`; return; }
  const m=r.m;
  draw('plotHr',[
    {x:m.xs,y:m.ys,mode:'markers',marker:{size:4,color:ink},name:'Measured (mean)',xaxis:'x',yaxis:'y'},
    {x:m.xs,y:m.fit,mode:'lines',line:{color:acc,width:2},name:`Voigt: W<sub>L</sub> = ${f(m.wL,4)} nm`,xaxis:'x',yaxis:'y'},
    {x:m.xs,y:m.ys.map((v,i)=>v-m.fit[i]),mode:'lines',line:{color:mut,width:1},showlegend:false,xaxis:'x',yaxis:'y2'}],
    baseLayout({grid:{rows:2,columns:1,pattern:'shared'},xaxis:{...baseLayout().xaxis,title:{text:'Wavelength (nm)'}},yaxis:{...baseLayout().yaxis,domain:[0.3,1]},yaxis2:{...baseLayout().yaxis,domain:[0,0.22],title:{text:'Resid.'}}}));
  $('hrKpis').innerHTML=[
    ['nₑ (fit to mean)', r.upper? `< ${fe(r.limit)} <small>cm⁻³ (upper limit)</small>` : `${fe(r.ne,2)} <small>± ${fe(r.sne)} cm⁻³</small>`,1],
    ['nₑ replicates', r.neR.length? `${fe(mean(r.neR),2)} <small>± ${fe(sdev(r.neR))}</small>`:'—'],
    ['W<sub>Stark</sub>', `${f(r.wStark,4)} <small>nm</small>`], ['W<sub>vdW</sub>', `${f(r.wV,4)} <small>nm</small>`], ['T<sub>g</sub> used', `${f(r.Tg,0)} <small>K (${r.src})</small>`], ['Fit R²', f(m.r2,4)]
  ].map(([k,v,mm])=>`<div class="kpi ${mm?'main':''}"><span class="k">${k}</span><span class="v">${v}</span></div>`).join('');
  $('hrTbl').innerHTML=`<thead><tr><th>Width component</th><th class="num">FWHM (nm)</th><th>How obtained</th></tr></thead><tbody>
    <tr><td>Instrument W<sub>I</sub></td><td class="num">${f(ST.hr.wI,4)}</td><td>input (measure with a narrow lamp line)</td></tr>
    <tr><td>Doppler W<sub>D</sub></td><td class="num">${f(r.wD,4)}</td><td>7.16×10⁻⁷ · 696.54 · √(T<sub>g</sub>/39.95)</td></tr>
    <tr><td>Gaussian W<sub>G</sub></td><td class="num">${f(m.wG,4)}</td><td>${ST.hr.freeG?'fitted':'√(W<sub>I</sub>² + W<sub>D</sub>²), fixed'}</td></tr>
    <tr><td>Lorentzian W<sub>L</sub></td><td class="num">${f(m.wL,4)}</td><td>Voigt fit (replicate SD ${f(r.reps.length>1?sdev(r.reps.map(q=>q.wL)):NaN,4)})</td></tr>
    <tr><td>van der Waals W<sub>vdW</sub></td><td class="num">${f(r.wV,4)}</td><td>C·P/T<sub>g</sub><sup>0.7</sup>, C = ${f(ST.hr.C,3)} nm·K<sup>0.7</sup> (±${ST.hr.Cunc}%)</td></tr>
    <tr><td>Stark W<sub>Stark</sub></td><td class="num">${f(r.wStark,4)}</td><td>W<sub>L</sub> − W<sub>vdW</sub> → nₑ = 10¹⁶ · W<sub>Stark</sub>/${ST.hr.wS}</td></tr>
    <tr><td>Total Voigt W<sub>V</sub></td><td class="num">${f(m.wV,4)}</td><td>0.5346 W<sub>L</sub> + √(0.2166 W<sub>L</sub>² + W<sub>G</sub>²)</td></tr></tbody>`;
  const st=[]; if(r.upper) st.push(`<div class="msg warn">The Stark width (${f(r.wStark,4)} nm) is smaller than twice its uncertainty (${f(r.limit*ST.hr.wS/1e16,4)} nm, from replicate scatter, the ±${ST.hr.Cunc}% vdW coefficient and the T<sub>g</sub> uncertainty), so only an upper limit is meaningful.</div>`);
  if(ST.hr.wI>0.06) st.push('<div class="msg warn">W<sub>I</sub> > 0.06 nm: the Lorentzian part is poorly constrained. A higher-resolution grating or a narrower slit is recommended.</div>');
  st.push(`<div class="msg info">This nₑ is used in the McWhirter test on tab 3 (unless you typed a value there).</div>`);
  $('hrStatus').innerHTML=st.join('');
}
$('hrPng').onclick=()=>savePlot('plotHr',safe(cond()?.name||'ar696')+'_Ar696_fit.png');
$('hrCsv').onclick=()=>{ const c=cond(), r=c&&c.res.hr; if(!r||!r.m) return toast('Fit first'); const m=r.m; const rows=[['wavelength_nm','measured_mean','voigt_fit']]; m.xs.forEach((x,i)=>rows.push([x,m.ys[i],m.fit[i]])); rows.push([]); rows.push(['W_L_nm',m.wL,'W_G_nm',m.wG,'W_vdW_nm',r.wV,'W_Stark_nm',r.wStark,'ne_cm-3',r.ne,'sigma_ne',r.sne,'Tg_K',r.Tg]); saveFile(safe(c.name)+'_Ar696_broadening.csv',csv(rows)); };

// ===================== RESULTS TAB =====================
function ratiosFor(c){ const P=proc(c); if(!P) return null; const out={}; const ref=P.nets.map(y=>windowArea(P.x,y,749.95,750.85));
  for(const F of RATIO_FEATURES){ if(F.key==='Ar') continue; if(F.lo<P.x[0]||F.hi>P.x[P.x.length-1]){ out[F.key]=null; continue; } const v=P.nets.map((y,i)=>windowArea(P.x,y,F.lo,F.hi)/ref[i]); out[F.key]={m:mean(v), s:sdev(v)}; } return out; }
function renderRes(){
  if(!ST.conds.length){ $('resTbl').innerHTML='<tbody><tr><td>No conditions loaded.</td></tr></tbody>'; draw('plotResT',[],baseLayout()); draw('plotResR',[],baseLayout()); return; }
  const rows=ST.conds.map(c=>{ const B=boltzFor(c); if(B) c.res.boltz=B; const L=B?lteChecks(B,{ne:neForLTE(c)}):null; return {c,B,L,R:c.res.rot,H:c.res.hr,Q:ratiosFor(c)}; });
  const pm=(m,s,d=0)=>Number.isFinite(m)? `${f(m,d)}${Number.isFinite(s)?' ± '+f(s,d):''}`:'—';
  $('resTbl').innerHTML=`<thead><tr><th>Condition</th><th class="num">T<sub>exc</sub> fit (K)</th><th class="num">T<sub>exc</sub> reps (K)</th><th>LTE check</th><th class="num">T<sub>rot</sub> N₂ (K)</th><th class="num">nₑ (cm⁻³)</th>${RATIO_FEATURES.filter(F=>F.key!=='Ar').map(F=>`<th class="num">${F.label.split(' ')[0]==='N₂'?F.label.split(' ').slice(0,2).join(' '):F.label.split(' ')[0]} / Ar750</th>`).join('')}</tr></thead><tbody>${rows.map(r=>`<tr><td>${esc(r.c.name)}</td>
    <td class="num">${r.B?pm(r.B.T,r.B.sT):'—'}</td><td class="num">${r.B&&r.B.Treps&&r.B.Treps.length?pm(mean(r.B.Treps),sdev(r.B.Treps)):'—'}</td>
    <td>${r.L? (r.L.verdict==='consistent'?'<span style="color:var(--good)">consistent</span>':'<span style="color:var(--bad)">non-LTE signs</span>'):'—'}</td>
    <td class="num">${r.R&&r.R.reps?pm(r.R.T,r.R.sdT):'—'}</td><td class="num">${r.H&&r.H.m? (r.H.upper? '< '+fe(r.H.limit) : fe(r.H.ne,2)) :'—'}</td>
    ${RATIO_FEATURES.filter(F=>F.key!=='Ar').map(F=>`<td class="num">${r.Q&&r.Q[F.key]?pm(r.Q[F.key].m,r.Q[F.key].s,3):'—'}</td>`).join('')}</tr>`).join('')}</tbody>`;
  const names=rows.map(r=>r.c.name), acc=css('--accent'), sec='#c77a12';
  draw('plotResT',[
    {x:names,y:rows.map(r=>r.B&&r.B.Treps&&r.B.Treps.length?mean(r.B.Treps):null),error_y:{type:'data',array:rows.map(r=>r.B&&r.B.Treps?sdev(r.B.Treps):0),visible:true},type:'bar',name:'T<sub>exc</sub> (Ar I)',marker:{color:acc},offsetgroup:1},
    {x:names,y:rows.map(r=>r.R&&r.R.reps?r.R.T:null),error_y:{type:'data',array:rows.map(r=>r.R?r.R.sdT:0),visible:true},type:'bar',name:'T<sub>rot</sub> (N₂)',marker:{color:sec},yaxis:'y2',offsetgroup:2}],
    baseLayout({barmode:'group',yaxis:{...baseLayout().yaxis,title:{text:'T<sub>exc</sub> (K)'}},yaxis2:{...baseLayout().yaxis,title:{text:'T<sub>rot</sub> (K)'},overlaying:'y',side:'right',showgrid:false},margin:{l:62,r:62,t:30,b:70}}));
  const keys=RATIO_FEATURES.filter(F=>F.key!=='Ar');
  const pal=['#c77a12','#0e8fa3','#b04a91','#d4453b','#e06b5f','#2f9e6b'];
  draw('plotResR',keys.map((F,k)=>({x:names,y:rows.map(r=>r.Q&&r.Q[F.key]?r.Q[F.key].m:null),error_y:{type:'data',array:rows.map(r=>r.Q&&r.Q[F.key]?r.Q[F.key].s:0),visible:true},type:'bar',name:F.label.replace(' (reference)',''),marker:{color:pal[k%pal.length]}})),
    baseLayout({barmode:'group',yaxis:{...baseLayout().yaxis,title:{text:'Intensity / Ar I 750.4'},type:'log'},margin:{l:62,r:16,t:30,b:70}}));
  ST._rows=rows;
}
$('resRun').onclick=()=>{ $('resRun').textContent='Computing…'; setTimeout(()=>{ ST.conds.forEach(c=>{ fitRot(c); if(c.hr.plasma.length) fitHr(c); }); renderRes(); $('resRun').textContent='Recompute all (incl. fits)'; toast('All conditions recomputed'); },20); };
$('resCsv').onclick=()=>{ if(!ST._rows) renderRes(); const rows=[['condition','Texc_fit_K','Texc_fit_sigma_K','Texc_reps_mean_K','Texc_reps_sd_K','lines_used','Ek_span_eV','R2','LTE_check','Trot_reps_mean_K','Trot_reps_sd_K','Trot_fit_K','ne_cm-3','ne_sigma','ne_is_upper_limit',...RATIO_FEATURES.filter(F=>F.key!=='Ar').flatMap(F=>[F.key+'_over_Ar750',F.key+'_sd'])]];
  for(const r of ST._rows){ const B=r.B; rows.push([r.c.name,B?B.T:NaN,B?B.sT:NaN,B&&B.Treps?mean(B.Treps):NaN,B&&B.Treps?sdev(B.Treps):NaN,B?B.use.length:NaN,B?B.Espan:NaN,B&&B.fit?B.fit.r2:NaN,r.L?r.L.verdict:'',r.R&&r.R.reps?r.R.T:NaN,r.R&&r.R.reps?r.R.sdT:NaN,r.R&&r.R.mean?r.R.mean.T:NaN,r.H&&r.H.m?r.H.ne:NaN,r.H?r.H.sne:NaN,r.H?(r.H.upper?'yes':'no'):'',...RATIO_FEATURES.filter(F=>F.key!=='Ar').flatMap(F=>r.Q&&r.Q[F.key]?[r.Q[F.key].m,r.Q[F.key].s]:[NaN,NaN])]); }
  saveFile('OES_summary.csv',csv(rows)); };

// ===================== tabs & refresh =====================
const TABS=['data','spec','boltz','rot','hr','res'];
$('tabs').addEventListener('click',e=>{ const b=e.target.closest('[data-t]'); if(!b) return; ST.tab=b.dataset.t; for(const x of $('tabs').children) x.setAttribute('aria-selected',String(x===b)); renderTab(); });
function renderTab(){ for(const t of TABS) $('t-'+t).hidden = t!==ST.tab; ({data:renderData,spec:renderSpec,boltz:renderBoltz,rot:renderRot,hr:renderHr,res:renderRes})[ST.tab](); if(hasPlotly()) setTimeout(()=>document.querySelectorAll('#t-'+ST.tab+' .js-plotly-plot').forEach(g=>Plotly.Plots.resize(g)),30); }
function refreshAll(){ renderCondSel(); renderTab(); }
const rerenderTheme=()=>{ if(ST.tab!=='data') renderTab(); };
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',rerenderTheme);
new MutationObserver(rerenderTheme).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});

// start in a working state with the labelled synthetic example
$('loadEx').click(); ST.tab='spec'; for(const x of $('tabs').children) x.setAttribute('aria-selected',String(x.dataset.t==='spec')); renderTab();
})();
