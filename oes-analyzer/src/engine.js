// ===================== OES engine (pure functions) =====================
const KB_EV = 8.617333262e-5;   // Boltzmann constant, eV/K
const C2 = 1.4387769;           // hc/k, cm·K

// ---------- small numerics ----------
function mean(a){ let s=0; for(const v of a) s+=v; return a.length? s/a.length : NaN; }
function sdev(a){ if(a.length<2) return NaN; const m=mean(a); let s=0; for(const v of a) s+=(v-m)*(v-m); return Math.sqrt(s/(a.length-1)); }
function median(a){ const b=Array.from(a).filter(Number.isFinite).sort((p,q)=>p-q); if(!b.length) return NaN; const k=b.length>>1; return b.length%2? b[k] : 0.5*(b[k-1]+b[k]); }
function interp(x,y,xq){ // x ascending
  const n=x.length; if(xq<=x[0]) return y[0]; if(xq>=x[n-1]) return y[n-1];
  let lo=0, hi=n-1; while(hi-lo>1){ const m=(lo+hi)>>1; if(x[m]<=xq) lo=m; else hi=m; }
  const t=(xq-x[lo])/(x[hi]-x[lo]); return y[lo]+(y[hi]-y[lo])*t; }
function idxAt(x,xq){ let lo=0, hi=x.length-1; if(xq<=x[0]) return 0; if(xq>=x[hi]) return hi; while(hi-lo>1){ const m=(lo+hi)>>1; if(x[m]<=xq) lo=m; else hi=m; } return (xq-x[lo] < x[hi]-xq) ? lo : hi; }
// Edlén (1966) refractive index of standard air, λ in nm (vacuum)
function nAir(lnm){ const s2=Math.pow(1000/lnm,2); return 1+1e-8*(8342.54+2406147/(130-s2)+15998/(38.9-s2)); }
const vac2air = l => l/nAir(l);

// weighted linear regression y = a + b x ; returns {a,b,sa,sb,r2,n,chi2red}
function linfit(x,y,w){
  const n=x.length; if(n<2) return null; w = w || x.map(()=>1);
  let S=0,Sx=0,Sy=0,Sxx=0,Sxy=0; for(let i=0;i<n;i++){ S+=w[i]; Sx+=w[i]*x[i]; Sy+=w[i]*y[i]; Sxx+=w[i]*x[i]*x[i]; Sxy+=w[i]*x[i]*y[i]; }
  const D=S*Sxx-Sx*Sx; if(Math.abs(D)<1e-300) return null;
  const b=(S*Sxy-Sx*Sy)/D, a=(Sxx*Sy-Sx*Sxy)/D;
  let ss=0, st=0; const ym=Sy/S; for(let i=0;i<n;i++){ const r=y[i]-a-b*x[i]; ss+=w[i]*r*r; st+=w[i]*(y[i]-ym)*(y[i]-ym); }
  const dof=Math.max(1,n-2), s2=ss/dof; // scale errors by reduced chi2 (robust when weights are relative)
  return { a, b, sa:Math.sqrt(s2*Sxx/D), sb:Math.sqrt(s2*S/D), r2: st>0? 1-ss/st : NaN, n, chi2red:s2 };
}
// solve small linear system (Gaussian elimination)
function solve(A,b){ const n=b.length; const M=A.map((r,i)=>[...r,b[i]]);
  for(let c=0;c<n;c++){ let p=c; for(let r=c+1;r<n;r++) if(Math.abs(M[r][c])>Math.abs(M[p][c])) p=r; [M[c],M[p]]=[M[p],M[c]];
    if(Math.abs(M[c][c])<1e-300) return null; for(let r=c+1;r<n;r++){ const f=M[r][c]/M[c][c]; for(let k=c;k<=n;k++) M[r][k]-=f*M[c][k]; } }
  const x=new Array(n).fill(0); for(let r=n-1;r>=0;r--){ let s=M[r][n]; for(let k=r+1;k<n;k++) s-=M[r][k]*x[k]; x[r]=s/M[r][r]; } return x; }
// least squares y ≈ Σ c_k basis_k  (basis: array of arrays)
function lstsq(basis,y){ const m=basis.length, n=y.length; const A=[...Array(m)].map(()=>new Array(m).fill(0)), b=new Array(m).fill(0);
  for(let i=0;i<n;i++){ for(let p=0;p<m;p++){ b[p]+=basis[p][i]*y[i]; for(let q=p;q<m;q++) A[p][q]+=basis[p][i]*basis[q][i]; } }
  for(let p=0;p<m;p++) for(let q=0;q<p;q++) A[p][q]=A[q][p]; return solve(A,b); }
// Nelder–Mead minimiser
function nelderMead(f,x0,step,opts={}){ const n=x0.length, maxIt=opts.maxIt||400, tol=opts.tol||1e-9;
  let S=[x0.slice()]; for(let i=0;i<n;i++){ const p=x0.slice(); p[i]+=step[i]; S.push(p); }
  let F=S.map(f);
  for(let it=0;it<maxIt;it++){
    const ord=F.map((v,i)=>i).sort((a,b)=>F[a]-F[b]); S=ord.map(i=>S[i]); F=ord.map(i=>F[i]);
    if(Math.abs(F[n]-F[0])<=tol*(Math.abs(F[0])+1e-30)) break;
    const c=new Array(n).fill(0); for(let i=0;i<n;i++) for(let k=0;k<n;k++) c[k]+=S[i][k]/n;
    const xr=c.map((v,k)=>v+(v-S[n][k])), fr=f(xr);
    if(fr<F[0]){ const xe=c.map((v,k)=>v+2*(v-S[n][k])), fe=f(xe); if(fe<fr){S[n]=xe;F[n]=fe;} else {S[n]=xr;F[n]=fr;} }
    else if(fr<F[n-1]){ S[n]=xr; F[n]=fr; }
    else { const xc=c.map((v,k)=>v+0.5*(S[n][k]-v)), fc=f(xc);
      if(fc<F[n]){ S[n]=xc; F[n]=fc; } else { for(let i=1;i<=n;i++){ S[i]=S[i].map((v,k)=>S[0][k]+0.5*(v-S[0][k])); F[i]=f(S[i]); } } }
  }
  let bi=0; for(let i=1;i<=n;i++) if(F[i]<F[bi]) bi=i; return {x:S[bi], f:F[bi]};
}

// ---------- file parsing ----------
function parseSpectrum(text){
  const xs=[], ys=[];
  for(const raw of text.split(/\r?\n/)){
    const s=raw.trim(); if(!s) continue;
    let parts=s.split(/[\t;]+/).map(t=>t.trim()).filter(Boolean), commaDec=true;
    if(parts.length<2){ parts=s.split(/\s+/); }
    if(parts.length<2){ parts=s.split(','); commaDec=false; }
    if(parts.length<2) continue;
    const num=t=>{ t=t.replace(/^"|"$/g,''); if(commaDec && /^-?\d+,\d+([eE][-+]?\d+)?$/.test(t)) t=t.replace(',','.'); return /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(t)? parseFloat(t) : NaN; };
    const a=num(parts[0]), b=num(parts[1]);
    if(Number.isFinite(a)&&Number.isFinite(b)){ xs.push(a); ys.push(b); }
  }
  if(xs.length<20) throw new Error('fewer than 20 numeric two-column rows were found');
  const idx=xs.map((_,i)=>i).sort((p,q)=>xs[p]-xs[q]);
  const x=[], y=[]; for(const i of idx){ if(x.length && xs[i]===x[x.length-1]) continue; x.push(xs[i]); y.push(ys[i]); }
  return { x:Float64Array.from(x), y:Float64Array.from(y) };
}

// ---------- line database ----------
// Ar I lines for the Boltzmann plot. Air wavelengths (nm), A_ki (s^-1), g_k, E_k (eV), lower level (Paschen).
// Values as tabulated in the NIST Atomic Spectra Database (editable in the app; check against NIST before publishing).
const AR_BOLTZ_DEFAULT = [
  // 5p → 4s (Paschen 3p → 1s)
  {l:404.442,A:3.13e5,g:5,E:14.688,low:'1s5',up:'3p5 5p',src:'verify'},
  {l:415.859,A:1.40e6,g:5,E:14.529,low:'1s5',up:'3p5 5p',src:'NIST'},
  {l:416.418,A:2.88e5,g:3,E:14.529,low:'1s5',up:'3p5 5p',src:'verify'},
  {l:418.188,A:5.61e5,g:1,E:14.693,low:'1s3',up:'3p5 5p',src:'verify'},
  {l:419.071,A:2.80e5,g:5,E:14.506,low:'1s5',up:'3p5 5p',src:'verify'},
  {l:419.103,A:5.39e5,g:3,E:14.682,low:'1s3',up:'3p5 5p',src:'verify'},
  {l:419.832,A:2.57e6,g:1,E:14.576,low:'1s4',up:'3p5 5p',src:'NIST'},
  {l:420.068,A:9.67e5,g:7,E:14.499,low:'1s5',up:'3p5 5p',src:'NIST'},
  {l:425.118,A:1.11e5,g:3,E:14.464,low:'1s4',up:'3p5 5p',src:'verify'},
  {l:425.936,A:3.98e6,g:1,E:14.738,low:'1s2',up:'3p5 5p',src:'NIST'},
  {l:426.629,A:3.12e5,g:5,E:14.528,low:'1s4',up:'3p5 5p',src:'verify'},
  {l:427.217,A:7.97e5,g:3,E:14.525,low:'1s4',up:'3p5 5p',src:'NIST'},
  {l:430.010,A:3.77e5,g:5,E:14.506,low:'1s4',up:'3p5 5p',src:'NIST'},
  {l:433.356,A:5.68e5,g:5,E:14.688,low:'1s2',up:'3p5 5p',src:'verify'},
  {l:433.534,A:3.87e5,g:3,E:14.686,low:'1s2',up:'3p5 5p',src:'verify'},
  {l:434.517,A:2.97e5,g:3,E:14.682,low:'1s2',up:'3p5 5p',src:'verify'},
  // higher levels (6d, 7s, 5d) → 4p (Paschen 2p)
  {l:545.165,A:4.70e5,g:7,E:15.180,low:'2p10',up:'7s',src:'Ohno 2006'},
  {l:549.587,A:1.69e6,g:9,E:15.331,low:'2p9',up:'6d',src:'Ohno 2006'},
  {l:555.870,A:1.42e6,g:5,E:14.955,low:'2p8',up:'6d',src:'Ohno 2006'},
  {l:560.673,A:2.20e6,g:3,E:14.742,low:'2p10',up:'6d',src:'Ohno 2006'},
  {l:565.070,A:3.20e6,g:1,E:15.003,low:'2p6',up:'6d',src:'Ohno 2006'},
  {l:603.213,A:2.46e6,g:9,E:15.130,low:'2p9',up:'5d',src:'verify'},
  {l:641.631,A:1.16e6,g:5,E:14.839,low:'2p10',up:'5d',src:'verify'},
  {l:675.283,A:1.93e6,g:5,E:14.743,low:'2p8',up:'5d',src:'verify'},
  {l:687.129,A:2.78e6,g:3,E:14.711,low:'2p8',up:'5d',src:'verify'},
  {l:703.025,A:2.67e6,g:5,E:14.839,low:'2p7',up:'5d',src:'verify'},
  // 4p → 4s (Paschen 2p → 1s)
  {l:696.543,A:6.39e6,g:3,E:13.328,low:'1s5',up:'2p2',src:'NIST'},
  {l:706.722,A:3.80e6,g:5,E:13.302,low:'1s5',up:'2p3',src:'NIST'},
  {l:714.704,A:6.25e5,g:3,E:13.283,low:'1s5',up:'2p4',src:'NIST'},
  {l:727.294,A:1.83e6,g:3,E:13.328,low:'1s4',up:'2p2',src:'NIST'},
  {l:738.398,A:8.47e6,g:5,E:13.302,low:'1s4',up:'2p3',src:'NIST'},
  {l:750.387,A:4.45e7,g:1,E:13.480,low:'1s2',up:'2p1',src:'NIST'},
  {l:751.465,A:4.02e7,g:1,E:13.273,low:'1s4',up:'2p5',src:'NIST'},
  {l:763.511,A:2.45e7,g:5,E:13.172,low:'1s5',up:'2p6',src:'NIST'},
  {l:772.376,A:5.18e6,g:3,E:13.153,low:'1s5',up:'2p7',src:'NIST'},
  {l:772.421,A:1.17e7,g:3,E:13.328,low:'1s3',up:'2p2',src:'verify'},
  {l:794.818,A:1.86e7,g:3,E:13.283,low:'1s3',up:'2p4',src:'NIST'},
  {l:800.616,A:4.90e6,g:5,E:13.172,low:'1s4',up:'2p6',src:'NIST'},
  {l:801.479,A:9.28e6,g:5,E:13.095,low:'1s5',up:'2p8',src:'NIST'},
  {l:810.369,A:2.50e7,g:3,E:13.153,low:'1s4',up:'2p7',src:'NIST'},
  {l:811.531,A:3.31e7,g:7,E:13.076,low:'1s5',up:'2p9',src:'NIST'},
  {l:826.452,A:1.53e7,g:3,E:13.328,low:'1s2',up:'2p2',src:'NIST'},
  {l:840.821,A:2.23e7,g:5,E:13.302,low:'1s2',up:'2p3',src:'NIST'},
  {l:842.465,A:2.15e7,g:5,E:13.095,low:'1s4',up:'2p8',src:'NIST'},
  {l:852.144,A:1.39e7,g:3,E:13.283,low:'1s2',up:'2p4',src:'NIST'},
  {l:866.794,A:2.43e6,g:3,E:13.153,low:'1s3',up:'2p7',src:'verify'}
];
const METASTABLE = new Set(['1s5','1s3']);

// identification list: species, λ (air nm), kind ('atom' sharp line, 'head' band head degraded to violet, 'band' broad feature), weight (likelihood)
const ID_LINES = [
  ...AR_BOLTZ_DEFAULT.map(d=>({sp:'Ar I', l:d.l, kind:'atom', w:3})),
  ...[416.418,418.188,419.071,419.103,425.118,426.629,433.356,433.534,434.517,603.213,641.631,667.728,675.283,687.129,703.025,710.748,731.172,735.329,737.212,741.234,743.536,794.818,866.794,867.032,935.422,965.778].map(l=>({sp:'Ar I', l, kind:'atom', w:2})),
  ...[434.806,442.600,454.505,460.957,472.687,476.487,480.602,487.986].map(l=>({sp:'Ar II', l, kind:'atom', w:1})),
  {sp:'H I (Hα)', l:656.279, kind:'atom', w:4}, {sp:'H I (Hβ)', l:486.133, kind:'atom', w:4}, {sp:'H I (Hγ)', l:434.047, kind:'atom', w:2},
  {sp:'O I', l:777.417, kind:'atom', w:4, note:'777.19/777.42/777.54 triplet'}, {sp:'O I', l:844.636, kind:'atom', w:3, note:'844.62/844.64/844.68 triplet'},
  {sp:'Na I', l:588.995, kind:'atom', w:1}, {sp:'Na I', l:589.592, kind:'atom', w:1},
  {sp:'OH A–X (0,0)', l:306.36, kind:'head', w:3, note:'R₂ band head'}, {sp:'OH A–X (0,0)', l:308.90, kind:'band', w:4, note:'Q₁ branch peak'},
  {sp:'OH A–X (0,0)', l:307.84, kind:'band', w:2, note:'R branch'}, {sp:'OH A–X (1,1)', l:312.20, kind:'head', w:2},
  {sp:'N₂ SPS (2,1)', l:313.60, kind:'head', w:2}, {sp:'N₂ SPS (1,0)', l:315.93, kind:'head', w:3},
  {sp:'N₂ SPS (2,2)', l:330.99, kind:'head', w:1}, {sp:'N₂ SPS (1,1)', l:333.89, kind:'head', w:2}, {sp:'N₂ SPS (0,0)', l:337.13, kind:'head', w:5},
  {sp:'N₂ SPS (2,3)', l:350.05, kind:'head', w:1}, {sp:'N₂ SPS (1,2)', l:353.67, kind:'head', w:2}, {sp:'N₂ SPS (0,1)', l:357.69, kind:'head', w:4},
  {sp:'N₂ SPS (2,4)', l:371.05, kind:'head', w:2}, {sp:'N₂ SPS (1,3)', l:375.54, kind:'head', w:3}, {sp:'N₂ SPS (0,2)', l:380.49, kind:'head', w:4},
  {sp:'N₂ SPS (2,5)', l:394.30, kind:'head', w:1}, {sp:'N₂ SPS (1,4)', l:399.84, kind:'head', w:2}, {sp:'N₂ SPS (0,3)', l:405.94, kind:'head', w:2},
  {sp:'N₂⁺ FNS (0,0)', l:391.44, kind:'head', w:4}, {sp:'N₂⁺ FNS (0,1)', l:427.81, kind:'head', w:2}, {sp:'N₂⁺ FNS (1,0)', l:358.21, kind:'head', w:1}
];

// ---------- processing of one condition ----------
// files: {x,y}; returns grid, nets (per replicate), mean, sd, bg mean/sd, noise σ, saturation info
function processCondition(cond, opts){
  const P=cond.plasma, B=cond.bg; if(!P.length) return null;
  const grid=P[0].x;
  const on=f=>Float64Array.from(grid, v=>interp(f.x,f.y,v));
  const pr=P.map(on), br=B.map(on);
  const scale = (cond.expP>0 && cond.expB>0) ? cond.expP/cond.expB : 1;
  const n=grid.length, bgMean=new Float64Array(n), bgSD=new Float64Array(n);
  if(br.length){ for(let i=0;i<n;i++){ const v=br.map(r=>r[i]); bgMean[i]=mean(v)*scale; bgSD[i]=br.length>1? sdev(v)*scale : 0; } }
  const nets=pr.map(r=>Float64Array.from(r,(v,i)=>v-bgMean[i]));
  const m=new Float64Array(n), s=new Float64Array(n);
  for(let i=0;i<n;i++){ const v=nets.map(r=>r[i]); m[i]=mean(v); s[i]=nets.length>1? sdev(v):0; }
  // noise: background replicate scatter if available, else robust high-frequency estimate of the mean spectrum
  let noise = br.length>1 ? median(bgSD) : NaN;
  if(!(noise>0)){ const d=[]; for(let i=1;i<n;i++) d.push(m[i]-m[i-1]); const md=median(d); noise=1.4826*median(d.map(v=>Math.abs(v-md)))/Math.SQRT2; }
  if(!(noise>0)) noise=1;
  const sat=opts.sat||Infinity; const satFiles=[];
  [...P.map((f,i)=>['plasma '+(i+1),f]), ...B.map((f,i)=>['background '+(i+1),f])].forEach(([nm,f])=>{ let c=0; for(const v of f.y) if(v>=sat) c++; if(c) satFiles.push(`${nm}: ${c} pixels`); });
  const step=median(Array.from(grid).slice(1).map((v,i)=>v-grid[i]));
  return { grid, nets, mean:m, sd:s, bgMean, bgSD, noise, satFiles, step, nRep:nets.length };
}

// wavelength shift from strong, isolated Ar I lines (measured − reference), nm
function autoShift(x,y){
  const refs=[696.543,706.722,738.398,750.387,763.511,794.818,811.531,842.465,852.144,912.297];
  const off=[]; let maxv=0; for(const v of y) if(v>maxv) maxv=v;
  for(const r of refs){ if(r<x[0]+1||r>x[x.length-1]-1) continue;
    const i0=idxAt(x,r-0.8), i1=idxAt(x,r+0.8); let ib=i0; for(let i=i0;i<=i1;i++) if(y[i]>y[ib]) ib=i;
    if(y[ib]<0.02*maxv) continue;
    let sw=0, sx=0; for(let i=Math.max(0,ib-3);i<=Math.min(x.length-1,ib+3);i++){ const w=Math.max(0,y[i]-0.3*y[ib]); sw+=w; sx+=w*x[i]; }
    if(sw>0) off.push(sx/sw-r); }
  return off.length>=2 ? {shift:median(off), n:off.length, spread:sdev(off)} : {shift:0, n:off.length, spread:NaN};
}

// peak finding on (x,y) with noise σ
function findPeaks(x,y,noise,opts){
  const thr=opts.thr*noise, n=x.length, step=(x[n-1]-x[0])/(n-1);
  const sm=Float64Array.from(y,(v,i)=>(y[Math.max(0,i-1)]+2*v+y[Math.min(n-1,i+1)])/4);
  const win=Math.max(3,Math.round(1.5/step));
  const out=[];
  for(let i=1;i<n-1;i++){
    if(!(sm[i]>sm[i-1]&&sm[i]>=sm[i+1]&&sm[i]>thr)) continue;
    let lmin=sm[i], rmin=sm[i]; for(let k=Math.max(0,i-win);k<i;k++) lmin=Math.min(lmin,sm[k]); for(let k=i+1;k<=Math.min(n-1,i+win);k++) rmin=Math.min(rmin,sm[k]);
    const prom=sm[i]-Math.max(lmin,rmin); if(prom<thr) continue;
    // sub-pixel centre by parabola
    const a=sm[i-1], b=sm[i], c=sm[i+1], den=a-2*b+c, dx= den!==0 ? 0.5*(a-c)/den : 0;
    out.push({i, l:x[i]+dx*step, h:y[i], prom, snr:y[i]/noise});
  }
  out.sort((p,q)=>q.h-p.h);
  const kept=[]; for(const p of out){ if(kept.every(k=>Math.abs(k.l-p.l)>=opts.sep)) kept.push(p); }
  kept.sort((p,q)=>p.l-q.l);
  for(const p of kept){
    let best=null, score=-Infinity;
    for(const L of ID_LINES){ const d=p.l-L.l; const tol = L.kind==='atom'? opts.tol : opts.tol*1.6; if(Math.abs(d)>tol) continue;
      const sc=L.w-3*Math.abs(d)/tol; if(sc>score){ score=sc; best=L; } }
    p.id=best;
  }
  return kept;
}

// integrated line intensity with a local linear baseline; returns per-replicate areas
function lineArea(x,y,l0,hw){
  const i0=idxAt(x,l0-hw), i1=idxAt(x,l0+hw); if(i1-i0<2) return NaN;
  const e=Math.max(1,Math.min(3,Math.floor((i1-i0)/4)));
  let bl=0, nl=0, br=0, nr=0; for(let k=Math.max(0,i0-e);k<=i0;k++){ bl+=y[k]; nl++; } for(let k=i1;k<=Math.min(x.length-1,i1+e);k++){ br+=y[k]; nr++; }
  bl/=nl; br/=nr; const xa=x[i0], xb=x[i1];
  let A=0; for(let k=i0;k<i1;k++){ const f=(xx,yy)=>yy-(bl+(br-bl)*(xx-xa)/(xb-xa)); A+=0.5*(f(x[k],y[k])+f(x[k+1],y[k+1]))*(x[k+1]-x[k]); }
  return A;
}

// Boltzmann plot. lines: [{l,A,g,E,low,include}], spectra: array of y (replicates) on x (shift-corrected), resp: fn(λ)->relative response or null
function boltzmann(x, reps, noise, lines, opts){
  const rows=[], mode=opts.mode||'area', lamF=d=>opts.energy? d.l : 1;
  const tol=(opts.tol||0.25);
  for(const d of lines){
    const base={...d, meta:METASTABLE.has(d.low)};
    base.blend=ID_LINES.filter(L=>L.l!==d.l && Math.abs(L.l-d.l)<opts.hw*1.05 && (L.sp==='Ar I'||L.sp.startsWith('N₂')||L.sp.startsWith('H'))).map(L=>`${L.sp} ${L.l.toFixed(2)}`);
    if(mode==='manual'){
      const m=opts.manual&&opts.manual[d.l];
      if(!m||!(m.I>0)){ rows.push({...base, missing:true, why:'no manual value'}); continue; }
      const yv=Math.log(m.I*lamF(d)/(d.g*d.A));
      rows.push({...base, I:m.I, sI:m.sI>0?m.sI:NaN, snr:NaN, y:yv, yr:[], sy:m.sI>0? m.sI/m.I : NaN, manual:true}); continue;
    }
    if(d.l<x[0]+opts.hw+0.2 || d.l>x[x.length-1]-opts.hw-0.2){ rows.push({...base, missing:true, why:'outside measured range'}); continue; }
    let center=d.l, pk=null;
    if(mode==='peakH'||mode==='peakA'){
      pk=(opts.peaks||[]).filter(p=>Math.abs(p.l-d.l)<=tol).sort((p,q)=>Math.abs(p.l-d.l)-Math.abs(q.l-d.l))[0];
      if(!pk){ rows.push({...base, missing:true, why:'not found by peak finder'}); continue; }
      center=pk.l;
    }
    const R = opts.resp ? opts.resp(d.l) : 1;
    const heightOf=y=>{ const i0=idxAt(x,center-opts.hw), i1=idxAt(x,center+opts.hw), ic=idxAt(x,center); let m=-Infinity; for(let k=Math.max(0,ic-2);k<=Math.min(x.length-1,ic+2);k++) m=Math.max(m,y[k]); const bl=0.5*(y[i0]+y[i1]); return m-bl; };
    const vals=reps.map(y=>(mode==='peakH'? heightOf(y) : lineArea(x,y,center,opts.hw))/R);
    const I=mean(vals), sI=sdev(vals);
    const ic=idxAt(x,center); const pkv=mean(reps.map(y=>{ let m=-Infinity; for(let k=Math.max(0,ic-3);k<=Math.min(x.length-1,ic+3);k++) m=Math.max(m,y[k]); return m; }));
    const yv = I>0 ? Math.log(I*lamF(d)/(d.g*d.A)) : NaN;
    const yr = vals.map(a=> a>0 ? Math.log(a*lamF(d)/(d.g*d.A)) : NaN);
    rows.push({...base, I, sI, snr:pkv/noise, y:yv, yr, sy:(I>0&&sI>0)? sI/I : NaN, center});
  }
  const use=rows.filter(r=>r.include && !r.missing && Number.isFinite(r.y));
  const res={rows, use};
  if(use.length>=2){
    const w=use.map(r=> Number.isFinite(r.sy)? 1/Math.max(r.sy,0.02)**2 : 1/0.05**2);
    const f=linfit(use.map(r=>r.E), use.map(r=>r.y), w);
    if(f && f.b<0){ res.fit=f; res.T=-1/(KB_EV*f.b); res.sT=f.sb/(KB_EV*f.b*f.b); }
    else if(f){ res.fit=f; res.T=NaN; }
    const nrep=reps.length, Tr=[];
    for(let k=0;k<nrep;k++){ const ix=use.map((r,i)=>i).filter(i=>use[i].yr && Number.isFinite(use[i].yr[k])); if(ix.length<2) continue; const g=linfit(ix.map(i=>use[i].E), ix.map(i=>use[i].yr[k]), ix.map(i=>w[i])); if(g && g.b<0) Tr.push(-1/(KB_EV*g.b)); }
    res.Treps=Tr;
    const Es=use.map(r=>r.E); res.Espan=Math.max(...Es)-Math.min(...Es);
  }
  return res;
}

// ---------- N2 C3Πu – B3Πg (0,0) simplified rotational model ----------
// Singlet-like Π–Π approximation (triplet fine structure unresolved). Constants: Huber & Herzberg / Lofthus & Krupenie.
const N2SPS = { Bu:1.8154, Du:6.4e-6, Bl:1.6285, Dl:5.9e-6, headAir:337.13, Jmax:70 };
(function(){ // set band origin so that the P-branch head falls at 337.13 nm (air)
  const F=(B,D,J)=>B*J*(J+1)-D*J*J*(J+1)*(J+1);
  const headVac=N2SPS.headAir*nAir(N2SPS.headAir), nuHead=1e7/headVac;
  let mn=Infinity; for(let J=2;J<N2SPS.Jmax;J++){ const d=F(N2SPS.Bu,N2SPS.Du,J-1)-F(N2SPS.Bl,N2SPS.Dl,J); if(d<mn) mn=d; }
  N2SPS.nu0=nuHead-mn;
})();
function n2Lines(T){
  const {Bu,Du,Bl,Dl,nu0,Jmax}=N2SPS; const F=(B,D,J)=>B*J*(J+1)-D*J*J*(J+1)*(J+1);
  const out=[];
  for(let Jl=1;Jl<=Jmax;Jl++){
    const add=(Ju,S)=>{ if(Ju<1||S<=0) return; const Eu=F(Bu,Du,Ju), nu=nu0+Eu-F(Bl,Dl,Jl); const lv=1e7/nu; out.push({l:vac2air(lv), I:Math.pow(nu,4)*S*Math.exp(-C2*Eu/T)}); };
    add(Jl+1,(Jl+2)*Jl/(Jl+1));            // R
    add(Jl,(2*Jl+1)/(Jl*(Jl+1)));          // Q
    if(Jl>=2) add(Jl-1,(Jl+1)*(Jl-1)/Jl); // P
  }
  return out;
}
function n2Sim(x,T,fwhm,shift){
  const L=n2Lines(T), s=fwhm/(2*Math.sqrt(2*Math.LN2)), out=new Float64Array(x.length);
  const cut=4*fwhm;
  for(const ln of L){ const c=ln.l+shift; const i0=idxAt(x,c-cut), i1=idxAt(x,c+cut);
    for(let i=i0;i<=i1;i++){ const d=(x[i]-c)/s; out[i]+=ln.I*Math.exp(-0.5*d*d); } }
  let mx=0; for(const v of out) mx=Math.max(mx,v); if(mx>0) for(let i=0;i<out.length;i++) out[i]/=mx;
  return out;
}
// fit one spectrum in window [lo,hi]; fixedW: number or null
function fitN2(x,y,lo,hi,init,fixedW){
  const i0=idxAt(x,lo), i1=idxAt(x,hi); const xs=Array.from(x.slice(i0,i1+1)), ys=Array.from(y.slice(i0,i1+1));
  if(xs.length<12) return null;
  const xc=(lo+hi)/2, one=xs.map(()=>1), lin=xs.map(v=>v-xc);
  const model=(T,w,sh)=>{ const sim=Array.from(n2Sim(xs,T,w,sh)); const c=lstsq([sim,one,lin],ys); if(!c) return null; const fit=sim.map((v,i)=>c[0]*v+c[1]+c[2]*lin[i]); return {fit,c,sim}; };
  const chi=(T,w,sh)=>{ const m=model(T,w,sh); if(!m||m.c[0]<=0) return 1e300; let s=0; for(let i=0;i<ys.length;i++){ const r=ys[i]-m.fit[i]; s+=r*r; } return s; };
  const T0=init.T||400, W0=fixedW||init.w||0.15;
  const tr=p=>({T:Math.exp(p[0]), w:fixedW? fixedW : Math.exp(p[1]), sh:fixedW? p[1] : p[2]});
  const f=p=>{ const q=tr(p); if(q.T<120||q.T>8000||q.w<0.01||q.w>3||Math.abs(q.sh)>0.6) return 1e300; return chi(q.T,q.w,q.sh); };
  const x0= fixedW? [Math.log(T0),0] : [Math.log(T0),Math.log(W0),0];
  const st= fixedW? [0.4,0.05] : [0.4,0.3,0.05];
  let r=nelderMead(f,x0,st,{maxIt:500}); r=nelderMead(f,r.x,st.map(v=>v/4),{maxIt:400});
  const q=tr(r.x), m=model(q.T,q.w,q.sh);
  const dof=Math.max(1,ys.length-(fixedW?5:6)), s2=r.f/dof;
  // χ² profile in T (other parameters fixed at best) → 1σ where Δχ² = s2
  const prof=[]; const Ts=[]; const tmin=Math.max(130,q.T*0.4), tmax=q.T*2.2; for(let k=0;k<40;k++) Ts.push(tmin*Math.pow(tmax/tmin,k/39));
  for(const T of Ts) prof.push({T, chi2:chi(T,q.w,q.sh)/s2});
  const cmin=r.f/s2;
  // 1σ from Δχ² = 1 (χ² scaled by reduced χ²): walk outwards from the optimum on each side
  const cross=dir=>{ let prevT=q.T, prevD=0; for(let k=1;k<=120;k++){ const T=q.T*Math.pow(1+dir*0.004,k); const d=chi(T,q.w,q.sh)/s2-cmin; if(d>=1){ const t=(1-prevD)/(d-prevD); return Math.abs(prevT+(T-prevT)*t-q.T); } prevT=T; prevD=d; } return NaN; };
  const sT=0.5*(cross(1)+cross(-1));
  const ssr=r.f, sst=ys.reduce((a,v)=>a+(v-mean(ys))**2,0);
  return { T:q.T, w:q.w, shift:q.sh, xs, ys, fit:m.fit, scale:m.c[0], r2:1-ssr/sst, sT, profile:prof, cmin };
}
// Gaussian + linear baseline fit to one line → FWHM (instrument width estimate)
function fitGaussLine(x,y,l0,half=0.9){
  const i0=idxAt(x,l0-half), i1=idxAt(x,l0+half); const xs=Array.from(x.slice(i0,i1+1)), ys=Array.from(y.slice(i0,i1+1)); if(xs.length<6) return null;
  const one=xs.map(()=>1);
  const g=(mu,w)=>xs.map(v=>Math.exp(-4*Math.LN2*((v-mu)/w)**2));
  const f=p=>{ const w=Math.exp(p[1]); if(w<0.005||w>3) return 1e300; const G=g(p[0],w), c=lstsq([G,one],ys); if(!c||c[0]<=0) return 1e300; let s=0; for(let i=0;i<ys.length;i++){ const r=ys[i]-c[0]*G[i]-c[1]; s+=r*r; } return s; };
  let ib=0; for(let i=1;i<ys.length;i++) if(ys[i]>ys[ib]) ib=i;
  const r=nelderMead(f,[xs[ib],Math.log(0.2)],[0.05,0.3],{maxIt:400});
  return {mu:r.x[0], fwhm:Math.exp(r.x[1])};
}

// integrated band / line intensity in a window with linear baseline (for ratios)
const RATIO_FEATURES = [
  {key:'OH',  label:'OH A–X 306–312 nm', lo:306.0, hi:312.6},
  {key:'N2',  label:'N₂ SPS (0,0) 337 nm', lo:334.6, hi:337.6},
  {key:'N2p', label:'N₂⁺ FNS (0,0) 391 nm', lo:389.6, hi:391.8},
  {key:'Hb',  label:'Hβ 486.1 nm', lo:485.5, hi:486.8},
  {key:'Ha',  label:'Hα 656.3 nm', lo:655.6, hi:657.0},
  {key:'O',   label:'O I 777 nm', lo:776.7, hi:778.2},
  {key:'Ar',  label:'Ar I 750.4 nm (reference)', lo:749.95, hi:750.85}
];
function windowArea(x,y,lo,hi){ const c=(lo+hi)/2, hw=(hi-lo)/2; return lineArea(x,y,c,hw); }

if(typeof module!=='undefined') module.exports={mean,sdev,median,interp,idxAt,nAir,vac2air,linfit,lstsq,nelderMead,parseSpectrum,AR_BOLTZ_DEFAULT,ID_LINES,processCondition,autoShift,findPeaks,lineArea,boltzmann,N2SPS,n2Lines,n2Sim,fitN2,fitGaussLine,RATIO_FEATURES,windowArea,KB_EV};

// ===================== LTE consistency checks =====================
// McWhirter criterion: ne ≥ 1.6e12 · Te^0.5 · ΔE^3  (cm^-3; Te in K, ΔE in eV)
const sciU = v => { if(!(v>0)) return '—'; const ex=Math.floor(Math.log10(v)), m=v/Math.pow(10,ex); const sup=String(ex).split('').map(c=>'⁰¹²³⁴⁵⁶⁷⁸⁹'['0123456789'.indexOf(c)]||'⁻').join(''); return m.toFixed(1)+'×10'+sup; };
const mcwhirter = (TeK, dE) => 1.6e12*Math.sqrt(TeK)*dE*dE*dE;
function lteChecks(res, opts){
  const out={ reasons:[], notes:[] };
  if(!res.fit){ out.verdict='insufficient'; return out; }
  const use=res.use;
  out.r2=res.fit.r2; out.chi2red=res.fit.chi2red; out.T=res.T;
  const g4=use.filter(r=>r.E<14), g5=use.filter(r=>r.E>=14);
  const fitT=rows=>{ if(rows.length<3) return null; const Es=rows.map(r=>r.E); if(Math.max(...Es)-Math.min(...Es)<0.12) return null; const f=linfit(Es, rows.map(r=>r.y)); return f && f.b<0 ? {T:-1/(KB_EV*f.b), sT:f.sb/(KB_EV*f.b*f.b), r2:f.r2} : {T:NaN}; };
  out.T4=fitT(g4); out.T5=fitT(g5);
  if(g4.length && g5.length){ const mE=rs=>mean(rs.map(r=>r.E)), mY=rs=>mean(rs.map(r=>r.y)); const b=(mY(g5)-mY(g4))/(mE(g5)-mE(g4)); out.T45 = b<0 ? -1/(KB_EV*b) : NaN; }
  // lines from the same upper level must agree after dividing by gA (branching-ratio test)
  const groups={}; for(const r of use){ (groups[r.up+'@'+r.E.toFixed(3)] ||= []).push(r); }
  out.sameUpper=Object.values(groups).filter(g=>g.length>1).map(g=>{ const ys=g.map(r=>r.y); return { up:g[0].up, E:g[0].E, lines:g.map(r=>r.l), dy:Math.max(...ys)-Math.min(...ys), low:g.map(r=>r.low) }; });
  // flags
  if(res.Espan<0.5) out.reasons.push(`The upper levels span only ${res.Espan.toFixed(2)} eV, so the slope (and T) is poorly constrained.${g5.length?'':' Add lines from higher levels (5p→4s at 404–435 nm, or 5d/6d→4p at 545–703 nm).'}`);
  if(out.r2<0.9) out.reasons.push(`The points scatter about the line (R² = ${out.r2.toFixed(2)}), so the populations do not follow a single Boltzmann distribution.`);
  const bad=out.sameUpper.filter(s=>s.dy>0.3);
  if(bad.length) out.reasons.push(`Lines from the same upper level disagree by up to ${Math.max(...bad.map(b=>b.dy)).toFixed(2)} in ln units (${bad.map(b=>b.lines.map(l=>l.toFixed(1)).join('/')).join('; ')}). That points to self-absorption of lines ending on metastable levels, or to a calibration error, not to temperature.`);
  if(out.T4 && out.T45 && Number.isFinite(out.T4.T) && Number.isFinite(out.T45) && Math.abs(out.T4.T/out.T45-1)>0.3) out.reasons.push(`The 4p levels alone give ${Math.round(out.T4.T)} K, while the step from 4p to the higher levels (5p, 5d, 6d; E ≥ 14 eV) gives ${Math.round(out.T45)} K. Different groups of levels give different temperatures, which is a typical non-LTE signature.`);
  if(out.T4 && !Number.isFinite(out.T4.T)) out.reasons.push('The 4p levels alone give a positive slope (no physical temperature).');
  // McWhirter
  const Te = Number.isFinite(res.T)? res.T : 1e4;
  out.mcwFull = mcwhirter(Te, 11.83);      // ground ↔ first resonance level (1s4)
  out.mcwPartial = mcwhirter(Te, 1.6);      // 4s ↔ 4p gap (partial LTE of excited levels)
  if(opts.ne>0){
    out.ne=opts.ne;
    if(opts.ne<out.mcwPartial) out.reasons.push(`nₑ = ${sciU(opts.ne)} cm⁻³ is below the McWhirter limit even for partial LTE among the excited levels (${sciU(out.mcwPartial)} cm⁻³). Collisions cannot hold the 2p levels in Boltzmann balance, so radiative and metastable processes dominate.`);
    else if(opts.ne<out.mcwFull) out.notes.push(`nₑ = ${sciU(opts.ne)} cm⁻³ satisfies partial LTE for the excited levels (≥ ${sciU(out.mcwPartial)}) but not full LTE (≥ ${sciU(out.mcwFull)}). T_exc may describe the excited levels, but it need not equal Tₑ.`);
  } else out.notes.push(`Full LTE in argon needs nₑ ≳ ${sciU(out.mcwFull)} cm⁻³ and partial LTE of the excited levels needs ≳ ${sciU(out.mcwPartial)} cm⁻³ (McWhirter, at T = ${Math.round(Te)} K). Enter nₑ from the Ar 696.5 nm broadening tab to test this.`);
  out.verdict = out.reasons.length ? 'nonLTE' : 'consistent';
  return out;
}

// ===================== Line broadening: Ar I 696.5 nm =====================
// pseudo-Voigt (Thompson–Cox–Hastings), area-normalised, FWHMs fG, fL
function pvoigt(x,mu,fG,fL){ fG=Math.max(fG,1e-6); fL=Math.max(fL,0);
  const f=Math.pow(fG**5+2.69269*fG**4*fL+2.42843*fG**3*fL**2+4.47163*fG**2*fL**3+0.07842*fG*fL**4+fL**5,0.2);
  const q=fL/f, eta=1.36603*q-0.47719*q*q+0.11116*q*q*q;
  const d=x-mu, G=Math.sqrt(4*Math.LN2/Math.PI)/f*Math.exp(-4*Math.LN2*d*d/(f*f)), L=(2/(Math.PI*f))/(1+4*d*d/(f*f));
  return eta*L+(1-eta)*G; }
const voigtFWHM=(fG,fL)=>0.5346*fL+Math.sqrt(0.2166*fL*fL+fG*fG);
// fit Lorentzian width with Gaussian width fixed (wG) or free (wG=null)
function fitVoigtLine(x,y,l0,half,wG){
  const i0=idxAt(x,l0-half), i1=idxAt(x,l0+half); const xs=Array.from(x.slice(i0,i1+1)), ys=Array.from(y.slice(i0,i1+1)); if(xs.length<8) return null;
  const one=xs.map(()=>1), lin=xs.map(v=>v-l0);
  let ib=0; for(let i=1;i<ys.length;i++) if(ys[i]>ys[ib]) ib=i;
  const shape=(mu,g,l)=>xs.map(v=>pvoigt(v,mu,g,l));
  const ev=p=>{ const mu=p[0], l=Math.exp(p[1]), g= wG ? wG : Math.exp(p[2]); if(l<1e-5||l>2||g<1e-4||g>2) return {s:1e300}; const S=shape(mu,g,l), c=lstsq([S,one,lin],ys); if(!c||c[0]<=0) return {s:1e300};
    let s=0; const fit=S.map((v,i)=>c[0]*v+c[1]+c[2]*lin[i]); for(let i=0;i<ys.length;i++){ const r=ys[i]-fit[i]; s+=r*r; } return {s,fit,c,g,l,mu}; };
  const f=p=>ev(p).s;
  const p0= wG ? [xs[ib],Math.log(0.02)] : [xs[ib],Math.log(0.02),Math.log(0.03)];
  const st= wG ? [0.005,0.6] : [0.005,0.6,0.4];
  let r=nelderMead(f,p0,st,{maxIt:600}); r=nelderMead(f,r.x,st.map(v=>v/3),{maxIt:600});
  const e=ev(r.x); const sst=ys.reduce((a,v)=>a+(v-mean(ys))**2,0);
  return { mu:e.mu, wL:e.l, wG:e.g, wV:voigtFWHM(e.g,e.l), xs, ys, fit:e.fit, r2:1-e.s/sst };
}
// Coulomb-approximation mean-square radius (units a0²) for a level with binding energy (Eion−E) eV and orbital l
const r2coul=(Eion,E,l)=>{ const n2=13.6057/(Eion-E); return n2/2*(5*n2+1-3*l*(l+1)); };
const AR696 = { l:696.543, M:39.948, Eion:15.7596, Eup:13.328, lup:1, Elow:11.548, llow:0 };
function vdwC(lam, R2){ return 7.5e-7*lam*lam*Math.pow(R2,0.4); } // nm·K^0.7, Ar perturbers at 1 atm (Yubero et al.)

if(typeof module!=='undefined') Object.assign(module.exports,{lteChecks,mcwhirter,pvoigt,voigtFWHM,fitVoigtLine,r2coul,AR696,vdwC});

// ===================== Synthetic example data (clearly labelled in the UI) =====================
function rng(seed){ return function(){ seed|=0; seed=seed+0x6D2B79F5|0; let t=Math.imul(seed^seed>>>15,1|seed); t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
function gaussRand(r){ let u=0,v=0; while(!u) u=r(); while(!v) v=r(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); }
const STARK696 = { wS:0.0112 };  // nm FWHM per 1e16 cm^-3 at Te ≈ 1e4 K (default, editable)
// Ar I 696.543 nm Stark FWHM models. ne in cm^-3, Te in eV, W in nm.
// kunze : Kunze, Introduction to Plasma Spectroscopy, eq. 10.26  (W = 8.50e-25 nm · Te^-0.37 · ne[m^-3])
// lu    : linear, W = wS · ne/1e16  (wS = 0.0112 nm reproduces Lu et al., Sci. Rep. 4, 7488, 2014)
// jap   : fit used in Ananthanarasimhan, Gangwar et al., J. Appl. Phys. 129, 223301 (2021): ne[m^-3] = 0.99·10^23.83·W  (Griem formula, Te 1.7–2.4 eV)
// griem : Griem/Kunze 9.57 & 9.62 with user w (electron half-halfwidth, nm, at 1e16 cm^-3) and A (ion parameter at 1e16 cm^-3)
const STARK_MODELS = {
  kunze:{ name:'Kunze (10.26), Tₑ-dependent', W:(ne,o)=>8.50e-25*Math.pow(o.Te,-0.37)*ne*1e6 },
  lu:   { name:'Linear w<sub>S</sub> (Lu et al. 2014 default)', W:(ne,o)=>o.wS*ne/1e16 },
  jap:  { name:'Griem formula, JAP 2021 fit (Tₑ 1.7–2.4 eV)', W:(ne)=>ne*1e6/(0.99*Math.pow(10,23.83)) },
  griem:{ name:'Griem formula with w, A from tables', W:(ne,o)=>{ const x=ne/1e16, A=o.gA*Math.pow(x,0.25), R=8.34e-5*Math.pow(o.Te,-0.5)*Math.pow(ne*1e6,1/6); return 2*o.gw*x*(1+1.75*A*(1-0.75*R)); } }
};
function starkW(model,ne,o){ return STARK_MODELS[model].W(ne,o); }
function starkNe(model,W,o){ if(!(W>0)) return 0; // invert monotonic W(ne) by bisection in log space
  let lo=1e10, hi=1e20; for(let k=0;k<100;k++){ const mid=Math.sqrt(lo*hi); if(starkW(model,mid,o)>W) hi=mid; else lo=mid; } return Math.sqrt(lo*hi); }
// van der Waals FWHM (nm) of Ar 696.5 nm in Ar: Lee et al. 1975 / Konjević: 8.18e-12 λ²(αR²)^0.4 (Tg/μ)^0.3 N0  (cgs; λ in cm, α cm³, R² in a0², N0 cm^-3)
function vdwLee(lamNm,R2a0,Tg,N0cm3,alpha=16.54e-25,mu=19.97){ const lam=lamNm*1e-7; return 8.18e-12*lam*lam*Math.pow(alpha*R2a0,0.4)*Math.pow(Tg/mu,0.3)*N0cm3*1e7; }
function makeExample(variant){
  const cfg = variant===0
    ? { name:'Example A — 20 W, 2 slm (synthetic)', T4:8200, T5:8200, Trot:360, Tg:360, sa:{763.511:0.82, 801.479:0.9, 811.531:0.75}, oh:1.0, n2:0.55, n2p:0.08, ha:0.35, hb:0.06, o:0.12, ne:4e15, seed:11 }
    : { name:'Example B — 35 W, 2 slm (synthetic)', T4:4200, T5:9000, Trot:430, Tg:430, sa:{763.511:0.7, 801.479:0.85, 696.543:0.85, 811.531:0.6}, oh:1.4, n2:0.8, n2p:0.15, ha:0.5, hb:0.09, o:0.2, ne:9e15, seed:29 };
  const r=rng(cfg.seed);
  const x=Float64Array.from({length:Math.round((810-200)/0.04)+1},(_,i)=>200+i*0.04);
  const W=0.15, sg=W/(2*Math.sqrt(2*Math.LN2));
  const prof=new Float64Array(x.length);
  const addLine=(l,area,w=W)=>{ const s=w/(2*Math.sqrt(2*Math.LN2)); const i0=idxAt(x,l-4*w), i1=idxAt(x,l+4*w); for(let i=i0;i<=i1;i++){ const d=(x[i]-l)/s; prof[i]+=area/(s*Math.sqrt(2*Math.PI))*Math.exp(-0.5*d*d); } };
  // Ar I lines: energy units, I ∝ gA/λ · exp(−E/kT); 4p and 5p groups may follow different T (non-LTE example)
  for(const d of AR_BOLTZ_DEFAULT){
    const T = d.E<14 ? cfg.T4 : cfg.T5;
    let I = d.g*d.A/d.l*Math.exp(-(d.E-13.0)/(KB_EV*T));
    I *= (cfg.sa[d.l]||1);
    addLine(d.l, I*1.2e-2);
  }
  // molecular bands from the simplified N2 model
  const band=(T,head,scale)=>{ const L=n2Lines(T); let s=0; for(const q of L) s+=q.I; const sh=head-337.13; for(const q of L) addLine(q.l+sh, scale*q.I/s); };
  band(cfg.Trot,337.13,5000*cfg.n2); band(cfg.Trot,357.69,6000*cfg.n2); band(cfg.Trot,380.49,2600*cfg.n2); band(cfg.Trot,315.93,3000*cfg.n2);
  band(cfg.Trot,333.89,900*cfg.n2); band(cfg.Trot,353.67,1500*cfg.n2); band(cfg.Trot,375.54,1500*cfg.n2); band(cfg.Trot,371.05,600*cfg.n2); band(cfg.Trot,313.60,700*cfg.n2);
  band(cfg.Trot*1.1,391.44,9000*cfg.n2p);
  // OH A–X (0,0): illustrative envelope only (not a physical model)
  for(const [l,a] of [[306.36,0.25],[306.7,0.35],[307.0,0.4],[307.3,0.45],[307.6,0.5],[307.84,0.55],[308.15,0.7],[308.33,0.8],[308.52,0.9],[308.9,1.6],[309.05,1.3],[309.24,1.0],[309.53,0.7],[309.9,0.5],[310.2,0.4],[310.6,0.32],[311.0,0.25],[311.5,0.18],[312.2,0.25],[312.6,0.15],[313.2,0.1]]) addLine(l, 380*cfg.oh*a);
  addLine(656.279,5200*cfg.ha,0.17); addLine(486.133,5200*cfg.hb,0.17); addLine(434.047,5200*cfg.hb*0.4,0.17);
  addLine(777.194,4000*cfg.o*0.42); addLine(777.417,4000*cfg.o*0.33); addLine(777.539,4000*cfg.o*0.25);
  // weak continuum
  for(let i=0;i<x.length;i++) prof[i]+=4+3*Math.exp(-Math.pow((x[i]-450)/180,2));
  const ambient=new Float64Array(x.length); // room-light lines appear in the background too
  const addAmb=(l,a)=>{ const i0=idxAt(x,l-0.6), i1=idxAt(x,l+0.6); for(let i=i0;i<=i1;i++){ const d=(x[i]-l)/sg; ambient[i]+=a/(sg*Math.sqrt(2*Math.PI))*Math.exp(-0.5*d*d); } };
  addAmb(435.833,25); addAmb(546.074,40); addAmb(611.6,10);
  const plasma=[], bg=[];
  for(let k=0;k<3;k++){ const amp=1+0.025*gaussRand(r); const y=Float64Array.from(prof,(v,i)=>Math.max(0,amp*v+ambient[i]+12+3*gaussRand(r)+Math.sqrt(Math.max(v,0))*0.6*gaussRand(r))); plasma.push({name:`exampleA_plasma_${k+1}.txt`.replace('A',variant?'B':'A'), x, y}); }
  for(let k=0;k<3;k++){ const y=Float64Array.from(ambient,(v,i)=>Math.max(0,v+12+3*gaussRand(r))); bg.push({name:`example${variant?'B':'A'}_background_${k+1}.txt`, x, y}); }
  // high-resolution Ar I 696.5 nm record
  const xh=Float64Array.from({length:501},(_,i)=>695.8+i*0.003);
  const wI=0.025, wD=7.16e-7*696.543*Math.sqrt(cfg.Tg/39.948), wG=Math.hypot(wI,wD);
  const R2=r2coul(AR696.Eion,AR696.Eup,1)-r2coul(AR696.Eion,AR696.Elow,0), wL=vdwC(696.543,R2)/Math.pow(cfg.Tg,0.7)+STARK696.wS*cfg.ne/1e16;
  const hp=[], hb=[];
  for(let k=0;k<3;k++){ const amp=1+0.02*gaussRand(r); hp.push({name:`example${variant?'B':'A'}_Ar696_hires_${k+1}.txt`, x:xh, y:Float64Array.from(xh,v=>Math.max(0,amp*600*pvoigt(v,696.543+0.004,wG,wL)+15+2.5*gaussRand(r)))}); }
  for(let k=0;k<3;k++) hb.push({name:`example${variant?'B':'A'}_Ar696_bg_${k+1}.txt`, x:xh, y:Float64Array.from(xh,()=>15+2.5*gaussRand(r))});
  return { name:cfg.name, example:true, plasma, bg, expP:null, expB:null, hr:{plasma:hp, bg:hb, wI}, truth:cfg };
}
if(typeof module!=='undefined') Object.assign(module.exports,{makeExample,STARK696,STARK_MODELS,starkW,starkNe,vdwLee});
