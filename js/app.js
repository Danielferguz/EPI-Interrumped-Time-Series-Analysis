// EPI-ITS · navegación, datos, interrupciones, regresión segmentada, pronóstico, grupos y reporte.
import { $, $$, esc, f2, fp, opt, errBox, howto, card, table, download } from './state.js?v=2026.10.25-1';
import { readFile, toCSV } from './io.js?v=2026.10.25-1';
import { bindFigs, figWrap } from './figs.js?v=2026.10.25-1';
import { fitITS, fitITSIndividual, aggregateByPeriod, poolEffects, nwLag } from './its.js?v=2026.10.25-1';
import { METHODS, METHOD_EN, counterfactualEffects, rollingOrigin, fitMethod } from './forecast.js?v=2026.10.25-1';
import { itsPlot, panelPlot, diffPlot, acfPlot, forestPlot, accuracyPlot, SEG_COLORS } from './itsplots.js?v=2026.10.25-1';
import { timeAxis, makeTicks, UNITS, PERIOD } from './timeutil.js?v=2026.10.25-1';
import { EXAMPLES } from './examples.js?v=2026.10.25-1';
import { GUIDE_HTML, HELP_HTML } from './help.js?v=2026.10.25-1';
import { renderReport } from './report.js?v=2026.10.25-1';

export const S = { raw: [], cols: [], name: '', dtype: 'agg', map: {}, ax: null, tk: null, brk: [], lag: 0, opt: {}, fit: null, cache: {}, fc: null };
const num = v => { const x = typeof v === 'number' ? v : Number(String(v).replace(',', '.')); return Number.isFinite(x) && String(v).trim() !== '' ? x : NaN; };
const yl = () => S.opt.ylab || S.map.y;
const dash = v => (Number.isFinite(v) ? String(f2(v, Math.abs(v) < 10 ? 2 : 1)).replace('-', '−') : '–');
const pc = v => (Number.isFinite(v) ? `${(100 * v).toFixed(Math.abs(v) < 0.1 ? 1 : 0).replace('-', '−')}%` : '–');
const ci = (e, lo, hi, d) => `${d(e)} (${d(lo)}, ${d(hi)})`;

// ---------------------------------------------------------------- navegación
function showTab(t) {
  $$('.tabs button').forEach(b => b.classList.toggle('active', b.dataset.tab === t)); $$('.tab').forEach(s => s.classList.toggle('active', s.id === 'tab-' + t));
  if (t === 'diseno') renderDesign(); if (t === 'pronostico') renderFcForm(); if (t === 'grupos') renderGroups(); if (t === 'reporte') renderReport();
  window.scrollTo({ top: 0 });
}
export { showTab };
$$('.tabs button').forEach(b => b.addEventListener('click', () => showTab(b.dataset.tab)));
document.addEventListener('click', e => { const g = e.target.closest('[data-go]'); if (g) showTab(g.dataset.go); });
const gate = () => { $$('.need-data').forEach(e => (e.hidden = !!S.raw.length)); $$('.has-data').forEach(e => (e.hidden = !S.raw.length)); };
$('#guide').innerHTML = GUIDE_HTML; $('#help').innerHTML = HELP_HTML;

// ---------------------------------------------------------------- ejemplos y carga
$('#examples').innerHTML = EXAMPLES.map((x, k) => `<button class="ghost" data-ex="${k}"><em>${esc(x.tag)}</em><b>${esc(x.title)}</b><small>${esc(x.desc)}</small></button>`).join('');
$$('[data-ex]').forEach(b => b.addEventListener('click', async () => {
  const x = EXAMPLES[+b.dataset.ex]; b.disabled = true;
  try { const txt = await (await fetch('examples/' + x.file)).text(); await loadFile(new File([txt], x.file), x); showTab(x.go || 'resultados'); if ((x.go || 'resultados') === 'resultados') runFit(); } catch (e) { alert(e.message); }
  b.disabled = false;
}));
const drop = $('#drop'); ['dragover', 'dragenter'].forEach(ev => drop.addEventListener(ev, e => { e.preventDefault(); drop.classList.add('over'); })); ['dragleave', 'drop'].forEach(ev => drop.addEventListener(ev, () => drop.classList.remove('over')));
drop.addEventListener('drop', e => { e.preventDefault(); const f = e.dataTransfer.files[0]; if (f) loadFile(f).catch(err => ($('#load-msg').innerHTML = errBox(err))); });
$('#file').addEventListener('change', e => e.target.files[0] && loadFile(e.target.files[0]).catch(err => ($('#load-msg').innerHTML = errBox(err))));

async function loadFile(file, ex = null) {
  const { sheets } = await readFile(file), name = Object.keys(sheets)[0], raw = sheets[name].map(r => ({ ...r }));
  if (raw.length < 6) throw new Error('El archivo necesita al menos 6 filas.');
  Object.assign(S, { raw, cols: Object.keys(raw[0]), name: file.name.replace(/\.[^.]+$/, ''), fit: null, cache: {}, fc: null, ex, cur: undefined, gr: null });
  const g = guess(ex); S.dtype = g.dtype; S.map = g.map; S.opt = { family: 'gaussian', ac: 'nw', nwL: 0, K: 0, P: 0, display: 'trend', scale: 1, zero: false, level: 0.95, unit: 'auto', ...(ex?.opt || {}) };
  if (ex?.brk) { S.brk = []; S.lag = ex.lag || 0; }
  gate(); drawMap(); $('#load-msg').innerHTML = `<span class="pill">✓ ${esc(file.name)}: ${raw.length} filas × ${S.cols.length} columnas</span>`;
  applyTime(); if (ex?.brk) S.brk = ex.brk.map(b => ({ T: b.T, name: b.name, eff: b.eff || 'both' })); if (ex?.opt) Object.assign(S.opt, ex.opt); $('#dz-lag').value = S.lag; renderDesign();
}
function guess(ex) {
  if (ex) return { dtype: ex.dtype, map: ex.map };
  const cols = S.cols, isNum = c => S.raw.slice(0, 60).filter(r => Number.isFinite(num(r[c]))).length > 40, tCol = cols.find(c => /date|month|quarter|time|period|year|fecha|mes|trimestre|a[ñn]o/i.test(c)) || cols[0];
  const idCol = cols.find(c => /(^|_)(id|patient|paciente|subject)/i.test(c) && c !== tCol), grp = cols.find(c => c !== tCol && !isNum(c) && new Set(S.raw.map(r => r[c])).size <= 40 && c !== idCol);
  const ycol = cols.find(c => c !== tCol && c !== idCol && isNum(c)) || cols[1], ind = idCol && new Set(S.raw.map(r => r[idCol])).size < S.raw.length * 0.8;
  return { dtype: ind ? 'ind' : 'agg', map: { time: tCol, y: ycol, n: '', group: grp || '', id: ind ? idCol : '' } };
}
$$('[data-dt]').forEach(b => b.addEventListener('click', () => { S.dtype = b.dataset.dt; S.fit = null; S.cache = {}; drawMap(); }));
function drawMap() {
  $$('[data-dt]').forEach(b => b.classList.toggle('active', b.dataset.dt === S.dtype));
  const sel = (id, label, v, blank) => `<label>${label}<select id="${id}">${blank ? '<option value="">(ninguna)</option>' : ''}${S.cols.map(c => opt(c, c, c === v)).join('')}</select></label>`, m = S.map, ind = S.dtype === 'ind';
  $('#map-form').innerHTML = sel('m-time', 'Tiempo (fecha, periodo o número)', m.time) + sel('m-y', ind ? 'Resultado de cada medición' : 'Resultado (tasa, %, conteo, media…)', m.y) +
    (ind ? sel('m-id', 'Identificador de la persona', m.id) : sel('m-n', 'Denominador o población (opcional)', m.n, true)) + sel('m-g', ind ? 'Grupo (opcional, p. ej., centro)' : 'Grupo o panel (opcional, p. ej., región)', m.group, true) +
    `<label>Unidad de tiempo<select id="m-unit">${UNITS.map(([v, t]) => opt(v, t, v === S.opt.unit)).join('')}</select></label>`;
  ['m-time', 'm-y', 'm-n', 'm-id', 'm-g', 'm-unit'].forEach(id => $('#' + id)?.addEventListener('change', readMap)); readMap(true);
}
function readMap(quiet) {
  S.map = { time: $('#m-time').value, y: $('#m-y').value, n: $('#m-n')?.value || '', id: $('#m-id')?.value || '', group: $('#m-g').value }; S.opt.unit = $('#m-unit').value; S.fit = null; S.cache = {}; S.fc = null;
  try { applyTime(); $('#map-note').textContent = ''; } catch (e) { $('#map-note').innerHTML = errBox(e); return; }
  const ax = S.ax, groups = S.map.group ? [...new Set(S.raw.map(r => r[S.map.group]))].length : 0;
  $('#map-note').innerHTML = `Unidad detectada: <b>${esc(UNITS.find(u => u[0] === ax.unit)?.[1] || ax.unit)}</b> · ${ax.N} periodos (${esc(ax.label(0))} a ${esc(ax.label(ax.N - 1))})${groups ? ` · ${groups} grupos` : ''}${S.dtype === 'ind' ? ` · ${new Set(S.raw.map(r => r[S.map.id])).size} personas` : ''}.`;
  if (!S.opt.P && S.opt.K) S.opt.P = PERIOD[ax.unit] || 0;
  const t = S.raw.slice(0, 8), cols = S.cols; $('#d-prev').innerHTML = `<table class="data"><thead><tr>${cols.map(c => `<th>${esc(c)}</th>`).join('')}</tr></thead><tbody>${t.map(r => `<tr>${cols.map(c => `<td>${esc(r[c])}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
  quickPlot(); if (quiet !== true) renderDesign();
}
function applyTime() { S.ax = timeAxis(S.raw.map(r => r[S.map.time]), S.opt.unit || 'auto'); S.tk = makeTicks(S.ax); }
function quickPlot() {
  try { const sr = series(S.map.group ? groupNames()[0] : null); $('#d-plot').innerHTML = figWrap(itsPlot({ title: S.map.group ? `Vista previa · ${groupNames()[0]}` : 'Vista previa de la serie', ylab: S.map.y, obs: sr.obs, ticks: S.tk.ticks, minor: S.tk.minor, breaks: [], segs: [], cfs: [], segNames: ['Observed'], zero: S.opt.zero, H: 340, xlab: '' }), 'vista_previa'); bindFigs($('#d-plot')); } catch { $('#d-plot').innerHTML = ''; }
}

// ---------------------------------------------------------------- series por grupo
export const groupNames = () => (S.map.group ? [...new Set(S.raw.map(r => String(r[S.map.group])))].sort((a, b) => a.localeCompare(b, 'es', { numeric: true })) : []);
const scaleOf = () => (S.opt.scale > 0 ? S.opt.scale : 1);
// serie agregada de un grupo (null = todo): y, n y valor observado en la escala de la figura
export function series(g) {
  const N = S.ax.N, y = Array(N).fill(NaN), n = Array(N).fill(NaN), cnt = Array(N).fill(0), fam = S.opt.family, ind = S.dtype === 'ind', all = g == null || g === '__all';
  const rows = S.raw.map((r, i) => [r, S.ax.idx[i]]).filter(([r]) => all || String(r[S.map.group]) === g);
  if (ind) { const obs = rows.map(([r, t]) => ({ t, y: num(r[S.map.y]), id: r[S.map.id] })).filter(r => Number.isFinite(r.y)); const a = aggregateByPeriod(obs, N, fam, scaleOf()); return { y: a.y, n: a.n, obs: a.y, se: a.se, N, rows: obs }; }
  for (const [r, t] of rows) { const v = num(r[S.map.y]), d = S.map.n ? num(r[S.map.n]) : NaN; if (!Number.isFinite(v)) continue; if (fam === 'gaussian') { const w = Number.isFinite(d) && d > 0 ? d : 1; y[t] = (Number.isFinite(y[t]) ? y[t] * cnt[t] + v * w : v * w) / (cnt[t] + w); cnt[t] += w; if (Number.isFinite(d)) n[t] = (n[t] || 0) + d; } else { y[t] = (Number.isFinite(y[t]) ? y[t] : 0) + v; n[t] = (Number.isFinite(n[t]) ? n[t] : 0) + (Number.isFinite(d) ? d : 0); } }
  const sc = scaleOf(), obs = y.map((v, t) => (fam === 'gaussian' ? v : Number.isFinite(v) && n[t] > 0 ? v / n[t] * sc : NaN));
  return { y, n: n.map(v => (Number.isFinite(v) && v > 0 ? v : NaN)), obs, N };
}
const brkIdx = () => S.brk.map(b => b.T), effOf = () => S.brk.map(b => b.eff || 'both');
function fitOpts() {
  const o = S.opt, P = o.K > 0 ? (o.P || PERIOD[S.ax.unit] || 0) : 0;
  return { breaks: brkIdx(), effects: effOf(), lag: S.lag, labels: S.brk.map(b => b.name), season: o.K > 0 && P > 1 ? { P, K: o.K } : null, family: o.family, ac: o.ac, nwL: o.nwL > 0 ? o.nwL : undefined, level: o.level, scale: o.family === 'gaussian' ? 1 : scaleOf(), display: o.display };
}
function runOne(g) {
  const key = String(g), cached = S.cache[key]; if (cached) return cached;
  const sr = series(g), fo = fitOpts(), fit = S.dtype === 'ind' ? fitITSIndividual(sr.rows, { ...fo, N: sr.N, family: fo.family }) : fitITS({ y: sr.y, n: sr.n }, fo);
  return (S.cache[key] = { g, sr, fit });
}

// ---------------------------------------------------------------- interrupciones y modelo
$('#brk-add').addEventListener('click', () => { if (S.brk.length >= 4) return; const N = S.ax.N, last = S.brk.length ? S.brk[S.brk.length - 1].T : Math.floor(N * 0.4), T = Math.min(N - 4, Math.max(last + 6, Math.floor(N * (0.4 + 0.15 * S.brk.length)))); S.brk.push({ T, name: `Interruption ${S.brk.length + 1}`, eff: 'both' }); renderDesign(); });
$('#dz-lag').addEventListener('input', e => { S.lag = Math.max(0, Math.floor(+e.target.value || 0)); S.cache = {}; S.fit = null; renderDesign(true); });
function renderDesign(keep) {
  if (!S.raw.length || !S.ax) return; gate();
  if (!S.brk.length) S.brk.push({ T: Math.floor(S.ax.N * 0.5), name: 'Intervention', eff: 'both' });
  const lab = S.ax.label;
  $('#brk-list').innerHTML = S.brk.map((b, j) => `<div class="brk" data-j="${j}"><span class="pill" style="background:${SEG_COLORS[(j + 1) % SEG_COLORS.length]};color:#fff">${j + 1}</span><label>Nombre (aparece en la figura)<input data-k="name" value="${esc(b.name)}"></label><label>Primer periodo posterior<select data-k="T">${Array.from({ length: S.ax.N - 3 }, (_, i) => i + 2).map(i => opt(i, lab(i), i === b.T)).join('')}</select></label><label>Efecto a estimar<select data-k="eff">${[['both', 'Nivel y pendiente'], ['level', 'Solo nivel'], ['slope', 'Solo pendiente']].map(([v, t]) => opt(v, t, v === b.eff)).join('')}</select></label>${S.brk.length > 1 ? `<button class="ghost danger" data-del="${j}" aria-label="Quitar">✕</button>` : ''}</div>`).join('');
  $$('#brk-list [data-k]').forEach(el => el.addEventListener('change', () => { const j = +el.closest('.brk').dataset.j; S.brk[j][el.dataset.k] = el.dataset.k === 'T' ? +el.value : el.value; S.brk.sort((a, b) => a.T - b.T); S.cache = {}; S.fit = null; S.fc = null; renderDesign(); }));
  $$('#brk-list [data-del]').forEach(el => el.addEventListener('click', () => { S.brk.splice(+el.dataset.del, 1); S.cache = {}; S.fit = null; S.fc = null; renderDesign(); }));
  const o = S.opt, ind = S.dtype === 'ind', fams = ind ? [['gaussian', 'Continuo (modelo mixto, intercepto aleatorio)'], ['binomial', 'Binario 0/1 (GEE logística)'], ['poisson', 'Conteo (GEE Poisson)']] : [['gaussian', 'Continuo, porcentaje o tasa (mínimos cuadrados)'], ['poisson', 'Conteo con población (Poisson, razón de tasas)'], ['binomial', 'Eventos sobre denominador (binomial, razón de odds)']];
  if (!$('#dz-form').dataset.k || !keep) $('#dz-form').innerHTML = `<label>Tipo de resultado<select id="o-family">${fams.map(([v, t]) => opt(v, t, v === o.family)).join('')}</select></label>
    ${ind ? '' : `<label>Errores autocorrelacionados<select id="o-ac">${[['nw', 'Newey–West (errores robustos)'], ['ar1', 'AR(1) de Prais–Winsten (solo continuo)'], ['none', 'Ninguno (errores clásicos)']].map(([v, t]) => opt(v, t, v === o.ac)).join('')}</select></label><label>Rezagos de Newey–West (0 = automático)<input id="o-nwL" type="number" min="0" value="${o.nwL}"></label>`}
    <label>Estacionalidad (armónicos de Fourier)<select id="o-K">${[[0, 'Ninguna'], [1, '1 armónico'], [2, '2 armónicos'], [3, '3 armónicos']].map(([v, t]) => opt(v, t, +v === o.K)).join('')}</select></label><label>Periodos por ciclo (12 mensual, 4 trimestral)<input id="o-P" type="number" min="0" value="${o.P || PERIOD[S.ax.unit] || 0}"></label>
    <label>Multiplicador para tasas / porcentajes<input id="o-scale" type="number" min="0" step="any" value="${o.scale}" title="Binomial: 100 = %; Poisson: 100000 = por 100 000"></label><label>Nivel de confianza<select id="o-level">${[0.9, 0.95, 0.99].map(v => opt(v, `${v * 100} %`, v === o.level)).join('')}</select></label>
    <label>Línea ajustada<select id="o-display">${[['trend', 'Tendencia (sin estacionalidad)'], ['full', 'Con estacionalidad']].map(([v, t]) => opt(v, t, v === o.display)).join('')}</select></label><label class="inline"><input type="checkbox" id="o-zero"${o.zero ? ' checked' : ''}> Eje Y desde 0</label>`;
  $('#dz-form').dataset.k = '1';
  $$('#dz-form select, #dz-form input').forEach(el => (el.onchange = readOpts)); explain();
}
function readOpts() {
  const o = S.opt, g = id => $('#' + id); o.family = g('o-family').value; o.ac = g('o-ac')?.value || 'nw'; o.nwL = +(g('o-nwL')?.value || 0); o.K = +g('o-K').value; o.P = +g('o-P').value; o.scale = +g('o-scale').value || 1; o.level = +g('o-level').value; o.display = g('o-display').value; o.zero = g('o-zero').checked;
  if (o.family === 'binomial' && S.dtype === 'agg' && o.scale === 1) { o.scale = 100; g('o-scale').value = 100; } if (o.family === 'gaussian') o.scale = 1;
  S.cache = {}; S.fit = null; S.fc = null; explain(); renderDesign(true);
}
function explain() {
  const o = S.opt, nb = S.brk.length; if (!$('#dz-explain')) return;
  const t = [`<b>Modelo.</b> ${nb} interrupción${nb > 1 ? 'es' : ''}: ${nb ? 'tendencia previa + ' : ''}${S.brk.map(b => `${b.eff === 'slope' ? 'cambio de pendiente' : b.eff === 'level' ? 'cambio de nivel' : 'cambio de nivel y de pendiente'} en «${esc(b.name)}»`).join(', ')}${S.lag ? `; se excluyen ${S.lag} periodo(s) de transición tras cada una` : ''}.`];
  t.push(`<b>Resultado.</b> ${{ gaussian: S.dtype === 'ind' ? 'Continuo con intercepto aleatorio por persona (la correlación entre mediciones de la misma persona).' : 'Continuo: los efectos son diferencias absolutas en las unidades del resultado.', poisson: 'Conteo: los efectos son razones de tasas (exp β), con la población como desplazamiento.', binomial: 'Binomial: los efectos son razones de odds (exp β); la figura muestra proporciones.' }[o.family]}`);
  if (S.dtype === 'agg') t.push(`<b>Autocorrelación.</b> ${{ nw: 'Newey–West corrige los errores estándar sin cambiar los coeficientes; es la opción segura por defecto.', ar1: 'Prais–Winsten modela un AR(1) en los errores y reestima los coeficientes (solo continuo).', none: 'Sin corrección: solo es razonable si el test de Ljung–Box de los residuos no muestra autocorrelación.' }[o.ac]}`);
  if (o.K > 0) t.push(`<b>Estacionalidad.</b> ${o.K} par(es) seno–coseno con ciclo de ${o.P || PERIOD[S.ax.unit] || '?'} periodos.`);
  $('#dz-explain').innerHTML = t.map(x => `<p>${x}</p>`).join('');
}
$('#dz-run').addEventListener('click', () => { showTab('resultados'); runFit(); });

// ---------------------------------------------------------------- resultados de la regresión segmentada
const eff = (e, fam) => (fam === 'gaussian' ? `${ci(e.abs.est, e.abs.lo, e.abs.hi, dash)} · ${pc(e.rel.est)}` : `${pc(e.rel.est)} (${pc(e.rel.lo)}, ${pc(e.rel.hi)})`);
function noteOf(ef, fam, fit) { const e = ef.last, text = fam === 'gaussian' ? `${dash(e.abs.est)} (${dash(e.abs.lo)}, ${dash(e.abs.hi)})` : `${pc(e.rel.est)} (${pc(e.rel.lo)}, ${pc(e.rel.hi)})`; return { j: ef.j, t: ef.t1, yo: e.obs, yc: e.cf, text }; }
export function segPlot(r, extra = {}) {
  const fit = r.fit, D = fit.D, names = ['Before ' + D.labels[0], ...D.labels.map((l, j) => `After ${l}`)];
  return itsPlot({ obs: r.sr.obs, se: S.dtype === 'ind' ? r.sr.se.map(v => v) : null, ticks: S.tk.ticks, minor: S.tk.minor, ticksFull: Array.from({ length: S.ax.N }, (_, i) => S.ax.label(i)), breaks: D.breaks.map((T, j) => ({ T, label: D.labels[j] })), lag: S.lag, segs: fit.curves.seg, cfs: fit.curves.cf, notes: fit.effects.map(e => noteOf(e, fit.family, fit)), segNames: names, ylab: yl(), zero: S.opt.zero, xlab: '', ...extra });
}
function runFit() {
  const out = $('#r-out'); out.innerHTML = '<div class="card"><p class="muted">Calculando…</p></div>';
  try {
    const gs = groupNames(), g = S.cur ?? (gs.length ? gs[0] : null); S.cur = g; const r = runOne(g); S.fit = r; renderResults(r);
    $('#r-ctl').innerHTML = `<div class="row gap wrap"><b>Resultados</b>${gs.length ? `<label class="inline">Grupo <select id="r-g"><option value="">Todos juntos</option>${gs.map(x => opt(x, x, x === g)).join('')}</select></label>` : ''}<span class="grow"></span><button class="ghost" data-go="diseno">Cambiar el modelo</button></div>`;
    $('#r-g')?.addEventListener('change', e => { S.cur = e.target.value || null; if (S.cur === null && gs.length) S.cur = '__all'; runFit(); });
  } catch (e) { out.innerHTML = card('No se pudo ajustar el modelo', errBox(e)); console.error(e); }
}
function renderResults(r) {
  const { fit, sr } = r, fam = fit.family, D = fit.D, exp = fam !== 'gaussian', nm = { gaussian: 'diferencia', poisson: 'razón de tasas', binomial: 'razón de odds' }[fam];
  const kp = [[fit.n, 'periodos/observaciones usados'], [D.p, 'parámetros'], fit.rho != null ? [f2(fit.rho, 2), 'ρ de Prais–Winsten'] : fit.L ? [fit.L, 'rezagos de Newey–West'] : null, fit.diag ? [fp(fit.diag.lb.p), `Ljung–Box (${fit.diag.lb.lag} rezagos)`] : null, fit.diag ? [f2(fit.diag.dw, 2), 'Durbin–Watson'] : null, fit.r2 != null ? [f2(fit.r2, 2), 'R²'] : null, fit.icc != null ? [f2(fit.icc, 2), 'CCI (persona)'] : null, fit.nclust ? [fit.nclust, 'personas'] : null].filter(Boolean);
  const coefRows = fit.coef.map(c => [esc(c.name), dash(c.est), `${dash(c.lo)} a ${dash(c.hi)}`, fp(c.p), ...(exp ? [c.exp ? `${f2(c.exp.est, 3)} (${f2(c.exp.lo, 3)} a ${f2(c.exp.hi, 3)})` : ''] : [])]);
  const effRows = fit.effects.flatMap(e => [['Inmediato', e.first, e.t0], [`Promedio (${e.nWin} periodos)`, e.mean, null], ['Al final del tramo', e.last, e.t1]].map(([w, c, t]) => [esc(e.label), `${w}${t != null ? ` · ${esc(S.ax.label(t))}` : ''}`, dash(c.obs), dash(c.cf), ci(c.abs.est, c.abs.lo, c.abs.hi, dash), `${pc(c.rel.est)} (${pc(c.rel.lo)}, ${pc(c.rel.hi)})`, fp(c.p)]));
  const sent = fit.effects.map(e => `<li><b>${esc(e.label)}:</b> ${sentence(e, fam, fit)}</li>`).join('');
  $('#r-out').innerHTML = card(`Regresión segmentada${r.g && r.g !== '__all' ? ` · ${esc(r.g)}` : ''}`, `${howto('Los puntos son lo observado en cada periodo; las líneas, el ajuste por tramo con su IC; la línea gris discontinua es la tendencia del tramo anterior prolongada (contrafactual).', 'Estimar cuánto cambió el nivel y la pendiente tras cada interrupción y qué diferencia hay respecto a no haber intervenido.', `La caja de cada interrupción resume la diferencia entre el ajuste y el contrafactual al final del tramo. ${exp ? 'Los efectos relativos son razones (exp β) menos 1.' : ''}`)}
    <div class="kpis">${kp.map(([v, l]) => `<div class="kpi"><b>${v}</b><span>${l}</span></div>`).join('')}</div>
    <div class="plots single">${figWrap(segPlot(r), 'serie_interrumpida')}</div>
    <h3>Cómo se lee</h3><ul>${sent}</ul>
    <h3>Efecto de cada interrupción frente a continuar la tendencia previa</h3>${table(['Interrupción', 'Momento', 'Ajustado', 'Contrafactual', 'Diferencia (IC)', 'Cambio relativo (IC)', 'P'], effRows)}
    <h3>Coeficientes</h3>${table(['Parámetro', 'Estimación', 'IC', 'P', ...(exp ? [`exp(β) · ${nm}`] : [])], coefRows)}
    ${fit.diag ? `<h3>Diagnóstico de los residuos</h3>${howto('Autocorrelación (ACF) y autocorrelación parcial (PACF) de los residuos con la banda ±1,96/√n.', 'Comprobar si los errores están correlacionados en el tiempo: si lo están, los errores estándar clásicos son demasiado pequeños.', fit.diag.lb.p < 0.05 ? 'Hay autocorrelación (Ljung–Box P &lt; 0,05): use Newey–West o AR(1), no «ninguno». Si persisten picos estacionales, añada armónicos de Fourier.' : 'Sin evidencia de autocorrelación residual (Ljung–Box P ≥ 0,05).')}<div class="plots single">${figWrap(acfPlot({ acf: fit.diag.acf, pacf: fit.diag.pacf, band: fit.diag.band, lb: fit.diag.lb, dw: fit.diag.dw }), 'acf_pacf_residuos')}</div>` : ''}`);
  bindFigs($('#r-out'));
}
function sentence(e, fam, fit) {
  const a = e.mean, f = fit.family, rel = pc(a.rel.est), ab = `${dash(a.abs.est)} (IC ${dash(a.abs.lo)} a ${dash(a.abs.hi)})`, dirw = a.abs.est < 0 ? 'menor' : 'mayor', sig = a.p < 0.05 ? 'estadísticamente distinto' : 'sin evidencia clara de diferencia';
  const lv = fit.coef[e.level], sl = fit.coef[e.slope], lt = lv ? `Cambio inmediato de nivel: ${fam === 'gaussian' ? dash(lv.est) : `${f2(lv.exp.est, 2)} (${f2(lv.exp.lo, 2)} a ${f2(lv.exp.hi, 2)})`}.` : '', st = sl ? ` Cambio de pendiente por periodo: ${fam === 'gaussian' ? dash(sl.est) : f2(sl.exp.est, 3)} (P ${fp(sl.p)}).` : '';
  return `${lt}${st} En promedio, tras el periodo de transición, el valor fue ${pc(Math.abs(a.rel.est))} ${dirw} que el contrafactual (${ab}; ${sig}).`;
}

// ---------------------------------------------------------------- pronóstico
const FC = { sel: ['linear', 'ets', 'arima'], train: 'previous', K: 2, show: null };
function renderFcForm() {
  if (!S.raw.length) return; gate(); const P = S.opt.P || PERIOD[S.ax.unit] || 0;
  $('#fc-form').innerHTML = `<div class="wide"><b>Métodos a comparar</b><div class="row gap wrap">${METHODS.map(([id, t]) => `<label class="inline"><input type="checkbox" data-m="${id}"${FC.sel.includes(id) ? ' checked' : ''}> ${esc(t)}</label>`).join('')}</div></div>
    <label>Datos de entrenamiento<select id="fc-train">${[['previous', 'Solo el tramo inmediatamente anterior'], ['all', 'Todo lo anterior a la interrupción']].map(([v, t]) => opt(v, t, v === FC.train)).join('')}</select></label><label>Armónicos de Fourier (estacionalidad)<select id="fc-K">${[0, 1, 2, 3].map(v => opt(v, v ? `${v}` : 'Ninguno', v === FC.K)).join('')}</select></label><label>Periodos por ciclo<input id="fc-P" type="number" min="0" value="${P}"></label>`;
  $$('#fc-form [data-m]').forEach(c => (c.onchange = () => (FC.sel = $$('#fc-form [data-m]').filter(x => x.checked).map(x => x.dataset.m)))); $('#fc-train').onchange = e => (FC.train = e.target.value); $('#fc-K').onchange = e => (FC.K = +e.target.value);
}
$('#fc-run').addEventListener('click', async () => { $('#fc-msg').textContent = 'Calculando…'; await new Promise(r => setTimeout(r, 30)); try { runForecast(); $('#fc-msg').textContent = ''; } catch (e) { $('#fc-msg').textContent = ''; $('#fc-out').innerHTML = card('No se pudo pronosticar', errBox(e)); console.error(e); } });
function runForecast() {
  if (!FC.sel.length) throw new Error('Elija al menos un método.'); const g = S.cur ?? (groupNames()[0] ?? null), sr = series(g), obs = sr.obs, P = +($('#fc-P')?.value || 0), K = FC.K, level = S.opt.level, minTrain = 8;
  if (!S.brk.length) throw new Error('Defina al menos una interrupción.'); const T0 = S.brk[0].T; if (T0 < 10) throw new Error('La primera interrupción necesita al menos 10 periodos previos para pronosticar.');
  const base = { P: P > 1 ? P : 0, K: P > 1 ? K : 0, level }, pre = obs.slice(0, T0), h = Math.max(3, Math.min(6, Math.floor(pre.length / 4)));
  const res = FC.sel.map(m => { let acc, ce = null, err = null; try { acc = rollingOrigin(pre, m, { h, ...base }); ce = counterfactualEffects(obs, { breaks: brkIdx(), lag: S.lag, spec: { method: m, ...base }, train: FC.train, minTrain }); } catch (e) { err = e.message; acc = { rmse: NaN, mae: NaN, cover: NaN, n: 0 }; } return { m, name: METHODS.find(x => x[0] === m)[1], acc, ce, err }; });
  const ok = res.filter(r => r.ce && Number.isFinite(r.acc.rmse)), best = (ok.length ? ok : res.filter(r => r.ce)).sort((a, b) => (a.acc.rmse || 1e9) - (b.acc.rmse || 1e9))[0]; if (!best) throw new Error(res.map(r => `${r.name}: ${r.err}`).join(' · '));
  FC.show = FC.show && res.find(r => r.m === FC.show && r.ce) ? FC.show : best.m; S.fc = { g, res, best: best.m, obs, h, sr }; renderForecast();
}
function fcPlots(r, obs, sr) {
  const cfs = r.ce.map(e => ({ j: e.j, pts: e.t.map((t, i) => ({ t, est: e.fc.fc[i], lo: e.fc.lo[i], hi: e.fc.hi[i] })) })), breaks = S.brk.map(b => ({ T: b.T, label: b.name }));
  const notes = r.ce.map(e => { const t = e.next - 1; return e.last ? { j: e.j, t: e.last ? t : t, yo: obs[t] ?? e.last.obs, yc: e.last.cf, text: `${dash(e.last.abs.est)} (${dash(e.last.abs.lo)}, ${dash(e.last.abs.hi)})` } : null; }).filter(Boolean);
  const main = itsPlot({ obs, se: null, ticks: S.tk.ticks, minor: S.tk.minor, ticksFull: Array.from({ length: S.ax.N }, (_, i) => S.ax.label(i)), breaks, lag: S.lag, segs: [], cfs, notes, forecast: true, segNames: ['Before ' + breaks[0].label, ...breaks.map(b => `After ${b.label}`)], cfLabel: `Forecast (${METHOD_EN[r.m]})`, bandLabel: '95% prediction interval', ylab: yl(), zero: S.opt.zero, xlab: '' });
  const diff = diffPlot({ title: 'Observed − forecast (counterfactual)', N: S.ax.N, ticks: S.tk.ticks, ylab: `Observed − forecast (${yl()})`, series: r.ce.map(e => ({ j: e.j, label: breaks[e.j].label, pts: e.t.map((t, i) => ({ t, d: obs[t] - e.fc.fc[i], lo: obs[t] - e.fc.hi[i], hi: obs[t] - e.fc.lo[i] })).filter(p => Number.isFinite(p.d) && p.t >= e.t0) })) });
  return { main, diff };
}
function renderForecast() {
  const F = S.fc, cur = F.res.find(r => r.m === FC.show), acc = F.res.map(r => ({ label: METHOD_EN[r.m], rmse: r.acc.rmse, cover: r.acc.cover })), p = fcPlots(cur, F.obs, F.sr);
  const rows = cur.ce.flatMap(e => [['Promedio del periodo', e.mean], ['Último periodo', e.last]].filter(([, c]) => c).map(([w, c]) => [esc(S.brk[e.j].name), w, esc(e.name), dash(c.obs), dash(c.cf ?? c.cf), ci(c.abs.est, c.abs.lo, c.abs.hi, dash), `${pc(c.rel.est)} (${pc(c.rel.lo)}, ${pc(c.rel.hi)})`]));
  const cmp = F.res.filter(r => r.ce).map(r => [esc(r.name), ...r.ce.map(e => `${dash(e.mean.abs.est)} (${dash(e.mean.abs.lo)}, ${dash(e.mean.abs.hi)}) · ${pc(e.mean.rel.est)}`)]);
  $('#fc-out').innerHTML = card('Resultado del pronóstico', `${howto('Arriba, la precisión de cada método pronosticando periodos que ya se conocen del tramo previo. Abajo, el contrafactual del método elegido y la diferencia con lo observado.', 'Contrastar si el efecto estimado depende del método usado para construir el contrafactual.', 'El IC de la diferencia es conservador (supone errores de pronóstico perfectamente correlacionados entre periodos). Si el método más preciso y la regresión segmentada coinciden, la conclusión es robusta.')}
    <div class="plots single">${figWrap(accuracyPlot({ rows: acc, best: METHOD_EN[F.best] }), 'precision_pronostico')}</div>
    <div class="row gap wrap mt"><label class="inline">Método mostrado <select id="fc-show">${F.res.filter(r => r.ce).map(r => opt(r.m, `${r.name}${r.m === F.best ? ' (más preciso)' : ''}`, r.m === FC.show)).join('')}</select></label></div>
    <div class="plots single">${figWrap(p.main, 'pronostico_contrafactual')}${figWrap(p.diff, 'diferencia_observado_pronostico')}</div>
    <h3>Efecto frente al pronóstico (${esc(cur.name)})</h3>${table(['Interrupción', 'Ventana', 'Modelo', 'Observado', 'Pronóstico', 'Diferencia (IC)', 'Cambio relativo (IC)'], rows)}
    <h3>Comparación de métodos (diferencia media tras cada interrupción)</h3>${table(['Método', ...S.brk.map(b => esc(b.name))], cmp)}${F.res.filter(r => r.err).map(r => `<p class="warn">${esc(r.name)}: ${esc(r.err)}</p>`).join('')}`);
  bindFigs($('#fc-out')); $('#fc-show').addEventListener('change', e => { FC.show = e.target.value; renderForecast(); });
}

// ---------------------------------------------------------------- grupos y paneles
const GP = { j: 0, par: 'level', labels: 'ends', cols: 4, shareY: true };
export function allGroups() { const gs = groupNames(); return gs.map(g => ({ g, ...safeRun(g) })); }
function safeRun(g) { try { return runOne(g); } catch (e) { return { err: e.message, sr: series(g) }; } }
function renderGroups() {
  const gs = groupNames(); gate(); if (!S.raw.length) return;
  if (!gs.length) { $('#g-ctl').innerHTML = '<h2>Grupos y paneles</h2><p class="muted">Los datos no tienen una columna de grupo. Elíjala en el paso 2 (p. ej., región o centro) para ver un panel por grupo, un forest de efectos y el efecto agrupado.</p>'; $('#g-out').innerHTML = ''; return; }
  const nb = S.brk.length;
  $('#g-ctl').innerHTML = `<h2>Grupos y paneles</h2><div class="howto"><div><b>Qué se ve</b>La serie de cada grupo en paneles pequeños (antes en negro, después en color) y el efecto de una interrupción en cada grupo con el efecto agrupado.</div><div><b>Para qué</b>Ver si el efecto es homogéneo entre regiones o centros y obtener un efecto resumen con un modelo de efectos aleatorios.</div><div><b>Cómo interpretar</b>Un I² alto indica que los efectos difieren entre grupos más de lo esperable por azar: el efecto agrupado debe leerse con cautela.</div></div>
    <div class="form"><label>Interrupción<select id="g-j">${S.brk.map((b, j) => opt(j, b.name, j === GP.j)).join('')}</select></label><label>Efecto<select id="g-par">${[['level', 'Cambio de nivel'], ['slope', 'Cambio de pendiente'], ['rel', 'Cambio relativo medio vs contrafactual']].map(([v, t]) => opt(v, t, v === GP.par)).join('')}</select></label><label>Valores en los paneles<select id="g-lab">${[['ends', 'Primero, último y en la interrupción'], ['all', 'Todos'], ['none', 'Ninguno']].map(([v, t]) => opt(v, t, v === GP.labels)).join('')}</select></label><label>Columnas<select id="g-cols">${[2, 3, 4, 5, 6].map(v => opt(v, v, v === GP.cols)).join('')}</select></label><label class="inline"><input type="checkbox" id="g-share"${GP.shareY ? ' checked' : ''}> Misma escala Y</label></div><div class="row gap mt wrap"><button class="primary" id="g-run">Ajustar todos los grupos</button><span class="muted" id="g-msg"></span></div>`;
  ['g-j', 'g-par', 'g-lab', 'g-cols', 'g-share'].forEach(id => ($('#' + id).onchange = () => { GP.j = +$('#g-j').value; GP.par = $('#g-par').value; GP.labels = $('#g-lab').value; GP.cols = +$('#g-cols').value; GP.shareY = $('#g-share').checked; if (S.gr) drawGroups(); }));
  $('#g-run').onclick = async () => { $('#g-msg').textContent = 'Ajustando…'; await new Promise(r => setTimeout(r, 30)); try { S.gr = allGroups(); drawGroups(); $('#g-msg').textContent = ''; } catch (e) { $('#g-msg').textContent = ''; $('#g-out').innerHTML = card('Error', errBox(e)); } };
  if (S.gr) drawGroups(); else $('#g-out').innerHTML = '';
}
export function groupRows() {
  const R = S.gr.filter(r => r.fit), j = Math.min(GP.j, S.brk.length - 1), fam = S.opt.family, exp = fam !== 'gaussian';
  return R.map(r => { const e = r.fit.effects[j], c = GP.par === 'level' ? r.fit.coef[e.level] : GP.par === 'slope' ? r.fit.coef[e.slope] : null; if (GP.par === 'rel') return { g: r.g, est: e.mean.rel.est, lo: e.mean.rel.lo, hi: e.mean.rel.hi, se: (e.mean.rel.hi - e.mean.rel.lo) / (2 * r.fit.q) }; if (!c) return { g: r.g, est: NaN, lo: NaN, hi: NaN, se: NaN }; return exp ? { g: r.g, est: c.exp.est, lo: c.exp.lo, hi: c.exp.hi, se: c.se, raw: c.est } : { g: r.g, est: c.est, lo: c.lo, hi: c.hi, se: c.se }; });
}
function drawGroups() {
  const fam = S.opt.family, exp = fam !== 'gaussian', rows = groupRows().filter(r => Number.isFinite(r.est)), j = Math.min(GP.j, S.brk.length - 1), bad = S.gr.filter(r => r.err);
  const panels = S.gr.filter(r => r.fit).map(r => ({ name: r.g, y: r.sr.obs, breaks: brkIdx(), fit: { segs: r.fit.curves.seg }, cf: r.fit.curves.cf[0]?.pts })), pooled = (() => { try { const o = runOne(null); return o.sr.obs.some(Number.isFinite) ? { name: 'Overall', y: o.sr.obs, breaks: brkIdx(), fit: { segs: o.fit.curves.seg }, cf: o.fit.curves.cf[0]?.pts } : null; } catch { return null; } })();
  const pl = panelPlot({ panels: [...panels, ...(pooled ? [pooled] : [])], cols: GP.cols, ticks: S.tk.ticks, ylab: yl(), labels: GP.labels, shareY: GP.shareY, dec: S.opt.family === 'gaussian' ? 1 : 1 });
  let fr = rows.map(r => ({ label: r.g, est: r.est, lo: r.lo, hi: r.hi })); const tr = exp && GP.par !== 'rel' ? Math.log : v => v;
  let poolTxt = '';
  if (rows.length >= 2) { const p = poolEffects(rows.map(r => (exp && GP.par !== 'rel' ? Math.log(r.est) : r.est)), rows.map(r => r.se)), back = exp && GP.par !== 'rel' ? Math.exp : v => v; fr.push({ label: 'Pooled (random effects)', est: back(p.est), lo: back(p.lo), hi: back(p.hi), pooled: true }); poolTxt = `<p class="sub">Efecto agrupado (efectos aleatorios, DerSimonian–Laird): <b>${dash(back(p.est))}</b> (IC ${dash(back(p.lo))} a ${dash(back(p.hi))}); I² = ${pc(p.I2)}; τ² = ${f2(p.tau2, 3)}; ${p.k} grupos.</p>`; }
  const ttl = `${S.brk[j].name}: ${{ level: 'level change', slope: 'slope change', rel: 'relative change vs counterfactual' }[GP.par]}`;
  const forest = forestPlot({ title: ttl, xlab: exp && GP.par !== 'rel' ? `${fam === 'poisson' ? 'Rate ratio' : 'Odds ratio'} (95% CI)` : GP.par === 'rel' ? 'Relative change (95% CI)' : `Difference in ${yl()} (95% CI)`, rows: fr, log: exp && GP.par !== 'rel', ref: exp && GP.par !== 'rel' ? 1 : 0, dec: GP.par === 'rel' ? 2 : 2 });
  $('#g-out').innerHTML = card('Paneles por grupo', `<div class="plots single">${figWrap(pl, 'paneles_por_grupo')}</div>`) + card('Efecto por grupo y efecto agrupado', `<div class="plots single">${figWrap(forest, 'forest_efectos_por_grupo')}</div>${poolTxt}${bad.length ? `<p class="warn">No se pudo ajustar: ${bad.map(b => `${esc(b.g)} (${esc(b.err)})`).join('; ')}</p>` : ''}${table(['Grupo', 'Estimación', 'IC'], rows.map(r => [esc(r.g), dash(r.est), `${dash(r.lo)} a ${dash(r.hi)}`]))}<div class="row gap mt"><button class="ghost" id="g-csv">Descargar tabla (CSV)</button></div>`);
  bindFigs($('#g-out')); $('#g-csv').onclick = () => download('efectos_por_grupo.csv', toCSV(rows.map(r => ({ group: r.g, estimate: r.est, ci_low: r.lo, ci_high: r.hi })), ['group', 'estimate', 'ci_low', 'ci_high']), 'text/csv');
}

gate();
window.__its = { S, series, runOne, FC, GP, fitOpts, showTab };
