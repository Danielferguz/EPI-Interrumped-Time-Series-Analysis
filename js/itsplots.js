// Figuras de series temporales interrumpidas (SVG, texto en inglés para publicación):
// serie segmentada con IC y contrafactual, pronóstico con intervalo y área de diferencia, paneles por grupo, diferencia observada − contrafactual,
// ACF/PACF de residuos, forest de efectos por grupo y precisión de los métodos de pronóstico.
export const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
export const SEG_COLORS = ['#2F6DB3', '#D1495B', '#2A9D8F', '#8E5BB5', '#E08E2B'], GREY = '#7A8286', INKC = 'currentColor';
const STYLE = '<style>text{font:12px system-ui,sans-serif;fill:currentColor}.h{font-weight:600}.m{fill:var(--muted,#666)}.s{font-size:11px}.t{font-size:13px;font-weight:600}.g{stroke:#8F9699;stroke-opacity:.22;stroke-width:1}.g2{stroke:#8F9699;stroke-opacity:.1;stroke-width:.8}</style>';
const f1 = (v, d = 1) => (Number.isFinite(v) ? (Math.abs(v) >= 1000 ? v.toFixed(0) : v.toFixed(d)).replace('-', '−') : '–');
export function niceTicks(lo, hi, n = 5) {
  const span = hi - lo || 1, step0 = Math.pow(10, Math.floor(Math.log10(span / n))), m = span / step0 > 5 * n / 1.2 ? 5 : span / step0 > 2 * n / 1.2 ? 2 : 1, st = step0 * m, out = [];
  for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + 1e-9 * span; v += st) out.push(+v.toFixed(10)); return out;
}
const lin = (lo, hi, a, b) => v => a + (v - lo) / ((hi - lo) || 1) * (b - a);
const path = (pts, sx, sy) => pts.map((p, i) => `${i ? 'L' : 'M'}${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`).join('');
const band = (pts, sx, sy) => (pts.length ? `M${pts.map(p => `${sx(p.t).toFixed(1)},${sy(p.hi).toFixed(1)}`).join('L')}L${[...pts].reverse().map(p => `${sx(p.t).toFixed(1)},${sy(p.lo).toFixed(1)}`).join('L')}Z` : '');
const textW = (t, fs = 12) => String(t).length * fs * 0.56;
// rotula las marcas del eje temporal: ticks = [{i, text}]
function xAxis(ticks, sx, y0, y1, minor, H0) {
  let s = ''; (minor || []).forEach(i => (s += `<line class="g2" x1="${sx(i)}" x2="${sx(i)}" y1="${y1}" y2="${y0}"/>`));
  ticks.forEach(k => (s += `<line class="g" x1="${sx(k.i)}" x2="${sx(k.i)}" y1="${y1}" y2="${y0}"/><line x1="${sx(k.i)}" x2="${sx(k.i)}" y1="${y0}" y2="${y0 + 4}" stroke="currentColor"/><text x="${sx(k.i) + 4}" y="${y0 + 15}" text-anchor="end" transform="rotate(-30 ${sx(k.i) + 4} ${y0 + 15})" class="s">${esc(k.text)}</text>`));
  return s;
}
const legendRow = (items, W, y) => { // items: [{label, color, kind: 'line'|'band'|'dash'}]; se reparte en varias filas centradas si no cabe
  const w = it => 26 + textW(it.label, 12) + 22, rows = [[]]; let cur = 0;
  for (const it of items) { if (cur + w(it) > W - 16 && rows[rows.length - 1].length) { rows.push([]); cur = 0; } rows[rows.length - 1].push(it); cur += w(it); }
  let s = ''; rows.forEach((row, q) => {
    let x = Math.max(8, (W - row.reduce((a, it) => a + w(it), 0)) / 2); const yy = y + q * 20;
    for (const it of row) {
      if (it.kind === 'band') s += `<rect x="${x}" y="${yy - 6}" width="22" height="12" fill="${it.color}" fill-opacity=".22"/>`;
      else s += `<line x1="${x}" x2="${x + 22}" y1="${yy}" y2="${yy}" stroke="${it.color}" stroke-width="2.8" ${it.kind === 'dash' ? 'stroke-dasharray="6 4"' : ''}/>${it.kind === 'line' ? `<circle cx="${x + 11}" cy="${yy}" r="3.6" fill="${it.color}" stroke="#fff" stroke-width="1"/>` : ''}`;
      s += `<text x="${x + 28}" y="${yy + 4}">${esc(it.label)}</text>`; x += w(it);
    }
  });
  return { s, rows: rows.length };
};

// ---------------------------------------------------------------- serie interrumpida: segmentos, contrafactual, efectos
// o: {title, ylab, obs[], se[]?, ticks, minor, breaks:[{T, label}], lag, segs:[{s, pts:[{t,est,lo,hi}]}], cfs:[{j, pts}], notes:[{j, t, yo, yc, text}],
//     segNames[], colors[], zero, W, H, forecast: bool (observado como línea y área de diferencia), N}
export function itsPlot(o) {
  const W = o.W || 780, H = o.H || 440, N = o.obs.length, col = o.colors || SEG_COLORS, nb = o.breaks.length, ml = 62, mr = 22, mb = 74, items = [...(o.segNames || []).slice(0, nb + 1).map((n, i) => ({ label: n, color: col[i % col.length], kind: 'line' })), ...((o.cfs || []).length ? [{ label: o.cfLabel || 'Counterfactual (previous trend)', color: GREY, kind: 'dash' }, { label: o.bandLabel || '95% CI', color: GREY, kind: 'band' }] : [])], leg = legendRow(items, W, o.title ? 38 : 16), mt = (o.title ? 56 : 40) + 6 + (leg.rows - 1) * 20;
  const all = [...o.obs.filter(Number.isFinite), ...(o.segs || []).flatMap(g => g.pts.flatMap(p => [p.lo, p.hi])), ...(o.cfs || []).flatMap(g => g.pts.flatMap(p => [p.lo, p.hi]))].filter(Number.isFinite);
  let lo = Math.min(...all), hi = Math.max(...all); const pad = (hi - lo) * 0.06 || 1; lo = o.zero ? Math.min(0, lo) : lo - pad; hi += pad; const yt = niceTicks(lo, hi, 6); if (!o.zero) lo = Math.min(lo, yt[0]); hi = Math.max(hi, yt[yt.length - 1]);
  const sx = lin(-0.5, N - 0.5, ml, W - mr), sy = lin(lo, hi, H - mb, mt);
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="plot" role="img">${STYLE}`;
  if (o.title) s += `<text x="8" y="16" class="t">${esc(o.title)}</text>`;
  s += leg.s;
  yt.filter(v => v >= lo - 1e-9 && v <= hi + 1e-9).forEach(v => (s += `<line class="g" x1="${ml}" x2="${W - mr}" y1="${sy(v)}" x2="${W - mr}" y2="${sy(v)}"/><text x="${ml - 7}" y="${sy(v) + 4}" text-anchor="end" class="s">${f1(v, Math.abs(yt[1] - yt[0]) < 1 ? 1 : 0)}</text><line x1="${ml - 4}" x2="${ml}" y1="${sy(v)}" y2="${sy(v)}" stroke="currentColor"/>`));
  s += xAxis(o.ticks, sx, H - mb, mt, o.minor, H);
  s += `<rect x="${ml}" y="${mt}" width="${W - ml - mr}" height="${H - mb - mt}" fill="none" stroke="currentColor" stroke-width=".8"/>`;
  // periodo de transición
  o.breaks.forEach(b => { if (o.lag > 0) s += `<rect x="${sx(b.T - 0.5)}" y="${mt}" width="${sx(b.T - 0.5 + o.lag) - sx(b.T - 0.5)}" height="${H - mb - mt}" fill="#8F9699" fill-opacity=".16"/>`; });
  // contrafactual: banda y línea punteada
  (o.cfs || []).forEach(g => { const c = GREY; s += `<path d="${band(g.pts, sx, sy)}" fill="${c}" fill-opacity=".2"/><path d="${path(g.pts.map(p => [p.t, p.est]), sx, sy)}" fill="none" stroke="${c}" stroke-width="2.4" stroke-dasharray="7 5"/>`; });
  // área entre lo observado y el contrafactual (modo pronóstico)
  if (o.forecast) (o.cfs || []).forEach(g => { const pts = g.pts.filter(p => Number.isFinite(o.obs[p.t])); if (pts.length > 1) s += `<path d="M${pts.map(p => `${sx(p.t).toFixed(1)},${sy(o.obs[p.t]).toFixed(1)}`).join('L')}L${[...pts].reverse().map(p => `${sx(p.t).toFixed(1)},${sy(p.est).toFixed(1)}`).join('L')}Z" fill="${col[(g.j + 1) % col.length]}" fill-opacity=".14"/>`; });
  // ajuste por segmento
  (o.segs || []).forEach(g => { const c = col[g.s % col.length]; s += `<path d="${band(g.pts, sx, sy)}" fill="${c}" fill-opacity=".17"/><path d="${path(g.pts.map(p => [p.t, p.est]), sx, sy)}" fill="none" stroke="${c}" stroke-width="3"/>`; });
  // interrupciones
  o.breaks.forEach((b, j) => {
    const x = sx(b.T - 0.5); s += `<line x1="${x}" x2="${x}" y1="${mt}" y2="${H - mb}" stroke="#222" stroke-width="1.6" stroke-dasharray="2 5" stroke-linecap="round"/>`;
    const tw = textW(b.label, 11) + 10, tx = Math.min(Math.max(x, ml + tw / 2), W - mr - tw / 2); s += `<rect x="${tx - tw / 2}" y="${mt + 4 + (j % 2) * 18}" width="${tw}" height="16" rx="3" fill="var(--card,#fff)" fill-opacity=".92" stroke="#222" stroke-width=".6"/><text x="${tx}" y="${mt + 16 + (j % 2) * 18}" text-anchor="middle" class="s h">${esc(b.label)}</text>`;
  });
  // observaciones
  const segOf = t => o.breaks.filter(b => t >= b.T).length;
  if (o.forecast) { const lines = {}; o.obs.forEach((v, t) => { if (Number.isFinite(v)) (lines[segOf(t)] ||= []).push([t, v]); }); Object.entries(lines).forEach(([k, pts]) => (s += `<path d="${path(pts, sx, sy)}" fill="none" stroke="${col[k % col.length]}" stroke-width="2.6"/>`)); }
  o.obs.forEach((v, t) => { if (!Number.isFinite(v)) return; const c = col[segOf(t) % col.length], e = o.se?.[t]; if (Number.isFinite(e)) s += `<line x1="${sx(t)}" x2="${sx(t)}" y1="${sy(v - 1.96 * e)}" y2="${sy(v + 1.96 * e)}" stroke="${c}" stroke-opacity=".32" stroke-width="1.2"/>`; s += `<circle cx="${sx(t)}" cy="${sy(v)}" r="${o.forecast ? 3 : 3.8}" fill="${c}" stroke="#fff" stroke-width="1"><title>${esc(o.ticksFull?.[t] ?? t)}: ${f1(v, 2)}</title></circle>`; });
  // efecto estimado: puente entre ajustado y contrafactual al final de cada segmento
  (o.notes || []).forEach(n => {
    const x = sx(n.t), y1 = sy(n.yo), y2 = sy(n.yc), c = col[(n.j + 1) % col.length]; s += `<path d="M${x},${y1}L${x},${y2}M${x - 5},${y1}h10M${x - 5},${y2}h10" stroke="${c}" stroke-width="2" fill="none"/>`;
    const tw = textW(n.text, 11.5) + 12, bx = x + 10 + tw > W - mr - 2 ? x - 10 - tw : x + 10, by = Math.abs(y1 - y2) < 34 ? Math.max(Math.min(y1, y2) - 34, mt + 4) : Math.min(Math.max((y1 + y2) / 2 - 11, mt + 40), H - mb - 26); s += `<rect x="${bx}" y="${by}" width="${tw}" height="22" rx="4" fill="var(--card,#fff)" fill-opacity=".94" stroke="${c}" stroke-width="1.2"/><text x="${bx + 6}" y="${by + 15}" class="s h" style="font-size:11.5px">${esc(n.text)}</text>`;
  });
  s += `<text x="${(ml + W - mr) / 2}" y="${H - 6}" text-anchor="middle" class="m">${esc(o.xlab || '')}</text><text transform="translate(14 ${(mt + H - mb) / 2}) rotate(-90)" text-anchor="middle">${esc(o.ylab || '')}</text>`;
  return s + '</svg>';
}

const thin = (a, k) => (a.length <= k ? a : a.filter((_, i) => i % Math.ceil(a.length / k) === 0));
// ---------------------------------------------------------------- paneles por grupo (como multiples pequeños)
// o: {panels:[{name, y[], breaks:[T..], fit?:{segs}}], cols, ticks, ylab, labels:'all'|'ends'|'none', shareY, colors, dec, scaleToZero}
export function panelPlot(o) {
  const P = o.panels, nc = Math.min(o.cols || 4, P.length), nr = Math.ceil(P.length / nc), W = o.W || 1000, gx = 8, gy = 10, ml = 54, mr = 8, hdr = 22, cw = (W - ml - mr - gx * (nc - 1)) / nc, ph = o.panelH || 96, xb = 40, H = 8 + nr * (hdr + ph + gy) + xb, N = P[0].y.length, pre = (o.colors || ['#222', '#2B3F8F'])[0], post = (o.colors || ['#222', '#2B3F8F'])[1], dec = o.dec ?? 1;
  const gl = (() => { const all = P.flatMap(p => p.y.filter(Number.isFinite)); return [Math.min(0, ...all), Math.max(...all)]; })();
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="plot" role="img">${STYLE}<text transform="translate(13 ${H / 2 - 10}) rotate(-90)" text-anchor="middle" class="t">${esc(o.ylab || '')}</text>`;
  P.forEach((p, k) => {
    const c = k % nc, r = Math.floor(k / nc), x0 = ml + c * (cw + gx), y0 = 8 + r * (hdr + ph + gy), fin = p.y.filter(Number.isFinite), lo = o.shareY === false ? Math.min(0, ...fin) : gl[0], hi = o.shareY === false ? Math.max(...fin) : gl[1], rng = hi - lo || 1, top = y0 + hdr, sx = lin(-0.5, N - 0.5, x0 + 6, x0 + cw - 6), sy = lin(lo - rng * 0.08, hi + rng * 0.22, top + ph, top);
    s += `<rect x="${x0}" y="${y0}" width="${cw}" height="${hdr}" fill="#D9D9D9" fill-opacity=".85" stroke="#8F9699" stroke-width=".7"/><text x="${x0 + cw / 2}" y="${y0 + 15}" text-anchor="middle" class="h">${esc(p.name)}</text><rect x="${x0}" y="${top}" width="${cw}" height="${ph}" fill="none" stroke="#8F9699" stroke-width=".7"/>`;
    niceTicks(lo, hi, 3).filter(v => v >= lo && v <= hi).forEach(v => { s += `<line class="g" x1="${x0}" x2="${x0 + cw}" y1="${sy(v)}" y2="${sy(v)}"/>`; if (c === 0) s += `<text x="${x0 - 4}" y="${sy(v) + 4}" text-anchor="end" class="s m">${f1(v, 0)}</text>`; });
    (thin(o.ticks || [], 5)).forEach(t => { s += `<line class="g2" x1="${sx(t.i)}" x2="${sx(t.i)}" y1="${top}" y2="${top + ph}"/>`; if (r === nr - 1 || k + nc >= P.length) s += `<text x="${sx(t.i)}" y="${top + ph + 13}" text-anchor="middle" class="s m">${esc(/\d{4}$/.test(t.text) ? t.text.slice(-4) : t.text)}</text>`; });
    const T = p.breaks?.[0] ?? N, pts = p.y.map((v, t) => [t, v]).filter(q => Number.isFinite(q[1])), preP = pts.filter(q => q[0] < T), postP = pts.filter(q => q[0] >= T - 1);
    if (p.fit) p.fit.segs.forEach(g => { s += `<path d="${path(g.pts.map(q => [q.t, q.est]), sx, sy)}" fill="none" stroke="${g.s === 0 ? pre : post}" stroke-opacity=".28" stroke-width="5" stroke-linecap="round"/>`; });
    s += `<path d="${path(preP, sx, sy)}" fill="none" stroke="${pre}" stroke-width="2"/><path d="${path(postP, sx, sy)}" fill="none" stroke="${post}" stroke-width="2"/>`;
    if (p.cf) s += `<path d="${path(p.cf.map(q => [q.t, q.est]), sx, sy)}" fill="none" stroke="${GREY}" stroke-width="1.6" stroke-dasharray="5 4"/>`;
    pts.forEach(([t, v]) => (s += `<circle cx="${sx(t)}" cy="${sy(v)}" r="2.4" fill="${t >= T ? post : pre}"/>`));
    const lab = o.labels || 'ends'; if (lab !== 'none') pts.forEach(([t, v], i) => { if (lab === 'ends' && !(i === 0 || i === pts.length - 1 || t === T)) return; const up = lab === 'all' ? i % 2 === 0 : v < (lo + hi) / 2 ? true : false; s += `<text x="${Math.min(Math.max(sx(t), x0 + 14), x0 + cw - 14)}" y="${sy(v) + (up ? -6 : 14)}" text-anchor="middle" class="s" style="font-size:10.5px;fill:${t >= T ? post : 'currentColor'}">${f1(v, dec)}</text>`; });
    (p.breaks || []).forEach(b => (s += `<line x1="${sx(b - 0.5)}" x2="${sx(b - 0.5)}" y1="${top}" y2="${top + ph}" stroke="#222" stroke-width="1" stroke-dasharray="2 4"/>`));
  });
  return s + '</svg>';
}

// ---------------------------------------------------------------- diferencia observada − contrafactual
// o: {title, ylab, series:[{j, label, pts:[{t, d, lo, hi}]}], ticks, N, colors, rel?}
export function diffPlot(o) {
  const W = o.W || 780, H = o.H || 330, ml = 62, mr = 22, mt = (o.title ? 34 : 18), mb = 62, N = o.N, col = o.colors || SEG_COLORS, all = o.series.flatMap(g => g.pts.flatMap(p => [p.lo, p.hi, p.d])).filter(Number.isFinite);
  let lo = Math.min(0, ...all), hi = Math.max(0, ...all); const pad = (hi - lo) * 0.08 || 1; lo -= pad; hi += pad; const sx = lin(-0.5, N - 0.5, ml, W - mr), sy = lin(lo, hi, H - mb, mt), yt = niceTicks(lo, hi, 5);
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="plot" role="img">${STYLE}${o.title ? `<text x="8" y="16" class="t">${esc(o.title)}</text>` : ''}`;
  yt.forEach(v => (s += `<line class="g" x1="${ml}" x2="${W - mr}" y1="${sy(v)}" y2="${sy(v)}"/><text x="${ml - 7}" y="${sy(v) + 4}" text-anchor="end" class="s">${f1(v, Math.abs(yt[1] - yt[0]) < 1 ? 1 : 0)}</text>`));
  s += xAxis(o.ticks, sx, H - mb, mt) + `<rect x="${ml}" y="${mt}" width="${W - ml - mr}" height="${H - mb - mt}" fill="none" stroke="currentColor" stroke-width=".8"/><line x1="${ml}" x2="${W - mr}" y1="${sy(0)}" y2="${sy(0)}" stroke="#222" stroke-width="1.4"/>`;
  o.series.forEach(g => { const c = col[(g.j + 1) % col.length]; s += `<path d="${band(g.pts.map(p => ({ t: p.t, lo: p.lo, hi: p.hi })), sx, sy)}" fill="${c}" fill-opacity=".2"/><path d="${path(g.pts.map(p => [p.t, p.d]), sx, sy)}" fill="none" stroke="${c}" stroke-width="2.8"/>`; const m = g.pts[Math.floor(g.pts.length / 2)]; if (m) s += `<text x="${sx(m.t)}" y="${sy(m.hi) - 6}" text-anchor="middle" class="s h" style="fill:${c}">${esc(g.label)}</text>`; });
  s += `<text x="${(ml + W - mr) / 2}" y="${H - 6}" text-anchor="middle" class="m">${esc(o.xlab || '')}</text><text transform="translate(14 ${(mt + H - mb) / 2}) rotate(-90)" text-anchor="middle">${esc(o.ylab || 'Observed − counterfactual')}</text>`;
  return s + '</svg>';
}

// ---------------------------------------------------------------- ACF y PACF de los residuos
export function acfPlot({ acf, pacf, band, title = 'Residual autocorrelation', lb, dw }) {
  const W = 760, H = 270, K = acf.length - 1, ml = 50, mr = 18, mt = 46, mb = 40, pw = (W - ml - mr - 24) / 2;
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="plot" role="img">${STYLE}<text x="8" y="16" class="t">${esc(title)}</text><text x="8" y="33" class="s m">${lb ? `Ljung–Box Q(${lb.lag}) = ${f1(lb.q, 2)}, P ${lb.p < 0.001 ? '< 0.001' : '= ' + lb.p.toFixed(3)}` : ''}${dw != null ? ` · Durbin–Watson = ${f1(dw, 2)}` : ''}</text>`;
  [['ACF', acf, 0], ['PACF', pacf, 1]].forEach(([name, v, c]) => {
    const x0 = ml + c * (pw + 24), sx = lin(0.5, K + 0.5, x0, x0 + pw), sy = lin(-1, 1, H - mb, mt);
    s += `<text x="${x0 + pw / 2}" y="${mt - 6}" text-anchor="middle" class="h">${name}</text><rect x="${x0}" y="${mt}" width="${pw}" height="${H - mb - mt}" fill="none" stroke="currentColor" stroke-width=".7"/>`;
    [-1, -0.5, 0, 0.5, 1].forEach(t => (s += `<line class="g" x1="${x0}" x2="${x0 + pw}" y1="${sy(t)}" y2="${sy(t)}"/>${c === 0 ? `<text x="${x0 - 6}" y="${sy(t) + 4}" text-anchor="end" class="s">${t.toFixed(1).replace('-', '−')}</text>` : ''}`));
    s += `<rect x="${x0}" y="${sy(band)}" width="${pw}" height="${sy(-band) - sy(band)}" fill="#2F6DB3" fill-opacity=".1"/><line x1="${x0}" x2="${x0 + pw}" y1="${sy(0)}" y2="${sy(0)}" stroke="currentColor"/>`;
    for (let k = 1; k <= K; k++) { const out = Math.abs(v[k]) > band; s += `<line x1="${sx(k)}" x2="${sx(k)}" y1="${sy(0)}" y2="${sy(v[k])}" stroke="${out ? '#D1495B' : '#2F6DB3'}" stroke-width="3.4"/>`; if (k % Math.ceil(K / 8) === 0 || k === 1) s += `<text x="${sx(k)}" y="${H - mb + 14}" text-anchor="middle" class="s">${k}</text>`; }
    s += `<text x="${x0 + pw / 2}" y="${H - 6}" text-anchor="middle" class="m">Lag (periods)</text>`;
  });
  return s + '</svg>';
}

// ---------------------------------------------------------------- forest de efectos por grupo
// o: {title, xlab, ref, rows:[{label, est, lo, hi, pooled?, group?}], log?, dec}
export function forestPlot(o) {
  const rows = o.rows, W = o.W || 700, rh = 24, ml = 190, mr = 150, mt = (o.title ? 40 : 22), mb = 44, H = mt + rows.length * rh + mb, tf = o.log ? Math.log : v => v, all = rows.flatMap(r => [r.lo, r.hi, r.est]).filter(Number.isFinite), ref = o.ref ?? (o.log ? 1 : 0);
  let lo = Math.min(...all, ref), hi = Math.max(...all, ref); const sp = tf(hi) - tf(lo) || 1, a = tf(lo) - sp * 0.06, b = tf(hi) + sp * 0.06, sx = v => ml + (tf(v) - a) / (b - a) * (W - ml - mr), dec = o.dec ?? 2;
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="plot" role="img">${STYLE}${o.title ? `<text x="8" y="16" class="t">${esc(o.title)}</text>` : ''}<text x="${W - mr + 10}" y="${mt - 8}" class="s h">${esc(o.colhead || 'Estimate (95% CI)')}</text>`;
  const ticks = o.log ? [0.25, 0.5, 0.75, 1, 1.5, 2, 4].filter(v => v >= lo * 0.8 && v <= hi * 1.25) : niceTicks(lo, hi, 5);
  ticks.forEach(t => (s += `<line class="g" x1="${sx(t)}" x2="${sx(t)}" y1="${mt}" y2="${H - mb}"/><line x1="${sx(t)}" x2="${sx(t)}" y1="${H - mb}" y2="${H - mb + 4}" stroke="currentColor"/><text x="${sx(t)}" y="${H - mb + 17}" text-anchor="middle" class="s">${+t.toPrecision(3)}</text>`));
  s += `<line x1="${sx(ref)}" x2="${sx(ref)}" y1="${mt}" y2="${H - mb}" stroke="#222" stroke-width="1.2" stroke-dasharray="4 3"/><line x1="${ml}" x2="${W - mr}" y1="${H - mb}" y2="${H - mb}" stroke="currentColor"/>`;
  rows.forEach((r, i) => {
    const y = mt + i * rh + rh / 2, bold = r.pooled; s += `<text x="${ml - 10}" y="${y + 4}" text-anchor="end" ${bold ? 'class="h"' : ''}>${esc(r.label)}</text>`;
    if (Number.isFinite(r.est)) { if (r.pooled) s += `<path d="M${sx(r.lo)},${y}L${sx(r.est)},${y - 8}L${sx(r.hi)},${y}L${sx(r.est)},${y + 8}Z" fill="#8B2E3B"/>`; else s += `<line x1="${sx(r.lo)}" x2="${sx(r.hi)}" y1="${y}" y2="${y}" stroke="#2F6DB3" stroke-width="2"/><rect x="${sx(r.est) - 4.5}" y="${y - 4.5}" width="9" height="9" fill="#2F6DB3"/>`; s += `<text x="${W - mr + 10}" y="${y + 4}" class="s">${f1(r.est, dec)} (${f1(r.lo, dec)}, ${f1(r.hi, dec)})</text>`; }
    if (r.pooled && i > 0) s += `<line x1="8" x2="${W - 8}" y1="${y - rh / 2}" y2="${y - rh / 2}" stroke="#8F9699" stroke-width=".7"/>`;
  });
  return s + `<text x="${(ml + W - mr) / 2}" y="${H - 8}" text-anchor="middle" class="m">${esc(o.xlab || '')}</text></svg>`;
}

// ---------------------------------------------------------------- precisión de los métodos de pronóstico
export function accuracyPlot({ rows, title = 'Forecast accuracy in the pre-interruption period', xlab = 'RMSE (rolling-origin cross-validation)', best = null }) {
  const W = 700, rh = 30, ml = Math.min(330, Math.max(...rows.map(r => textW(r.label, 12))) + 24), mr = 90, mt = 44, mb = 44, H = mt + rows.length * rh + mb, mx = Math.max(...rows.map(r => r.rmse).filter(Number.isFinite)) * 1.12 || 1, sx = v => ml + v / mx * (W - ml - mr);
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" class="plot" role="img">${STYLE}<text x="8" y="16" class="t">${esc(title)}</text><text x="8" y="33" class="s m">Lower is better · coverage = share of observations inside the 95% prediction interval</text>`;
  niceTicks(0, mx, 5).forEach(t => (s += `<line class="g" x1="${sx(t)}" x2="${sx(t)}" y1="${mt}" y2="${H - mb}"/><text x="${sx(t)}" y="${H - mb + 16}" text-anchor="middle" class="s">${+t.toPrecision(3)}</text>`));
  rows.forEach((r, i) => { const y = mt + i * rh + 4, isB = r.label === best; s += `<text x="${ml - 8}" y="${y + 14}" text-anchor="end" ${isB ? 'class="h"' : ''}>${esc(r.label)}</text>`; if (Number.isFinite(r.rmse)) s += `<rect x="${ml}" y="${y}" width="${sx(r.rmse) - ml}" height="${rh - 10}" fill="${isB ? '#2A9D8F' : '#9DB7C9'}"/><text x="${sx(r.rmse) + 6}" y="${y + 14}" class="s">${f1(r.rmse, 2)} · ${Number.isFinite(r.cover) ? Math.round(100 * r.cover) + '%' : ''}</text>`; else s += `<text x="${ml + 6}" y="${y + 14}" class="s m">n/a</text>`; });
  return s + `<line x1="${ml}" x2="${ml}" y1="${mt}" y2="${H - mb}" stroke="currentColor"/><text x="${(ml + W - mr) / 2}" y="${H - 8}" text-anchor="middle" class="m">${esc(xlab)}</text></svg>`;
}
