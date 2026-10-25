// Series temporales interrumpidas (ITS): regresión segmentada con 1 a 3 (o más) interrupciones, periodo de transición (rezago),
// estacionalidad de Fourier, familias gaussiana / Poisson / binomial agrupada, errores autocorrelacionados (Newey-West o AR(1) de Prais-Winsten)
// y datos individuales (modelo mixto o GEE). Sin dependencias del navegador: se usa igual en la app y en las pruebas.
import { symInv, mv, dot, zeros } from './linalg.js?v=2026.10.25-1';
import { qnorm, qt, pchisqU, p2z, p2t } from './dist.js?v=2026.10.25-1';
import { fitLMM, fitGEE } from './models.js?v=2026.10.25-1';

const sum = a => { let s = 0; for (const v of a) s += v; return s; };
const expit = z => 1 / (1 + Math.exp(-z)), logit = p => Math.log(p / (1 - p));
const isNum = v => typeof v === 'number' && isFinite(v);

// ---------------------------------------------------------------- diseño
// breaks: índices (0-based) del primer periodo tras cada interrupción, ascendentes. effects: 'both' | 'level' | 'slope' (o uno por interrupción).
export function buildDesign({ breaks = [], lag = 0, effects = 'both', season = null, labels = [] } = {}) {
  const B = [...breaks].sort((a, b) => a - b), eff = B.map((_, j) => (Array.isArray(effects) ? effects[j] : effects) || 'both'), names = ['Intercepto', 'Tendencia previa (por periodo)'], kind = [{ k: 'int' }, { k: 'trend' }];
  const lab = j => labels[j] || `Interrupción ${j + 1}`;
  B.forEach((T, j) => {
    if (eff[j] !== 'slope') { names.push(`Cambio de nivel · ${lab(j)}`); kind.push({ k: 'step', j }); }
    if (eff[j] !== 'level') { names.push(`Cambio de pendiente · ${lab(j)}`); kind.push({ k: 'slope', j }); }
  });
  const P = season?.P > 1 && season?.K > 0 ? season.P : 0, K = P ? Math.min(season.K, Math.floor(P / 2)) : 0;
  for (let k = 1; k <= K; k++) { if (2 * k === P) { names.push(`Estacional cos(π t)`); kind.push({ k: 'seas', f: k, c: true }); } else { names.push(`Estacional sen ${k}`, `Estacional cos ${k}`); kind.push({ k: 'seas', f: k, c: false }, { k: 'seas', f: k, c: true }); } }
  // fila de diseño en el periodo t; cf = j → anula las interrupciones j, j+1, … (contrafactual: se sigue el régimen anterior); noSeason → sin estacionalidad
  const row = (t, { cf = null, noSeason = false } = {}) => Float64Array.from(kind, d => {
    if (d.k === 'int') return 1; if (d.k === 'trend') return t;
    if (d.k === 'step') return t >= B[d.j] && (cf == null || d.j < cf) ? 1 : 0;
    if (d.k === 'slope') return t >= B[d.j] && (cf == null || d.j < cf) ? t - B[d.j] : 0;
    if (noSeason) return 0; const a = 2 * Math.PI * d.f * t / P; return d.c ? Math.cos(a) : Math.sin(a);
  });
  // segmento del periodo t (0 = antes de la primera interrupción) y si está en el periodo de transición
  const seg = t => B.filter(T => t >= T).length, inLag = t => B.some(T => t >= T && t < T + lag);
  return { breaks: B, effects: eff, names, kind, row, seg, inLag, p: kind.length, lag, P, K, labels: B.map((_, j) => lab(j)) };
}

// ---------------------------------------------------------------- estimación
// MCO con covarianza robusta de Newey-West (HAC) o clásica
function ols(X, y) {
  const n = X.length, p = X[0].length, H = zeros(p, p), b = new Float64Array(p);
  for (let i = 0; i < n; i++) for (let a = 0; a < p; a++) { b[a] += X[i][a] * y[i]; for (let c = 0; c <= a; c++) H[a][c] += X[i][a] * X[i][c]; }
  for (let a = 0; a < p; a++) for (let c = a + 1; c < p; c++) H[a][c] = H[c][a];
  const Hi = symInv(H); if (!Hi) throw new Error('Matriz de diseño singular: pocos periodos en algún segmento o columnas redundantes (revise el rezago, la estacionalidad y las interrupciones).');
  const beta = mv(Hi, b), res = Float64Array.from(y, (v, i) => v - dot(X[i], beta));
  return { beta, Hi, res, n, p };
}
// covarianza HAC: sándwich con núcleo de Bartlett sobre las puntuaciones S (una fila por observación) con los periodos idx
function hac(Hi, S, idx, L) {
  const p = Hi.length, M = zeros(p, p), pos = new Map(idx.map((t, i) => [t, i]));
  for (let i = 0; i < S.length; i++) for (let a = 0; a < p; a++) for (let c = 0; c < p; c++) M[a][c] += S[i][a] * S[i][c];
  for (let l = 1; l <= L; l++) {
    const w = 1 - l / (L + 1);
    for (let i = 0; i < S.length; i++) { const j = pos.get(idx[i] - l); if (j == null) continue; for (let a = 0; a < p; a++) for (let c = 0; c < p; c++) M[a][c] += w * (S[i][a] * S[j][c] + S[j][a] * S[i][c]); }
  }
  const HM = zeros(p, p), V = zeros(p, p);
  for (let a = 0; a < p; a++) for (let c = 0; c < p; c++) { let s = 0; for (let k = 0; k < p; k++) s += Hi[a][k] * M[k][c]; HM[a][c] = s; }
  for (let a = 0; a < p; a++) for (let c = 0; c < p; c++) { let s = 0; for (let k = 0; k < p; k++) s += HM[a][k] * Hi[k][c]; V[a][c] = s; }
  return V;
}
export const nwLag = n => Math.max(1, Math.floor(4 * Math.pow(n / 100, 2 / 9)));
// Prais-Winsten: AR(1) en los errores, iterado; los huecos de la serie reinician la recursión
function praisWinsten(X, y, idx, maxit = 100) {
  const n = X.length, p = X[0].length; let rho = 0, fit = null;
  const consecutive = i => i > 0 && idx[i] - idx[i - 1] === 1;
  for (let it = 0; it < maxit; it++) {
    const Xs = X.map((x, i) => { if (consecutive(i)) return Float64Array.from(x, (v, a) => v - rho * X[i - 1][a]); return Float64Array.from(x, v => v * Math.sqrt(1 - rho * rho)); });
    const ys = Float64Array.from(y, (v, i) => (consecutive(i) ? v - rho * y[i - 1] : v * Math.sqrt(1 - rho * rho)));
    fit = ols(Xs, ys); const e = Float64Array.from(y, (v, i) => v - dot(X[i], fit.beta));
    let num = 0, den = 0; for (let i = 1; i < n; i++) if (consecutive(i)) { num += e[i] * e[i - 1]; den += e[i - 1] * e[i - 1]; }
    const nr = den > 0 ? Math.max(-0.98, Math.min(0.98, num / den)) : 0; if (Math.abs(nr - rho) < 1e-8) { rho = nr; break; } rho = nr;
  }
  const Xs = X.map((x, i) => (consecutive(i) ? Float64Array.from(x, (v, a) => v - rho * X[i - 1][a]) : Float64Array.from(x, v => v * Math.sqrt(1 - rho * rho)))), ys = Float64Array.from(y, (v, i) => (consecutive(i) ? v - rho * y[i - 1] : v * Math.sqrt(1 - rho * rho)));
  fit = ols(Xs, ys); const sig2 = sum(fit.res.map(r => r * r)) / (n - p), V = fit.Hi.map(r => r.map(v => v * sig2));
  return { beta: fit.beta, V, rho, sigma2: sig2, innov: fit.res };
}
// MLG por mínimos cuadrados ponderados iterados: Poisson (con desplazamiento log n) o binomial agrupada (y eventos de m ensayos)
function irls(X, y, fam, m, off) {
  const n = X.length, p = X[0].length; let eta = Float64Array.from(y, (v, i) => (fam === 'poisson' ? Math.log(v + 0.5) : logit((v + 0.5) / (m[i] + 1))) + (fam === 'poisson' ? 0 : 0));
  if (fam === 'poisson' && off) eta = Float64Array.from(eta, (e, i) => e);
  let beta = new Float64Array(p), Hi = null, mu, w, conv = false;
  for (let it = 0; it < 200; it++) {
    mu = Float64Array.from(eta, (e, i) => (fam === 'poisson' ? Math.exp(e) : m[i] * expit(e))); w = Float64Array.from(mu, (u, i) => (fam === 'poisson' ? Math.max(u, 1e-10) : Math.max(m[i] * (u / m[i]) * (1 - u / m[i]), 1e-10)));
    const z = Float64Array.from(eta, (e, i) => e - (off ? off[i] : 0) + (y[i] - mu[i]) / w[i]), H = zeros(p, p), b = new Float64Array(p);
    for (let i = 0; i < n; i++) for (let a = 0; a < p; a++) { b[a] += w[i] * X[i][a] * z[i]; for (let c = 0; c <= a; c++) H[a][c] += w[i] * X[i][a] * X[i][c]; }
    for (let a = 0; a < p; a++) for (let c = a + 1; c < p; c++) H[a][c] = H[c][a];
    Hi = symInv(H); if (!Hi) throw new Error('Matriz de información singular (segmentos con pocos periodos, ceros estructurales o columnas redundantes).');
    const nb = mv(Hi, b), d = Math.max(...nb.map((v, a) => Math.abs(v - beta[a]))); beta = nb; eta = Float64Array.from(X.map(x => dot(x, beta)), (e, i) => e + (off ? off[i] : 0)); if (d < 1e-10) { conv = true; break; }
  }
  mu = Float64Array.from(eta, (e, i) => (fam === 'poisson' ? Math.exp(e) : m[i] * expit(e))); w = Float64Array.from(mu, (u, i) => (fam === 'poisson' ? Math.max(u, 1e-10) : Math.max(u * (1 - u / m[i]), 1e-10)));
  const pearson = Float64Array.from(y, (v, i) => (v - mu[i]) / Math.sqrt(w[i])), phi = sum(pearson.map(r => r * r)) / (n - p);
  return { beta, Hi, mu, w, eta, pearson, phi, conv, resp: Float64Array.from(y, (v, i) => v - mu[i]) };
}

// ---------------------------------------------------------------- autocorrelación
export function acf(e, maxLag) {
  const n = e.length, m = sum(e) / n, c0 = sum(e.map(v => (v - m) ** 2)); const out = [1];
  for (let k = 1; k <= maxLag; k++) { let s = 0; for (let i = k; i < n; i++) s += (e[i] - m) * (e[i - k] - m); out.push(s / c0); }
  return out;
}
export function pacf(r) { // Durbin-Levinson sobre las autocorrelaciones r[0..K]
  const K = r.length - 1, out = [1]; let phi = [];
  for (let k = 1; k <= K; k++) {
    let num = r[k], den = 1; for (let j = 0; j < k - 1; j++) { num -= phi[j] * r[k - 1 - j]; den -= phi[j] * r[j + 1]; }
    const a = num / den, nphi = phi.map((v, j) => v - a * phi[k - 2 - j]); nphi.push(a); phi = nphi; out.push(a);
  }
  return out;
}
export const durbinWatson = e => { let n = 0, d = 0; for (let i = 1; i < e.length; i++) d += (e[i] - e[i - 1]) ** 2; for (const v of e) n += v * v; return d / n; };
export function ljungBox(e, lag) { const n = e.length, r = acf(e, lag); let q = 0; for (let k = 1; k <= lag; k++) q += r[k] ** 2 / (n - k); q *= n * (n + 2); return { q, lag, p: pchisqU(q, lag) }; }
const diagnose = e => { const n = e.length, K = Math.max(2, Math.min(24, Math.floor(n / 3))), r = acf(e, K), lag = Math.max(1, Math.min(10, Math.floor(n / 5))); return { acf: r, pacf: pacf(r), dw: durbinWatson(e), lb: ljungBox(e, lag), band: 1.96 / Math.sqrt(n), n }; };

// ---------------------------------------------------------------- ajuste de la serie agregada
// ser: {y: [..], n: [..]|null} con NaN permitidos. o: {breaks, lag, effects, season, family, ac, nwL, level, scale, labels}
export function fitITS(ser, o = {}) {
  const family = o.family || 'gaussian', ac = o.ac || 'nw', level = o.level || 0.95, N = ser.y.length, D = buildDesign(o), p = D.p;
  const idx = []; for (let t = 0; t < N; t++) { const yv = ser.y[t]; if (!isNum(yv) || D.inLag(t)) continue; if (family !== 'gaussian' && !(isNum(ser.n?.[t]) && ser.n[t] > 0)) continue; idx.push(t); }
  if (idx.length <= p + 1) throw new Error(`Hay ${idx.length} periodos utilizables para ${p} parámetros: reduzca interrupciones, estacionalidad o rezago.`);
  // cada segmento necesita al menos 2 periodos si estima pendiente y 1 si solo nivel
  D.breaks.forEach((T, j) => { const cnt = idx.filter(t => D.seg(t) === j + 1).length, need = D.effects[j] === 'level' ? 1 : D.effects[j] === 'slope' ? 1 : 2; if (cnt < need) throw new Error(`El segmento tras «${D.labels[j]}» tiene ${cnt} periodos utilizables: no alcanza para estimar el efecto.`); });
  const X = idx.map(t => D.row(t)), y = Float64Array.from(idx, t => ser.y[t]), n = idx.length, L = o.nwL || nwLag(n);
  let beta, V, df, rho = null, resid, fitted, core = {};
  if (family === 'gaussian') {
    const f0 = ols(X, y), sig2 = sum(f0.res.map(r => r * r)) / (n - p); df = n - p; resid = f0.res; core.diag = diagnose(Array.from(f0.res));
    if (ac === 'ar1') { const pw = praisWinsten(X, y, idx); beta = pw.beta; V = pw.V; rho = pw.rho; resid = Float64Array.from(y, (v, i) => v - dot(X[i], beta)); core.innov = pw.innov; core.diagFinal = diagnose(Array.from(pw.innov)); }
    else if (ac === 'nw') { beta = f0.beta; V = hac(f0.Hi, X.map((x, i) => Float64Array.from(x, v => v * f0.res[i])), idx, L); V = V.map(r => r.map(v => v * n / (n - p))); }
    else { beta = f0.beta; V = f0.Hi.map(r => r.map(v => v * sig2)); }
    fitted = Float64Array.from(X, x => dot(x, beta)); core.sigma = Math.sqrt(sig2);
    const ym = sum(y) / n, tss = sum(y.map(v => (v - ym) ** 2)), rss = sum(f0.res.map(r => r * r)); core.r2 = 1 - rss / tss; core.aic = n * Math.log(rss / n) + 2 * (p + 1) + n * (1 + Math.log(2 * Math.PI));
  } else {
    const m = Float64Array.from(idx, t => ser.n[t]), off = family === 'poisson' ? Float64Array.from(m, v => Math.log(v)) : null, g = irls(X, y, family, m, off); df = Infinity; beta = g.beta;
    const S = X.map((x, i) => Float64Array.from(x, v => v * g.resp[i]));
    if (ac === 'none') V = g.Hi.map(r => r.map(v => v * Math.max(1, g.phi))); else V = hac(g.Hi, S, idx, L).map(r => r.map(v => v * n / (n - p)));
    resid = g.pearson; fitted = g.mu; core.diag = diagnose(Array.from(g.pearson)); core.phi = g.phi; core.conv = g.conv;
    core.dev = 2 * sum(Array.from(y, (v, i) => (family === 'poisson' ? (v > 0 ? v * Math.log(v / g.mu[i]) : 0) - (v - g.mu[i]) : (v > 0 ? v * Math.log(v / g.mu[i]) : 0) + (m[i] - v > 0 ? (m[i] - v) * Math.log((m[i] - v) / (m[i] - g.mu[i])) : 0))));
  }
  const q = isFinite(df) ? qt(1 - (1 - level) / 2, df) : qnorm(1 - (1 - level) / 2);
  const res = { family, ac, level, D, beta, V, df, q, idx, n, N, rho, resid, fitted, L: ac === 'nw' ? L : null, scale: o.scale ?? (family === 'binomial' ? 100 : 1), ...core };
  return summarize(res);
}

// ---------------------------------------------------------------- resultados a partir de beta, V y el diseño
const LINK = { gaussian: { f: e => e, d: () => 1 }, poisson: { f: e => Math.exp(e), d: e => Math.exp(e) }, binomial: { f: e => expit(e), d: e => expit(e) * (1 - expit(e)) } };
export function summarize(R) {
  const { D, beta, V, df, q, family } = R, sc = R.scale ?? 1, lk = LINK[family], se = beta.map((_, a) => Math.sqrt(V[a][a])), z = beta.map((b, a) => b / se[a]), pv = z.map(v => (isFinite(df) ? p2t(v, df) : p2z(v)));
  const link = family !== 'gaussian'; R.coef = D.names.map((name, a) => ({ name, est: beta[a], se: se[a], lo: beta[a] - q * se[a], hi: beta[a] + q * se[a], p: pv[a], exp: link ? { est: Math.exp(beta[a]), lo: Math.exp(beta[a] - q * se[a]), hi: Math.exp(beta[a] + q * se[a]) } : null }));
  const quad = (g, h = g) => { let s = 0; for (let a = 0; a < g.length; a++) for (let c = 0; c < h.length; c++) s += g[a] * V[a][c] * h[c]; return s; };
  // predicción puntual y con IC de f(η) (escala de la respuesta)
  const point = (t, opts) => { const u = D.row(t, opts), e = dot(u, beta), s = Math.sqrt(Math.max(0, quad(u))); return { est: sc * lk.f(e), lo: sc * lk.f(e - q * s), hi: sc * lk.f(e + q * s) }; };
  // contraste observado − contrafactual (periodos ts) para la interrupción j: diferencia absoluta y cambio relativo
  const contrast = (j, ts) => {
    const gA = new Float64Array(D.p), gB = new Float64Array(D.p); let A = 0, B = 0;
    for (const t of ts) { const u = D.row(t), v = D.row(t, { cf: j }), eu = dot(u, beta), ev = dot(v, beta); A += lk.f(eu); B += lk.f(ev); for (let a = 0; a < D.p; a++) { gA[a] += lk.d(eu) * u[a]; gB[a] += lk.d(ev) * v[a]; } }
    const k = ts.length; A /= k; B /= k; for (let a = 0; a < D.p; a++) { gA[a] /= k; gB[a] /= k; }
    const gd = gA.map((v, a) => v - gB[a]), sd = Math.sqrt(Math.max(0, quad(gd))), gr = gA.map((v, a) => (v * B - A * gB[a]) / (B * B)), sr = Math.sqrt(Math.max(0, quad(gr)));
    return { obs: sc * A, cf: sc * B, abs: { est: sc * (A - B), lo: sc * (A - B - q * sd), hi: sc * (A - B + q * sd) }, rel: { est: A / B - 1, lo: A / B - 1 - q * sr, hi: A / B - 1 + q * sr }, p: (() => { const zz = (A - B) / sd; return isFinite(zz) ? (isFinite(df) ? p2t(zz, df) : p2z(zz)) : NaN; })() };
  };
  R.effects = D.breaks.map((T, j) => {
    const next = D.breaks[j + 1] ?? R.N, t0 = Math.min(T + D.lag, next - 1), win = []; for (let t = t0; t < next; t++) win.push(t);
    return { j, label: D.labels[j], T, t0, t1: next - 1, nWin: win.length, first: contrast(j, [t0]), mean: contrast(j, win), last: contrast(j, [next - 1]), level: D.kind.findIndex(d => d.k === 'step' && d.j === j), slope: D.kind.findIndex(d => d.k === 'slope' && d.j === j) };
  });
  // curvas para las figuras: línea ajustada por segmento y contrafactual de cada interrupción
  const trend = R.display === 'full' ? {} : { noSeason: true };
  R.curves = { seg: [], cf: [] };
  for (let s = 0; s <= D.breaks.length; s++) {
    const a = s === 0 ? 0 : D.breaks[s - 1] + D.lag, b = (D.breaks[s] ?? R.N) - 1; if (a > b) continue; const pts = []; for (let t = a; t <= b; t++) pts.push({ t, ...point(t, trend) }); R.curves.seg.push({ s, pts });
  }
  D.breaks.forEach((T, j) => { const next = D.breaks[j + 1] ?? R.N, pts = []; for (let t = T; t < next; t++) pts.push({ t, ...point(t, { ...trend, cf: j }) }); R.curves.cf.push({ j, pts }); });
  R.point = point; R.contrast = contrast; return R;
}

// ---------------------------------------------------------------- datos individuales
// obs: [{t, y, id}] con t = índice del periodo (o tiempo relativo discreto). Gaussiano → modelo mixto con intercepto aleatorio; Poisson / binario → GEE intercambiable.
export function fitITSIndividual(obs, o = {}) {
  const family = o.family || 'gaussian', level = o.level || 0.95, D = buildDesign(o), rows = obs.filter(r => isNum(r.y) && isNum(r.t) && !D.inLag(r.t));
  const X = rows.map(r => D.row(r.t)), y = Float64Array.from(rows, r => r.y), ids = rows.map(r => r.id), nclust = new Set(ids).size, N = o.N ?? (Math.max(...obs.map(r => r.t)) + 1);
  if (rows.length <= D.p + 2) throw new Error('Muy pocas observaciones para el modelo.');
  let f, extra = {};
  if (family === 'gaussian') { f = fitLMM(X, y, ids); extra = { icc: f.icc, tau2: f.tau2, sigma2: f.sigma2, model: 'Modelo lineal mixto (intercepto aleatorio por persona)' }; }
  else { f = fitGEE(X, y, ids, family, 'exchangeable'); extra = { alpha: f.alpha, model: 'GEE (correlación intercambiable, covarianza robusta)' }; }
  const R = { family, level, D, beta: Float64Array.from(f.beta), V: f.V, df: Infinity, q: qnorm(1 - (1 - level) / 2), n: rows.length, N, nclust, scale: o.scale ?? (family === 'binomial' ? 100 : 1), display: o.display, individual: true, ...extra };
  return summarize(R);
}
// resumen por periodo de los datos individuales (para dibujar los puntos): media, error estándar y n
export function aggregateByPeriod(obs, N, family = 'gaussian', scale = 1) {
  const g = Array.from({ length: N }, () => ({ s: 0, ss: 0, n: 0 }));
  for (const r of obs) if (isNum(r.y) && r.t >= 0 && r.t < N) { const a = g[r.t]; a.s += r.y; a.ss += r.y * r.y; a.n++; }
  return { y: g.map(a => (a.n ? (family === 'gaussian' ? a.s / a.n : a.s / a.n * scale) : NaN)), se: g.map(a => (a.n > 1 ? Math.sqrt(Math.max(0, a.ss / a.n - (a.s / a.n) ** 2) / (a.n - 1)) * (family === 'gaussian' ? 1 : scale) : NaN)), n: g.map(a => a.n) };
}

// meta-análisis de efectos entre grupos (DerSimonian-Laird) para el forest de paneles
export function poolEffects(est, se, level = 0.95) {
  const k = est.length, w = se.map(s => 1 / (s * s)), sw = sum(w), mu0 = sum(est.map((e, i) => w[i] * e)) / sw, Q = sum(est.map((e, i) => w[i] * (e - mu0) ** 2));
  const tau2 = Math.max(0, (Q - (k - 1)) / (sw - sum(w.map(v => v * v)) / sw)), wr = se.map(s => 1 / (s * s + tau2)), mu = sum(est.map((e, i) => wr[i] * e)) / sum(wr), sm = Math.sqrt(1 / sum(wr)), zq = qnorm(1 - (1 - level) / 2);
  return { est: mu, se: sm, lo: mu - zq * sm, hi: mu + zq * sm, tau2, Q, I2: Q > 0 ? Math.max(0, (Q - (k - 1)) / Q) : 0, k };
}
