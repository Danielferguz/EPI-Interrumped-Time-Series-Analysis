// Pronóstico del contrafactual: qué habría ocurrido sin la interrupción, a partir del periodo previo.
// Métodos: tendencia lineal con estacionalidad de Fourier (MCO), suavizado exponencial aditivo de Holt-Winters (con tendencia amortiguada opcional),
// ARIMA(p,d,q) con estacionalidad de Fourier como regresores (suma de cuadrados condicional, selección por AICc), y referencias simples (media, deriva, estacional ingenuo).
// La comparación entre métodos usa validación de origen móvil en el periodo previo. Sin dependencias del navegador.
import { symInv, mv, dot, zeros } from './linalg.js?v=2026.10.25-1';
import { qnorm, qt } from './dist.js?v=2026.10.25-1';

const sum = a => { let s = 0; for (const v of a) s += v; return s; }, mean = a => sum(a) / a.length, expit = z => 1 / (1 + Math.exp(-z));

// ---------------------------------------------------------------- optimizador Nelder-Mead
export function nelderMead(f, x0, { step = 0.3, maxit = 600, tol = 1e-9 } = {}) {
  const n = x0.length; let S = [x0.slice(), ...x0.map((_, i) => { const x = x0.slice(); x[i] += step; return x; })], F = S.map(f);
  for (let it = 0; it < maxit; it++) {
    const o = S.map((_, i) => i).sort((a, b) => F[a] - F[b]); S = o.map(i => S[i]); F = o.map(i => F[i]);
    if (Math.abs(F[n] - F[0]) < tol * (1 + Math.abs(F[0]))) break;
    const c = Array.from({ length: n }, (_, j) => mean(S.slice(0, n).map(x => x[j]))), at = (k) => c.map((v, j) => v + k * (S[n][j] - v)), xr = at(-1), fr = f(xr);
    if (fr < F[0]) { const xe = at(-2), fe = f(xe); if (fe < fr) { S[n] = xe; F[n] = fe; } else { S[n] = xr; F[n] = fr; } }
    else if (fr < F[n - 1]) { S[n] = xr; F[n] = fr; }
    else { const xc = fr < F[n] ? at(-0.5) : at(0.5), fc = f(xc); if (fc < Math.min(fr, F[n])) { S[n] = xc; F[n] = fc; } else for (let i = 1; i <= n; i++) { S[i] = S[i].map((v, j) => S[0][j] + 0.5 * (v - S[0][j])); F[i] = f(S[i]); } }
  }
  const b = F.indexOf(Math.min(...F)); return { x: S[b], f: F[b] };
}

// rellena huecos por interpolación lineal (los métodos de suavizado y ARIMA necesitan una serie completa)
export function fillGaps(y) {
  const out = Array.from(y), n = out.length, ok = out.map(v => Number.isFinite(v)); if (ok.filter(Boolean).length < n * 0.8) throw new Error('Demasiados periodos sin dato en el segmento de entrenamiento (>20 %).');
  for (let i = 0; i < n; i++) if (!ok[i]) { let a = i - 1; while (a >= 0 && !ok[a]) a--; let b = i + 1; while (b < n && !ok[b]) b++; out[i] = a < 0 ? out[b] : b >= n ? out[a] : out[a] + (out[b] - out[a]) * (i - a) / (b - a); ok[i] = true; }
  return out;
}
const fourier = (t, P, K) => { const r = []; for (let k = 1; k <= K; k++) { if (2 * k === P) r.push(Math.cos(Math.PI * t)); else r.push(Math.sin(2 * Math.PI * k * t / P), Math.cos(2 * Math.PI * k * t / P)); } return r; };
const nK = (P, K) => { let c = 0; for (let k = 1; k <= K; k++) c += 2 * k === P ? 1 : 2; return c; };
const seasOK = (P, n, K) => P > 1 && K > 0 && n >= 2 * P;

// ---------------------------------------------------------------- tendencia lineal + Fourier
export function linearForecast(y0, { P = 0, K = 2, level = 0.95, trend = true } = {}) {
  const y = fillGaps(y0), n = y.length, useS = seasOK(P, n, K), row = t => [1, ...(trend ? [t] : []), ...(useS ? fourier(t, P, K) : [])].map(Number), X = y.map((_, t) => Float64Array.from(row(t))), p = X[0].length;
  const H = zeros(p, p), b = new Float64Array(p); for (let i = 0; i < n; i++) for (let a = 0; a < p; a++) { b[a] += X[i][a] * y[i]; for (let c = 0; c < p; c++) H[a][c] += X[i][a] * X[i][c]; }
  const Hi = symInv(H); if (!Hi) throw new Error('Tendencia lineal: segmento de entrenamiento demasiado corto.'); const beta = mv(Hi, b), res = y.map((v, i) => v - dot(X[i], beta)), s2 = sum(res.map(r => r * r)) / Math.max(1, n - p), q = qt(1 - (1 - level) / 2, Math.max(1, n - p));
  return { name: `Tendencia lineal${useS ? ` + Fourier (K=${K})` : ''}`, fitted: y.map((v, i) => v - res[i]), resid: res, sigma: Math.sqrt(s2), k: p, n, useS,
    predict(h) { const fc = [], se = []; for (let s = 0; s < h; s++) { const x = Float64Array.from(row(n + s)); let v = 0; for (let a = 0; a < p; a++) for (let c = 0; c < p; c++) v += x[a] * Hi[a][c] * x[c]; fc.push(dot(x, beta)); se.push(Math.sqrt(s2 * (1 + v))); } return { fc, se, lo: fc.map((v, i) => v - q * se[i]), hi: fc.map((v, i) => v + q * se[i]) }; } };
}

// ---------------------------------------------------------------- suavizado exponencial aditivo (ETS A,{N|A|Ad},{N|A})
function etsRun(y, P, par, trend, seasonal, init) {
  const { al, be, ga, ph } = par; let l = init.l, b = init.b, s = init.s.slice(), sse = 0; const fit = [], res = [];
  for (let t = 0; t < y.length; t++) {
    const si = seasonal ? s[t % P] : 0, f = l + (trend ? ph * b : 0) + si, e = y[t] - f; fit.push(f); res.push(e); sse += e * e;
    const ln = l + (trend ? ph * b : 0) + al * e; if (trend) b = ph * b + be * e; if (seasonal) s[t % P] = si + ga * e; l = ln;
  }
  return { l, b, s, sse, fit, res };
}
export function etsForecast(y0, { P = 0, level = 0.95, trend = 'damped', seasonal = true } = {}) {
  const y = fillGaps(y0), n = y.length, useS = seasonal && P > 1 && n >= 2 * P, useT = trend !== 'none', damped = trend === 'damped';
  const b0 = useT ? (useS ? (mean(y.slice(P, 2 * P)) - mean(y.slice(0, P))) / P : (y[Math.min(n - 1, 4)] - y[0]) / Math.min(n - 1, 4)) : 0, m0 = useS ? mean(y.slice(0, P)) : y[0];
  const init = { l: useS ? m0 - b0 * (P - 1) / 2 : y[0] - b0, b: b0, s: useS ? (() => { const r = y.slice(0, P).map((v, i) => v - m0 - b0 * (i - (P - 1) / 2)), c = mean(r); return r.map(v => v - c); })() : [] };
  const unpack = x => { const al = expit(x[0]), be = useT ? al * expit(x[1]) : 0, ga = useS ? (1 - al) * expit(x[2]) : 0, ph = damped ? 0.8 + 0.195 * expit(x[3]) : 1; return { al, be, ga, ph }; };
  // con los parámetros de suavizado fijos, los errores son lineales en los estados iniciales: se estiman por mínimos cuadrados y se reoptimizan los parámetros
  const refine = par => {
    const zero = Array(n).fill(0), dim = 1 + (useT ? 1 : 0) + (useS ? P : 0), cols = [], e0 = etsRun(y, P, par, useT, useS, { l: 0, b: 0, s: useS ? Array(P).fill(0) : [] }).res;
    for (let i = 0; i < dim; i++) { const st = { l: i === 0 ? 1 : 0, b: useT && i === 1 ? 1 : 0, s: useS ? Array.from({ length: P }, (_, k) => (i === (useT ? 2 : 1) + k ? 1 : 0)) : [] }; cols.push(etsRun(zero, P, par, useT, useS, st).res); }
    // e(y, x0) = e0 + Σ x_i c_i con c_i = respuesta a cada estado inicial: min ‖e0 + C x‖ → x = −(C'C)⁻¹ C' e0
    const H = zeros(dim, dim), g = new Float64Array(dim); for (let i = 0; i < dim; i++) { for (let j = 0; j <= i; j++) { let v = 0; for (let t = 0; t < n; t++) v += cols[i][t] * cols[j][t]; H[i][j] = H[j][i] = v; } let v = 0; for (let t = 0; t < n; t++) v += cols[i][t] * e0[t]; g[i] = -v; H[i][i] += 1e-8; }
    const Hi = symInv(H); if (!Hi) return null; const x = mv(Hi, g), st = { l: x[0], b: useT ? x[1] : 0, s: useS ? x.slice(useT ? 2 : 1, (useT ? 2 : 1) + P) : [] }; if (useS) { const c = mean(st.s); st.s = st.s.map(v => v - c); st.l += c; }
    return st;
  };
  let cur = init, par = null, best = null;
  for (let round = 0; round < 3; round++) {
    const obj = x => { const r = etsRun(y, P, unpack(x), useT, useS, cur).sse; return isFinite(r) ? r : 1e300; };
    best = [[0, 0, 0, 0], [-1.5, -1, -1, 0], [1, 1, -1, 1], [-2, 0, 1, -1], [0.5, -2, 0, 2]].map(x0 => nelderMead(obj, round ? best.x.map((v, i) => v + (x0[i] ? 0.3 * x0[i] : 0)) : x0, { step: 1, maxit: 1500 })).sort((u, v) => u.f - v.f)[0]; par = unpack(best.x);
    const st = refine(par); if (!st) break; cur = st;
  }
  const run = etsRun(y, P, par, useT, useS, cur), kPar = 1 + (useT ? 1 : 0) + (useS ? 1 : 0) + (damped ? 1 : 0), kSt = 1 + (useT ? 1 : 0) + (useS ? P - 1 : 0), k = kPar + kSt;
  const sigma2 = run.sse / Math.max(1, n - k), q = qnorm(1 - (1 - level) / 2), aicc = n * Math.log(run.sse / n) + 2 * k + 2 * k * (k + 1) / Math.max(1, n - k - 1);
  return { name: `Holt-Winters aditivo (${useT ? (damped ? 'tendencia amortiguada' : 'tendencia') : 'sin tendencia'}${useS ? ', estacional' : ''})`, par, fitted: run.fit, resid: run.res, sigma: Math.sqrt(sigma2), aicc, k, n, sse: run.sse,
    predict(h) {
      const fc = [], se = []; let cum = 0, phsum = 0;
      for (let s = 1; s <= h; s++) { phsum = damped ? phsum + Math.pow(par.ph, s) : s; fc.push(run.l + (useT ? phsum * run.b : 0) + (useS ? run.s[(n + s - 1) % P] : 0));
        const j = s - 1; let cj = 0; if (j > 0) { const ph = par.ph; cj = par.al + (useT ? par.be * (damped ? ph * (1 - Math.pow(ph, j)) / (1 - ph) : j) : 0) + (useS && j % P === 0 ? par.ga : 0); }
        cum += cj * cj; se.push(Math.sqrt(sigma2 * (1 + cum))); }
      return { fc, se, lo: fc.map((v, i) => v - q * se[i]), hi: fc.map((v, i) => v + q * se[i]) };
    } };
}
export function etsAuto(y, o = {}) { // elige tendencia (ninguna, aditiva o amortiguada) y estacionalidad por AICc
  const cands = ['none', 'additive', 'damped'].map(tr => { try { return etsForecast(y, { ...o, trend: tr }); } catch { return null; } }).filter(Boolean);
  return cands.sort((a, b) => a.aicc - b.aicc)[0];
}

// ---------------------------------------------------------------- ARIMA(p,d,q) con regresores de Fourier, estimación por suma de cuadrados condicional
const stationary = c => { if (c.length === 0) return true; if (c.length === 1) return Math.abs(c[0]) < 0.995; return Math.abs(c[1]) < 0.995 && c[0] + c[1] < 0.995 && c[1] - c[0] < 0.995; };
function diffN(y, d) { let w = y.slice(); for (let i = 0; i < d; i++) w = w.slice(1).map((v, j) => v - w[j]); return w; }
function cssResiduals(w, phi, theta, c) {
  const n = w.length, p = phi.length, q = theta.length, e = new Float64Array(n); let sse = 0, cnt = 0;
  for (let t = p; t < n; t++) { let v = w[t] - c; for (let i = 0; i < p; i++) v -= phi[i] * w[t - 1 - i]; for (let j = 0; j < q; j++) if (t - 1 - j >= p) v -= theta[j] * e[t - 1 - j]; e[t] = v; sse += v * v; cnt++; }
  return { e, sse, cnt };
}
export function adfStat(y) { // Dickey-Fuller aumentado (1 rezago, constante y tendencia): estadístico t de γ en Δy_t = a + bt + γ y_{t−1} + δ Δy_{t−1}
  const n = y.length; if (n < 12) return NaN; const rows = [], z = [];
  for (let t = 2; t < n; t++) { rows.push([1, t, y[t - 1], y[t - 1] - y[t - 2]]); z.push(y[t] - y[t - 1]); }
  const p = 4, H = zeros(p, p), b = new Float64Array(p); rows.forEach((r, i) => { for (let a = 0; a < p; a++) { b[a] += r[a] * z[i]; for (let c = 0; c < p; c++) H[a][c] += r[a] * r[c]; } }); const Hi = symInv(H); if (!Hi) return NaN;
  const beta = mv(Hi, b), res = z.map((v, i) => v - dot(rows[i], beta)), s2 = sum(res.map(r => r * r)) / (rows.length - p); return beta[2] / Math.sqrt(s2 * Hi[2][2]);
}
function fitArima(w, p, q, withC) {
  const unpack = x => ({ phi: x.slice(0, p), theta: x.slice(p, p + q), c: withC ? x[p + q] : 0 });
  const obj = x => { const { phi, theta, c } = unpack(x); if (!stationary(phi) || !stationary(theta)) return 1e300; const r = cssResiduals(w, phi, theta, c); return r.sse; };
  const m0 = mean(w), x0 = [...Array(p + q).fill(0), ...(withC ? [m0] : [])], r = nelderMead(obj, x0, { step: 0.2, maxit: 1200 }); return { ...unpack(r.x), ...cssResiduals(w, ...(({ phi, theta, c }) => [phi, theta, c])(unpack(r.x))) };
}
export function arimaForecast(y0, { P = 0, K = 2, level = 0.95, order = 'auto', d = 'auto', drift = true } = {}) {
  const y = fillGaps(y0), n = y.length, useS = seasOK(P, n, K), q0 = qnorm(1 - (1 - level) / 2);
  // 1) parte estacional por regresión de Fourier (sin tendencia); 2) ARIMA sobre el resto
  let sBeta = null, base = y.slice(); const sRow = t => [1, ...fourier(t, P, K)];
  if (useS) { const X = y.map((_, t) => Float64Array.from(sRow(t))), pp = X[0].length, H = zeros(pp, pp), b = new Float64Array(pp); X.forEach((x, i) => { for (let a = 0; a < pp; a++) { b[a] += x[a] * y[i]; for (let c = 0; c < pp; c++) H[a][c] += x[a] * x[c]; } }); sBeta = mv(symInv(H), b); base = y.map((v, i) => v - dot(X[i], sBeta) + sBeta[0]); }
  const dd = d === 'auto' ? (adfStat(base) < -3.41 ? 0 : 1) : d, w = diffN(base, dd), withC = dd === 0 || drift;
  const grid = order === 'auto' ? [0, 1, 2].flatMap(p => [0, 1, 2].map(q => [p, q])) : [[order[0], order[2]]];
  let best = null; for (const [p, q] of grid) { if (w.length < p + q + 6) continue; let f; try { f = fitArima(w, p, q, withC); } catch { continue; } const m = f.cnt, k = p + q + (withC ? 1 : 0) + 1, aicc = m * Math.log(f.sse / m) + 2 * k + 2 * k * (k + 1) / Math.max(1, m - k - 1); if (!best || aicc < best.aicc) best = { ...f, p, q, aicc, k }; }
  if (!best) throw new Error('ARIMA: segmento de entrenamiento demasiado corto.');
  const { phi, theta, c, e, p, q } = best, sigma2 = best.sse / best.cnt, Phi = (() => { let poly = [1, ...phi.map(v => -v)]; for (let i = 0; i < dd; i++) { const nx = Array(poly.length + 1).fill(0); poly.forEach((v, j) => { nx[j] += v; nx[j + 1] -= v; }); poly = nx; } return poly.slice(1).map(v => -v); })();
  // residuos de una etapa de la serie original (para ACF y precisión dentro de la muestra)
  const off = n - w.length, resid = Array(n).fill(NaN); for (let t = 0; t < w.length; t++) resid[t + off] = e[t];
  return { name: `ARIMA(${p},${dd},${q})${withC && dd > 0 ? ' con deriva' : ''}${useS ? ` + Fourier (K=${K})` : ''}`, order: [p, dd, q], phi, theta, c, sigma: Math.sqrt(sigma2), aicc: best.aicc, resid, fitted: y.map((v, i) => v - (resid[i] || 0)), n, k: best.k,
    predict(h) {
      const yy = base.slice(), ee = Array.from(e), off2 = n - e.length, fcb = []; const eAt = i => (i - off2 >= 0 && i - off2 < ee.length ? ee[i - off2] : 0);
      for (let s = 0; s < h; s++) { const t = n + s; let v = c; for (let i = 0; i < Phi.length; i++) v += Phi[i] * (yy[t - 1 - i] ?? 0); for (let j = 0; j < q; j++) v += theta[j] * (t - 1 - j < n ? eAt(t - 1 - j) : 0); yy.push(v); fcb.push(v); }
      const psi = [1]; for (let j = 1; j < h; j++) { let v = j <= q ? theta[j - 1] : 0; for (let i = 1; i <= Math.min(j, Phi.length); i++) v += Phi[i - 1] * psi[j - i]; psi.push(v); }
      const fc = fcb.map((v, s) => v + (useS ? dot(Float64Array.from(sRow(n + s)), sBeta) - sBeta[0] : 0)); let cum = 0; const se = psi.map(v => { cum += v * v; return Math.sqrt(sigma2 * cum); });
      return { fc, se, lo: fc.map((v, i) => v - q0 * se[i]), hi: fc.map((v, i) => v + q0 * se[i]) };
    } };
}

// ---------------------------------------------------------------- referencias simples
export function simpleForecast(y0, kind, { P = 0, level = 0.95 } = {}) {
  const y = fillGaps(y0), n = y.length, q = qnorm(1 - (1 - level) / 2); let name, fcf, res;
  if (kind === 'mean') { const m = mean(y); name = 'Media del periodo previo'; res = y.map(v => v - m); fcf = () => m; }
  else if (kind === 'drift') { const sl = (y[n - 1] - y[0]) / (n - 1); name = 'Deriva (pendiente primero-último)'; res = y.slice(1).map((v, i) => v - y[i] - sl); fcf = h => y[n - 1] + sl * h; }
  else { const m = P > 1 && n > P ? P : 1; name = m > 1 ? 'Estacional ingenuo (mismo periodo del ciclo anterior)' : 'Ingenuo (último valor)'; res = y.slice(m).map((v, i) => v - y[i]); fcf = h => y[n - m + ((h - 1) % m)]; }
  const s = Math.sqrt(sum(res.map(r => r * r)) / Math.max(1, res.length - 1)), growth = kind === 'mean' ? () => 1 : kind === 'drift' ? h => Math.sqrt(h * (1 + h / (n - 1))) : h => Math.sqrt(Math.ceil(h / (P > 1 ? P : 1)));
  return { name, sigma: s, resid: res, fitted: y.map((v, i) => v - (res[i - (n - res.length)] || 0)), n, k: 1, predict(h) { const fc = [], se = []; for (let i = 1; i <= h; i++) { fc.push(fcf(i)); se.push(s * growth(i)); } return { fc, se, lo: fc.map((v, i) => v - q * se[i]), hi: fc.map((v, i) => v + q * se[i]) }; } };
}

export const METHOD_EN = { linear: 'Linear trend + seasonality', ets: 'Holt-Winters (additive)', arima: 'ARIMA (automatic, AICc)', drift: 'Drift', snaive: 'Seasonal naive', mean: 'Pre-period mean' };
export const METHODS = [
  ['linear', 'Tendencia lineal + estacionalidad (MCO)'], ['ets', 'Holt-Winters aditivo (suavizado exponencial)'], ['arima', 'ARIMA automático (AICc)'], ['drift', 'Deriva'], ['snaive', 'Estacional ingenuo'], ['mean', 'Media del periodo previo'],
];
// ajusta el método elegido sobre y (segmento de entrenamiento) y devuelve el modelo con predict(h)
export function fitMethod(y, method, o = {}) {
  if (method === 'linear') return linearForecast(y, o); if (method === 'ets') return etsAuto(y, o); if (method === 'arima') return arimaForecast(y, o);
  return simpleForecast(y, method, o);
}

// ---------------------------------------------------------------- validación de origen móvil dentro del periodo previo
// Para cada origen, se ajusta con los datos hasta ahí y se pronostica h periodos; devuelve RMSE, MAE, sesgo y cobertura del intervalo.
export function rollingOrigin(y0, method, { h = 6, minTrain = null, maxOrigins = 8, ...o } = {}) {
  const y = Array.from(y0), n = y.length, start = minTrain ?? Math.max(8, Math.floor(n * 0.5)), last = n - h; if (last < start) return { n: 0, rmse: NaN, mae: NaN, bias: NaN, cover: NaN };
  const origins = []; const stepO = Math.max(1, Math.floor((last - start) / (maxOrigins - 1 || 1))); for (let s = start; s <= last; s += stepO) origins.push(s);
  const err = [], cov = []; for (const s of origins) { let m; try { m = fitMethod(y.slice(0, s), method, o); } catch { continue; } const f = m.predict(h); for (let i = 0; i < h; i++) { const a = y[s + i]; if (!Number.isFinite(a)) continue; err.push(a - f.fc[i]); cov.push(a >= f.lo[i] && a <= f.hi[i] ? 1 : 0); } }
  if (!err.length) return { n: 0, rmse: NaN, mae: NaN, bias: NaN, cover: NaN };
  return { n: err.length, origins: origins.length, rmse: Math.sqrt(mean(err.map(e => e * e))), mae: mean(err.map(Math.abs)), bias: mean(err), cover: mean(cov) };
}

// ---------------------------------------------------------------- efecto de cada interrupción respecto al pronóstico
// obs: serie completa; breaks/lag como en fitITS; spec: {method, P, K, level, ...}; train: 'previous' (solo el segmento previo) | 'all' (todo lo anterior)
export function counterfactualEffects(obs, { breaks, lag = 0, spec, train = 'previous', minTrain = 10 }) {
  const B = [...breaks].sort((a, b) => a - b), N = obs.length, level = spec.level || 0.95, z = qnorm(1 - (1 - level) / 2);
  return B.map((T, j) => {
    const next = B[j + 1] ?? N, prevStart = j === 0 ? 0 : B[j - 1] + (train === 'previous' ? lag : 0), from = train === 'previous' && T - prevStart >= minTrain ? prevStart : 0, trainY = obs.slice(from, T), h = next - T;
    const m = fitMethod(trainY, spec.method, spec), f = m.predict(h), t0 = T + lag, idx = []; for (let t = t0; t < next; t++) if (Number.isFinite(obs[t])) idx.push(t);
    const at = t => ({ obs: obs[t], cf: f.fc[t - T], lo: f.lo[t - T], hi: f.hi[t - T], se: f.se[t - T] });
    const A = idx.length ? mean(idx.map(t => obs[t])) : NaN, Bc = idx.length ? mean(idx.map(t => f.fc[t - T])) : NaN, seM = idx.length ? mean(idx.map(t => f.se[t - T])) : NaN; // IC conservador: errores de pronóstico perfectamente correlacionados
    const win = { obs: A, cf: Bc, abs: { est: A - Bc, lo: A - Bc - z * seM, hi: A - Bc + z * seM }, rel: { est: A / Bc - 1, lo: A / (Bc + z * seM) - 1, hi: A / (Bc - z * seM) - 1 } };
    const tl = idx[idx.length - 1], one = tl == null ? null : (() => { const r = at(tl); return { ...r, abs: { est: r.obs - r.cf, lo: r.obs - r.hi, hi: r.obs - r.lo }, rel: { est: r.obs / r.cf - 1, lo: r.obs / r.hi - 1, hi: r.obs / r.lo - 1 } }; })();
    return { j, T, t0, next, from, model: m, name: m.name, fc: f, mean: win, last: one, nTrain: trainY.length, t: Array.from({ length: h }, (_, i) => T + i) };
  });
}
