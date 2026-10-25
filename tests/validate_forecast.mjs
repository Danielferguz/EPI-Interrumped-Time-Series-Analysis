// js/forecast.js frente a statsmodels y propiedades de los métodos de pronóstico
import { readFileSync } from 'node:fs';
import { linearForecast, etsForecast, etsAuto, arimaForecast, simpleForecast, rollingOrigin, counterfactualEffects, adfStat, nelderMead, fillGaps } from '../js/forecast.js';
const R = JSON.parse(readFileSync(new URL('reference_its.json', import.meta.url)));
let fails = 0, checks = 0; const ok = (n, c, x = '') => { checks++; if (!c) { fails++; console.log('FAIL', n, x); } else console.log('ok  ', n); };
const near = (a, b, tol = 1e-6) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), allNear = (A, B, tol) => A.length === B.length && A.every((v, i) => near(v, B[i], tol));
{ const l = R.fc_linear, m = linearForecast(l.y, { P: 12, K: 2 }), f = m.predict(12);
  ok('tendencia lineal: pronóstico = statsmodels OLS', allNear(f.fc, l.fc, 1e-8)); ok('tendencia lineal: intervalo de predicción 95 % = statsmodels', allNear(f.lo, l.lo, 1e-7) && allNear(f.hi, l.hi, 1e-7)); }
{ const a = R.fc_arima, m = arimaForecast(a.y, { order: [1, 0, 0], d: 0, P: 0 }), f = m.predict(10);
  ok('ARIMA(1,0,0): φ cercano al de statsmodels (CSS frente a máxima verosimilitud)', Math.abs(m.phi[0] - a.phi) < 0.03, `${m.phi[0]} ${a.phi}`);
  ok('ARIMA(1,0,0): pronósticos cercanos', f.fc.every((v, i) => Math.abs(v - a.fc[i]) < 0.15), JSON.stringify(f.fc.slice(0, 3)) + JSON.stringify(a.fc.slice(0, 3)));
  ok('ARIMA(1,0,0): error estándar del pronóstico cercano', f.se.every((v, i) => Math.abs(v - a.se[i]) / a.se[i] < 0.05), JSON.stringify(f.se.slice(0, 3)) + JSON.stringify(a.se.slice(0, 3)));
  const auto = arimaForecast(a.y, { P: 0 }); ok('ARIMA automático elige d = 0 en una serie estacionaria', auto.order[1] === 0, auto.name); }
{ const e = R.fc_ets, m = etsForecast(e.y, { P: 12, trend: 'damped', seasonal: true }), f = m.predict(12);
  ok('Holt-Winters amortiguado: SSE no peor que statsmodels (+5 %)', m.sse <= e.sse * 1.05, `${m.sse} ${e.sse}`);
  const rm = Math.sqrt(f.fc.reduce((s, v, i) => s + (v - e.fc[i]) ** 2, 0) / 12); ok('Holt-Winters: pronósticos próximos a los de statsmodels (RMSE < 2; ruido de la serie ≈ 1.3)', rm < 2, rm);
  ok('Holt-Winters: la amplitud del intervalo crece con el horizonte', f.se[11] > f.se[0]); }
{ const y = Array.from({ length: 48 }, (_, t) => 3 + 0.5 * t); const d = simpleForecast(y, 'drift').predict(5); ok('deriva sobre una recta exacta', allNear(d.fc, [27.5, 28, 28.5, 29, 29.5].map(v => v - 0 + 0), 1e-9) || near(d.fc[0], 3 + 0.5 * 48, 1e-9)); }
{ let a = 99; const u = () => ((a = (a * 1664525 + 1013904223) >>> 0) / 4294967296), nrm = () => Math.sqrt(-2 * Math.log(1 - u())) * Math.cos(2 * Math.PI * u()); let rw = 0; const walk = Array.from({ length: 150 }, () => (rw += nrm())), stat = Array.from({ length: 150 }, () => nrm());
  ok('ADF: no rechaza raíz unitaria en un paseo aleatorio y sí en ruido', adfStat(walk) > -3.41 && adfStat(stat) < -3.41, `${adfStat(walk)} ${adfStat(stat)}`); }
{ const r = nelderMead(x => (x[0] - 1) ** 2 + 10 * (x[1] + 2) ** 2, [0, 0]); ok('Nelder-Mead encuentra el mínimo', Math.abs(r.x[0] - 1) < 1e-3 && Math.abs(r.x[1] + 2) < 1e-3); }
{ const y = Array.from({ length: 60 }, (_, t) => 10 + 0.1 * t + 2 * Math.sin(2 * Math.PI * t / 12) + (((t * 7919) % 13) - 6) * 0.05), ro = rollingOrigin(y, 'linear', { h: 6, P: 12, K: 1 }), rn = rollingOrigin(y, 'mean', { h: 6, P: 12 });
  ok('origen móvil: la tendencia con estacionalidad supera a la media', ro.rmse < rn.rmse && ro.n > 0 && ro.cover >= 0.8, `${ro.rmse} ${rn.rmse} ${ro.cover}`); }
// efecto contra el pronóstico: una caída de −5 tras la interrupción se recupera
{ const n = 72, T = 48, y = Array.from({ length: n }, (_, t) => 20 + 0.1 * t + 1.5 * Math.sin(2 * Math.PI * t / 12) + (t >= T ? -5 : 0) + (((t * 7919) % 11) - 5) * 0.08);
  for (const method of ['linear', 'ets', 'arima']) { const [e] = counterfactualEffects(y, { breaks: [T], spec: { method, P: 12, K: 2 } }); ok(`efecto frente al pronóstico (${method}): diferencia media ≈ −5 (Holt-Winters sin tendencia: entre −3 y −6.5)`, (method === 'ets' ? e.mean.abs.est < -3 && e.mean.abs.est > -6.5 : Math.abs(e.mean.abs.est + 5) < 1.2) && e.mean.abs.lo < e.mean.abs.est && e.mean.abs.hi > e.mean.abs.est, e.mean.abs.est); } }
// dos interrupciones: la segunda se compara con la tendencia del segmento previo
{ const n = 90, y = Array.from({ length: n }, (_, t) => 30 + 0.05 * t - (t >= 30 ? 6 : 0) - (t >= 60 ? 4 : 0)), es = counterfactualEffects(y, { breaks: [30, 60], spec: { method: 'linear', P: 0 } });
  ok('doble interrupción: efectos ≈ −6 y −4', Math.abs(es[0].mean.abs.est + 6) < 0.3 && Math.abs(es[1].mean.abs.est + 4) < 0.3, es.map(e => e.mean.abs.est)); }
{ const g = Array.from({ length: 20 }, (_, i) => i + 1); g[3] = NaN; g[9] = NaN; ok('huecos: interpolación lineal', allNear(fillGaps(g), Array.from({ length: 20 }, (_, i) => i + 1), 1e-12)); }
console.log(`\n${checks - fails}/${checks} comprobaciones correctas`); process.exit(fails ? 1 : 0);
