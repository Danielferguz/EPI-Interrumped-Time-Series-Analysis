// js/its.js frente a statsmodels (tests/reference_its.json)
import { readFileSync } from 'node:fs';
import { fitITS, fitITSIndividual, buildDesign, poolEffects } from '../js/its.js';
const R = JSON.parse(readFileSync(new URL('reference_its.json', import.meta.url)));
let fails = 0, checks = 0; const ok = (n, c, x = '') => { checks++; if (!c) { fails++; console.log('FAIL', n, x); } else console.log('ok  ', n); };
const near = (a, b, tol = 1e-6) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), allNear = (A, B, tol) => A.length === B.length && A.every((v, i) => near(v, B[i], tol));
const season = { P: 12, K: 2 }, g = R.gauss, base = { breaks: g.breaks, season, display: 'full' };
{ const f = fitITS({ y: g.y }, { ...base, ac: 'none' });
  ok('gaussiano: coeficientes = statsmodels OLS', allNear(Array.from(f.beta), g.beta, 1e-8), JSON.stringify(Array.from(f.beta).slice(0, 3)));
  ok('gaussiano: EE clásicos = statsmodels', allNear(f.coef.map(c => c.se), g.se_ols, 1e-8));
  ok('gaussiano: R² = statsmodels', near(f.r2, g.r2, 1e-9));
  ok('gaussiano: Durbin-Watson', near(f.diag.dw, g.dw, 1e-9));
  ok('gaussiano: Ljung-Box (10 rezagos) estadístico y P', near(f.diag.lb.q, g.lb10, 1e-8) && near(f.diag.lb.p, g.lb10p, 1e-6), `${f.diag.lb.q} ${g.lb10}`);
  ok('ACF de los residuos = statsmodels', allNear(f.diag.acf.slice(0, 13), g.acf, 1e-8));
  ok('PACF (Durbin-Levinson) = statsmodels ywm', allNear(f.diag.pacf.slice(0, 13), g.pacf, 1e-7)); }
{ const f = fitITS({ y: g.y }, { ...base, ac: 'nw', nwL: g.L });
  ok('gaussiano: Newey-West (Bartlett, corrección n/(n−p)) = statsmodels HAC', allNear(f.coef.map(c => c.se), g.se_hac, 1e-7), JSON.stringify(f.coef.map(c => c.se).slice(0, 3)) + ' ' + JSON.stringify(g.se_hac.slice(0, 3))); }
{ const f = fitITS({ y: g.y }, { ...base, ac: 'ar1' }), w = R.pw;
  ok('Prais-Winsten: ρ y coeficientes = implementación independiente', near(f.rho, w.rho, 1e-6) && allNear(Array.from(f.beta), w.beta, 1e-6), `${f.rho} ${w.rho}`);
  ok('Prais-Winsten: errores estándar', allNear(f.coef.map(c => c.se), w.se, 1e-6));
}
{ let a = 12345; const u = () => ((a = (a * 1664525 + 1013904223) >>> 0) / 4294967296), nrm = () => Math.sqrt(-2 * Math.log(1 - u())) * Math.cos(2 * Math.PI * u()); let e = 0; const y = Array.from({ length: 600 }, (_, t) => { e = 0.6 * e + nrm(); return 5 + 0.01 * t - 2 * (t >= 300) + e; });
  const f = fitITS({ y }, { breaks: [300], ac: 'ar1' }); ok('Prais-Winsten recupera ρ = 0.6 en una serie larga', Math.abs(f.rho - 0.6) < 0.08, f.rho); }
{ const p = R.pois, f = fitITS({ y: p.y, n: p.n }, { ...base, family: 'poisson', ac: 'nw', nwL: g.L });
  ok('Poisson con desplazamiento: coeficientes = statsmodels GLM', allNear(Array.from(f.beta), p.beta, 1e-6));
  ok('Poisson: HAC = statsmodels', allNear(f.coef.map(c => c.se), p.se_hac, 1e-5), JSON.stringify(f.coef.map(c => c.se).slice(0, 3)) + ' ' + JSON.stringify(p.se_hac.slice(0, 3)));
  const f0 = fitITS({ y: p.y, n: p.n }, { ...base, family: 'poisson', ac: 'none' }); ok('Poisson: EE modelo × √φ (cuasi-Poisson)', allNear(f0.coef.map(c => c.se), p.se_naive.map(s => s * Math.sqrt(Math.max(1, p.phi))), 1e-5)); }
{ const b = R.binom, f = fitITS({ y: b.y, n: b.n }, { ...base, family: 'binomial', ac: 'nw', nwL: g.L });
  ok('binomial agrupada: coeficientes = statsmodels GLM', allNear(Array.from(f.beta), b.beta, 1e-6));
  ok('binomial: HAC = statsmodels', allNear(f.coef.map(c => c.se), b.se_hac, 1e-5)); }
{ const m = R.lmm, obs = m.t.map((t, i) => ({ t, y: m.y[i], id: m.id[i] })), f = fitITSIndividual(obs, { breaks: m.breaks, family: 'gaussian', N: 36 });
  ok('modelo mixto individual: efectos fijos = statsmodels MixedLM (REML)', allNear(Array.from(f.beta), m.beta, 2e-3), JSON.stringify(Array.from(f.beta)) + ' ' + JSON.stringify(m.beta));
  ok('modelo mixto: EE = statsmodels', allNear(f.coef.map(c => c.se), m.se, 5e-3));
  const e = R.gee, f2 = fitITSIndividual(m.t.map((t, i) => ({ t, y: e.y[i], id: m.id[i] })), { breaks: m.breaks, family: 'poisson', N: 36 });
  ok('GEE Poisson individual: coeficientes y EE robustos = statsmodels', allNear(Array.from(f2.beta), e.beta, 5e-3) && allNear(f2.coef.map(c => c.se), e.se, 2e-2), JSON.stringify(f2.coef.map(c => c.se)) + ' ' + JSON.stringify(e.se)); }
{ const f = fitITS({ y: g.y }, { breaks: g.breaks, ac: 'none' }), e0 = f.effects[0], D = f.D, ib = D.kind.findIndex(d => d.k === 'step' && d.j === 0), isl = D.kind.findIndex(d => d.k === 'slope' && d.j === 0);
  ok('efecto inmediato = cambio de nivel', near(e0.first.abs.est, f.beta[ib], 1e-9));
  ok('efecto al final del segmento = nivel + pendiente × (t − T)', near(e0.last.abs.est, f.beta[ib] + f.beta[isl] * (e0.t1 - e0.T), 1e-9));
  const e1 = f.effects[1], jb = D.kind.findIndex(d => d.k === 'step' && d.j === 1); ok('2.ª interrupción: efecto = su cambio de nivel (relativo al régimen previo)', near(e1.first.abs.est, f.beta[jb], 1e-9));
  ok('IC del efecto inmediato = EE del coeficiente × t', near(e0.first.abs.hi - e0.first.abs.est, f.q * f.coef[ib].se, 1e-8)); }
{ const f = fitITS({ y: g.y }, { breaks: [30], lag: 3, ac: 'none' }); ok('rezago: se excluyen 3 periodos del ajuste', f.n === 72 - 3 && !f.idx.includes(30) && !f.idx.includes(32) && f.idx.includes(33)); }
{ const D = buildDesign({ breaks: [20, 40, 60] }), tr = [10, 0.2, -3, -0.1, -2, 0.05, -1, 0], y = Array.from({ length: 80 }, (_, t) => D.row(t).reduce((s, v, a) => s + v * tr[a], 0) + (((t * 7919) % 13) - 6) * 0.02);
  const f = fitITS({ y }, { breaks: [20, 40, 60], ac: 'none' }); ok('tres interrupciones: recupera niveles y pendientes', allNear(Array.from(f.beta), tr, 0.05), JSON.stringify(Array.from(f.beta).map(v => +v.toFixed(2)))); }
{ const p = poolEffects([-1, -1.2, -0.8, -1.1], [0.2, 0.25, 0.3, 0.2]); ok('meta-análisis de efectos: estimador agrupado razonable', Math.abs(p.est + 1.05) < 0.1 && p.I2 >= 0 && p.hi > p.est); }
console.log(`\n${checks - fails}/${checks} comprobaciones correctas`); process.exit(fails ? 1 : 0);
