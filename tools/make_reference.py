# Genera tests/reference_its.json: series simuladas y los resultados de statsmodels con los que se comprueba js/its.js.
# Uso: python3 tools/make_reference.py   (requiere numpy, pandas y statsmodels)
import json, numpy as np, pandas as pd, statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import durbin_watson
from statsmodels.tsa.stattools import acf, pacf
rng = np.random.default_rng(7); out = {}
N = 72; t = np.arange(N); breaks = [30, 54]
def design(t, breaks, season=True, lag=0):
    cols = [np.ones_like(t, dtype=float), t.astype(float)]
    for T in breaks: cols += [(t >= T).astype(float), np.where(t >= T, t - T, 0).astype(float)]
    if season:
        for k in (1, 2): cols += [np.sin(2 * np.pi * k * t / 12), np.cos(2 * np.pi * k * t / 12)]
    return np.column_stack(cols)
X = design(t, breaks)
# serie gaussiana con AR(1) y estacionalidad
e = np.zeros(N); 
for i in range(1, N): e[i] = 0.5 * e[i - 1] + rng.normal(0, 1)
beta = np.array([20, 0.10, -3.0, -0.12, -2.0, 0.05, 1.5, 0.5, 0.6, -0.3]); yg = X @ beta + e
m = sm.OLS(yg, X).fit(); L = int(np.floor(4 * (N / 100) ** (2 / 9)))
mh = sm.OLS(yg, X).fit(cov_type='HAC', cov_kwds={'maxlags': L, 'use_correction': True})
out['gauss'] = {'t': t.tolist(), 'y': yg.tolist(), 'breaks': breaks, 'L': L, 'beta': m.params.tolist(), 'se_ols': m.bse.tolist(), 'se_hac': mh.bse.tolist(), 'dw': float(durbin_watson(m.resid)),
                'lb10': float(acorr_ljungbox(m.resid, lags=[10], return_df=True)['lb_stat'].iloc[0]), 'lb10p': float(acorr_ljungbox(m.resid, lags=[10], return_df=True)['lb_pvalue'].iloc[0]),
                'acf': acf(m.resid, nlags=12, fft=False).tolist(), 'pacf': pacf(m.resid, nlags=12, method='ywm').tolist(), 'r2': float(m.rsquared)}
# Prais-Winsten independiente (numpy)
def pw(X, y):
    rho = 0.0
    for _ in range(200):
        Xs = np.vstack([X[0] * np.sqrt(1 - rho ** 2), X[1:] - rho * X[:-1]]); ys = np.concatenate([[y[0] * np.sqrt(1 - rho ** 2)], y[1:] - rho * y[:-1]])
        b = np.linalg.lstsq(Xs, ys, rcond=None)[0]; r = y - X @ b; nr = (r[1:] * r[:-1]).sum() / (r[:-1] ** 2).sum()
        if abs(nr - rho) < 1e-10: rho = nr; break
        rho = nr
    Xs = np.vstack([X[0] * np.sqrt(1 - rho ** 2), X[1:] - rho * X[:-1]]); ys = np.concatenate([[y[0] * np.sqrt(1 - rho ** 2)], y[1:] - rho * y[:-1]])
    b = np.linalg.lstsq(Xs, ys, rcond=None)[0]; res = ys - Xs @ b; s2 = (res ** 2).sum() / (len(y) - X.shape[1]); V = s2 * np.linalg.inv(Xs.T @ Xs)
    return b, np.sqrt(np.diag(V)), rho
b, s, r = pw(X, yg); out['pw'] = {'beta': b.tolist(), 'se': s.tolist(), 'rho': float(r)}
# Poisson con desplazamiento
pop = rng.integers(8000, 12000, N).astype(float); eta = np.log(pop) + X @ np.array([-5, 0.002, -0.25, -0.004, -0.15, 0.001, 0.10, 0.04, 0.05, -0.02]); yc = rng.poisson(np.exp(eta)).astype(float)
gp = sm.GLM(yc, X, family=sm.families.Poisson(), offset=np.log(pop)).fit(); gh = sm.GLM(yc, X, family=sm.families.Poisson(), offset=np.log(pop)).fit(cov_type='HAC', cov_kwds={'maxlags': L, 'use_correction': True})
out['pois'] = {'y': yc.tolist(), 'n': pop.tolist(), 'beta': gp.params.tolist(), 'se_naive': gp.bse.tolist(), 'se_hac': gh.bse.tolist(), 'phi': float(gp.pearson_chi2 / gp.df_resid)}
# binomial agrupada
trials = rng.integers(300, 600, N).astype(float); pr = 1 / (1 + np.exp(-(X @ np.array([-1.2, 0.004, -0.3, -0.005, -0.2, 0.003, 0.12, 0.05, 0.06, -0.03])))); yb = rng.binomial(trials.astype(int), pr).astype(float)
gb = sm.GLM(np.column_stack([yb, trials - yb]), X, family=sm.families.Binomial()).fit(); gbh = sm.GLM(np.column_stack([yb, trials - yb]), X, family=sm.families.Binomial()).fit(cov_type='HAC', cov_kwds={'maxlags': L, 'use_correction': True})
out['binom'] = {'y': yb.tolist(), 'n': trials.tolist(), 'beta': gb.params.tolist(), 'se_naive': gb.bse.tolist(), 'se_hac': gbh.bse.tolist()}
# individuales: intercepto aleatorio y GEE
P, T = 80, 36; ids = np.repeat(np.arange(P), T); tt = np.tile(np.arange(T), P); Xi = design(tt, [14, 26], season=False); bi = rng.normal(0, 1.2, P)
yi = Xi @ np.array([8, 0.05, -1.5, -0.06, -1.0, 0.03]) + bi[ids] + rng.normal(0, 1, len(ids)); keep = rng.random(len(ids)) < 0.8
Xk, yk, ik, tk = Xi[keep], yi[keep], ids[keep], tt[keep]
ml = sm.MixedLM(yk, Xk, groups=ik).fit(reml=True); out['lmm'] = {'t': tk.tolist(), 'y': yk.tolist(), 'id': ik.tolist(), 'breaks': [14, 26], 'beta': ml.fe_params.tolist(), 'se': ml.bse_fe.tolist(), 'tau2': float(np.asarray(ml.cov_re)[0, 0]), 'sigma2': float(ml.scale)}
cnt = rng.poisson(np.exp(Xi @ np.array([0.8, 0.003, -0.4, -0.005, -0.2, 0.002]) + bi[ids] * 0.3)); ck = cnt[keep].astype(float)
ge = sm.GEE(ck, Xk, groups=ik, family=sm.families.Poisson(), cov_struct=sm.cov_struct.Exchangeable()).fit(); out['gee'] = {'y': ck.tolist(), 'beta': ge.params.tolist(), 'se': ge.bse.tolist()}
json.dump(out, open('tests/reference_its.json', 'w')); print('ok', {k: len(v.get('y', [])) for k, v in out.items()})

# ---- pronóstico
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.holtwinters import ExponentialSmoothing
n = 72; tt = np.arange(n); seas = 2.0 * np.sin(2 * np.pi * tt / 12) + 1.0 * np.cos(2 * np.pi * tt / 12 * 2)
ar = np.zeros(n)
for i in range(1, n): ar[i] = 0.6 * ar[i - 1] + rng.normal(0, 1)
yf = 10 + 0.15 * tt + seas + ar
Xf = np.column_stack([np.ones(n), tt, np.sin(2*np.pi*tt/12), np.cos(2*np.pi*tt/12), np.sin(4*np.pi*tt/12), np.cos(4*np.pi*tt/12)])
mo = sm.OLS(yf, Xf).fit(); th = np.arange(n, n + 12); Xh = np.column_stack([np.ones(12), th, np.sin(2*np.pi*th/12), np.cos(2*np.pi*th/12), np.sin(4*np.pi*th/12), np.cos(4*np.pi*th/12)])
pr = mo.get_prediction(Xh).summary_frame(alpha=0.05)
out['fc_linear'] = {'y': yf.tolist(), 'fc': pr['mean'].tolist(), 'lo': pr['obs_ci_lower'].tolist(), 'hi': pr['obs_ci_upper'].tolist()}
ya = np.zeros(120)
for i in range(1, 120): ya[i] = 0.6 * ya[i - 1] + rng.normal(0, 1)
ya = ya + 5
am = ARIMA(ya, order=(1, 0, 0), trend='c').fit(); af = am.get_forecast(10); out['fc_arima'] = {'y': ya.tolist(), 'phi': float(am.params[1]), 'const': float(am.params[0]), 'sigma2': float(am.params[2]), 'fc': af.predicted_mean.tolist(), 'se': af.se_mean.tolist()}
es = ExponentialSmoothing(yf, trend='add', damped_trend=True, seasonal='add', seasonal_periods=12, initialization_method='heuristic').fit(); out['fc_ets'] = {'y': yf.tolist(), 'sse': float(es.sse), 'fc': es.forecast(12).tolist()}
json.dump(out, open('tests/reference_its.json', 'w')); print('ok forecast')
