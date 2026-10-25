// Genera un script de Python (statsmodels) que reproduce el análisis: regresión segmentada, efectos frente al contrafactual y pronóstico. Sin figuras.
export function pythonScript(cfg) {
  const c = JSON.stringify(cfg, null, 2).replace(/\\/g, '\\\\').replace(/'''/g, "\\'\\'\\'");
  return String.raw`# Análisis de series temporales interrumpidas (ITS) · reproduce los resultados de la herramienta
# Requiere: python >= 3.10, numpy, pandas, statsmodels, scipy
import json, warnings
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import norm, t as tdist
from statsmodels.stats.stattools import durbin_watson
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.holtwinters import ExponentialSmoothing
from statsmodels.tsa.arima.model import ARIMA
warnings.filterwarnings('ignore')

CFG = json.loads('''${c}''')

# ---------------------------------------------------------------- datos y eje de tiempo
df = pd.read_excel(CFG['file']) if CFG['file'].lower().endswith(('.xlsx', '.xls')) else pd.read_csv(CFG['file'])
M = CFG['map']; unit = CFG['unit']
def period_index(col, unit):
    if unit == 'index':
        v = pd.to_numeric(col); return (v - v.min()).round().astype(int).to_numpy()
    d = pd.to_datetime(col, errors='coerce')
    if unit == 'month': o = d.dt.year * 12 + d.dt.month - 1
    elif unit == 'quarter': o = d.dt.year * 4 + (d.dt.month - 1) // 3
    elif unit == 'year': o = d.dt.year
    elif unit == 'week': o = (d - pd.Timestamp('1970-01-01')).dt.days // 7
    else: o = (d - pd.Timestamp('1970-01-01')).dt.days
    return (o - o.min()).astype(int).to_numpy()
df['_t'] = period_index(df[M['time']], unit); N = int(df['_t'].max()) + 1
G = CFG.get('group')
if G is not None and G != '__all': df = df[df[M['group']].astype(str) == str(G)]
fam, scale = CFG['family'], float(CFG['scale'])

# ---------------------------------------------------------------- diseño
B = CFG['breaks']; EFF = CFG['effects']; LAG = CFG['lag']; P, K = CFG['season']['P'], CFG['season']['K']
def row(t, cf=None, no_season=False):
    r = [1.0, float(t)]
    for j, T in enumerate(B):
        on = t >= T and (cf is None or j < cf)
        if EFF[j] != 'slope': r.append(1.0 if on else 0.0)
        if EFF[j] != 'level': r.append(float(t - T) if on else 0.0)
    if P > 1 and K > 0 and not no_season:
        for k in range(1, K + 1):
            a = 2 * np.pi * k * t / P
            r += [np.cos(np.pi * t)] if 2 * k == P else [np.sin(a), np.cos(a)]
    return np.array(r)
def in_lag(t): return any(T <= t < T + LAG for T in B)
LINK = {'gaussian': (lambda e: e, lambda e: 1.0), 'poisson': (np.exp, np.exp), 'binomial': (lambda e: 1 / (1 + np.exp(-e)), lambda e: np.exp(-e) / (1 + np.exp(-e)) ** 2)}
f, fd = LINK[fam]

def newey_lag(n): return max(1, int(np.floor(4 * (n / 100) ** (2 / 9))))
def prais_winsten(X, y, idx):
    rho = 0.0
    def tr(rho):
        cons = np.r_[False, np.diff(idx) == 1]
        Xs = np.where(cons[:, None], X - rho * np.vstack([X[:1], X[:-1]]), X * np.sqrt(1 - rho ** 2))
        ys = np.where(cons, y - rho * np.r_[y[:1], y[:-1]], y * np.sqrt(1 - rho ** 2)); return Xs, ys, cons
    for _ in range(200):
        Xs, ys, cons = tr(rho); b = np.linalg.lstsq(Xs, ys, rcond=None)[0]; e = y - X @ b
        nr = (e[1:] * e[:-1])[cons[1:]].sum() / (e[:-1] ** 2)[cons[1:]].sum()
        if abs(nr - rho) < 1e-10: rho = nr; break
        rho = nr
    Xs, ys, cons = tr(rho); b = np.linalg.lstsq(Xs, ys, rcond=None)[0]; res = ys - Xs @ b
    V = (res ** 2).sum() / (len(y) - X.shape[1]) * np.linalg.inv(Xs.T @ Xs); return b, V, rho

# ---------------------------------------------------------------- ajuste
if CFG['dtype'] == 'agg':
    g = df.groupby('_t').agg(y=(M['y'], 'sum' if fam != 'gaussian' else 'mean'), n=((M['n'] if M['n'] else M['y']), 'sum')).reindex(range(N))
    y_all, n_all = g['y'].to_numpy(float), g['n'].to_numpy(float)
    ok = [t for t in range(N) if np.isfinite(y_all[t]) and not in_lag(t) and (fam == 'gaussian' or n_all[t] > 0)]
    X = np.array([row(t) for t in ok]); y = y_all[ok]; n = n_all[ok]; idx = np.array(ok)
    L = CFG['nwL'] or newey_lag(len(ok)); rho = None
    kw = {'cov_type': 'HAC', 'cov_kwds': {'maxlags': L, 'use_correction': True}} if CFG['ac'] == 'nw' else {}
    if fam == 'gaussian':
        if CFG['ac'] == 'ar1':
            beta, V, rho = prais_winsten(X, y, idx); res_ = sm.OLS(y, X).fit()
        else:
            res_ = sm.OLS(y, X).fit(**kw); beta, V = res_.params, res_.cov_params()
        resid = sm.OLS(y, X).fit().resid; q = tdist.ppf(1 - (1 - CFG['level']) / 2, len(y) - X.shape[1])
    else:
        mod = sm.GLM(y, X, family=sm.families.Poisson(), offset=np.log(n)) if fam == 'poisson' else sm.GLM(np.column_stack([y, n - y]), X, family=sm.families.Binomial())
        res_ = mod.fit(**kw) if CFG['ac'] == 'nw' else (mod.fit(scale='X2') if mod.fit().pearson_chi2 / mod.fit().df_resid > 1 else mod.fit())
        beta, V = res_.params, res_.cov_params(); resid = mod.fit().resid_pearson; q = norm.ppf(1 - (1 - CFG['level']) / 2)
    obs = y_all if fam == 'gaussian' else y_all / n_all * scale
else:
    d = df[[M['y'], M['id'], '_t']].dropna(); d = d[[not in_lag(t) for t in d['_t']]]
    X = np.array([row(t) for t in d['_t']]); y = d[M['y']].to_numpy(float); ids = d[M['id']].to_numpy()
    if fam == 'gaussian': res_ = sm.MixedLM(y, X, groups=ids).fit(reml=True); beta, V = res_.fe_params, np.asarray(res_.cov_params())[:X.shape[1], :X.shape[1]]
    else: res_ = sm.GEE(y, X, groups=ids, family=sm.families.Poisson() if fam == 'poisson' else sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit(); beta, V = res_.params, res_.cov_params()
    q = norm.ppf(1 - (1 - CFG['level']) / 2); resid = None; rho = None; N = N
    a = df.groupby('_t')[M['y']].mean().reindex(range(N)); obs = a.to_numpy(float) * (1 if fam == 'gaussian' else scale)
beta = np.asarray(beta); V = np.asarray(V); se = np.sqrt(np.diag(V))
names = ['Intercepto', 'Tendencia previa'] + [f'{"Nivel" if k == "step" else "Pendiente"} · {CFG["labels"][j]}' for j in range(len(B)) for k in (['step'] if EFF[j] == 'level' else ['slope'] if EFF[j] == 'slope' else ['step', 'slope'])] + ['Estacional'] * (len(beta) - 2 - sum(1 if e != 'both' else 2 for e in EFF))
coef = pd.DataFrame({'term': names[:len(beta)], 'estimate': beta, 'se': se, 'lo': beta - q * se, 'hi': beta + q * se})
if fam != 'gaussian': coef['exp_est'] = np.exp(coef['estimate']); coef['exp_lo'] = np.exp(coef['lo']); coef['exp_hi'] = np.exp(coef['hi'])
print('== Coeficientes =='); print(coef.round(4).to_string(index=False))
if resid is not None and len(resid) > 12:
    lb = acorr_ljungbox(resid, lags=[min(10, max(1, len(resid) // 5))], return_df=True); print(f'\nDurbin-Watson = {durbin_watson(resid):.2f}; Ljung-Box P = {lb["lb_pvalue"].iloc[0]:.3f}' + (f'; rho (Prais-Winsten) = {rho:.3f}' if rho is not None else ''))

# ---------------------------------------------------------------- efectos frente al contrafactual (método delta)
def contrast(j, ts):
    A = B_ = 0.0; gA = np.zeros(len(beta)); gB = np.zeros(len(beta))
    for t in ts:
        u, v = row(t), row(t, cf=j); eu, ev = u @ beta, v @ beta; A += f(eu); B_ += f(ev); gA += fd(eu) * u; gB += fd(ev) * v
    k = len(ts); A, B_, gA, gB = A / k, B_ / k, gA / k, gB / k; sd = np.sqrt((gA - gB) @ V @ (gA - gB)); gr = (gA * B_ - A * gB) / B_ ** 2; sr = np.sqrt(gr @ V @ gr)
    return dict(obs=scale * A, cf=scale * B_, diff=scale * (A - B_), diff_lo=scale * (A - B_ - q * sd), diff_hi=scale * (A - B_ + q * sd), rel=A / B_ - 1, rel_lo=A / B_ - 1 - q * sr, rel_hi=A / B_ - 1 + q * sr)
rows = []
for j, T in enumerate(B):
    nxt = B[j + 1] if j + 1 < len(B) else N; t0 = min(T + LAG, nxt - 1)
    for w, ts in (('inmediato', [t0]), ('promedio', list(range(t0, nxt))), ('final del tramo', [nxt - 1])): rows.append({'interrupcion': CFG['labels'][j], 'ventana': w, **contrast(j, ts)})
print('\n== Efecto frente a continuar la tendencia previa =='); print(pd.DataFrame(rows).round(3).to_string(index=False))

# ---------------------------------------------------------------- pronóstico del contrafactual (tendencia lineal, Holt-Winters, ARIMA)
FC = CFG.get('forecast')
if FC:
    def fourier(t, P, K): return np.column_stack([fn(2 * np.pi * k * np.arange(t[0], t[1]) / P) for k in range(1, K + 1) for fn in (np.sin, np.cos)]) if P > 1 and K > 0 else np.empty((t[1] - t[0], 0))
    ob = pd.Series(obs).interpolate(limit_direction='both').to_numpy(); out = []
    for j, T in enumerate(B):
        nxt = B[j + 1] if j + 1 < len(B) else N; a0 = (B[j - 1] + LAG if j and FC['train'] == 'previous' and T - (B[j - 1] + LAG) >= 8 else 0); tr = ob[a0:T]; h = nxt - T; ts = range(T + LAG, nxt)
        for m in FC['methods']:
            if m == 'linear':
                Xt = np.column_stack([np.ones(len(tr)), np.arange(a0, T), fourier((a0, T), P, K)]); fit = sm.OLS(tr, Xt).fit(); Xf = np.column_stack([np.ones(h), np.arange(T, nxt), fourier((T, nxt), P, K)]); pr = fit.get_prediction(Xf).summary_frame(alpha=1 - CFG['level']); fc, lo, hi = pr['mean'].to_numpy(), pr['obs_ci_lower'].to_numpy(), pr['obs_ci_upper'].to_numpy()
            elif m == 'ets':
                best = None
                for trd, dmp in ((None, False), ('add', False), ('add', True)):
                    try: e = ExponentialSmoothing(tr, trend=trd, damped_trend=dmp, seasonal='add' if P > 1 and len(tr) >= 2 * P else None, seasonal_periods=P if P > 1 and len(tr) >= 2 * P else None, initialization_method='estimated').fit()
                    except Exception: continue
                    if best is None or e.aicc < best.aicc: best = e
                fc = best.forecast(h); sg = np.sqrt(best.sse / len(tr)); lo, hi = fc - norm.ppf(0.975) * sg * np.sqrt(np.arange(1, h + 1)), fc + norm.ppf(0.975) * sg * np.sqrt(np.arange(1, h + 1))
            elif m == 'arima':
                d_ = 0 if adfuller(tr, regression='ct')[0] < -3.41 else 1; best = None; ex = fourier((a0, T), P, K); exf = fourier((T, nxt), P, K)
                for p in range(3):
                    for qq in range(3):
                        try: e = ARIMA(tr, exog=ex if ex.shape[1] else None, order=(p, d_, qq), trend='c' if d_ == 0 else 't').fit()
                        except Exception: continue
                        if best is None or e.aicc < best.aicc: best = e
                g_ = best.get_forecast(h, exog=exf if exf.shape[1] else None); fc = g_.predicted_mean; sgf = g_.se_mean; lo, hi = fc - norm.ppf(0.975) * sgf, fc + norm.ppf(0.975) * sgf
            else: continue
            idx_ = [t - T for t in ts if np.isfinite(obs[t])]; A_ = np.mean([obs[T + i] for i in idx_]); F_ = np.mean(np.asarray(fc)[idx_]); sm_ = np.mean((np.asarray(hi) - np.asarray(lo))[idx_] / (2 * norm.ppf(0.975)))
            out.append({'interrupcion': CFG['labels'][j], 'metodo': m, 'observado': A_, 'pronostico': F_, 'diferencia': A_ - F_, 'ic_inf': A_ - F_ - norm.ppf(0.975) * sm_, 'ic_sup': A_ - F_ + norm.ppf(0.975) * sm_, 'cambio_relativo': A_ / F_ - 1})
    print('\n== Efecto frente al pronóstico (IC conservador) =='); print(pd.DataFrame(out).round(3).to_string(index=False))
`;
}
