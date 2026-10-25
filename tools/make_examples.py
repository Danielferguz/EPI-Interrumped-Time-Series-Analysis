# Datos SIMULADOS para los ejemplos de EPI-ITS (ningún dato es real).
# Uso: python3 tools/make_examples.py   (requiere numpy y pandas)
import numpy as np, pandas as pd
def months(start, n): return pd.date_range(start, periods=n, freq='MS')
def ar1(rng, n, rho, sd):
    e = np.zeros(n)
    for i in range(1, n): e[i] = rho * e[i - 1] + rng.normal(0, sd)
    return e

# ---- 1-3. proteinuria (UACR) mensual en una cohorte de ERC: SGLT2i, luego finerenona, luego AR-GLP-1 (1, 2 y 3 interrupciones)
def uacr(rng, n, steps, name):
    t = np.arange(n); ln = np.log(450) + 0.004 * t + ar1(rng, n, 0.5, 0.05) + rng.normal(0, 0.03, n)
    for (T, lev, sl) in steps: ln += np.where(t >= T, lev + sl * (t - T), 0)
    npat = rng.integers(140, 260, n); d = pd.DataFrame({'month': months('2017-01-01', n).strftime('%Y-%m-%d'), 'uacr_gm': np.exp(ln).round(1), 'ln_uacr': ln.round(4), 'patients': npat})
    d.to_csv(f'examples/{name}.csv', index=False); return d
rng = np.random.default_rng(101)
uacr(rng, 72, [(30, -0.26, -0.005)], 'uacr_sglt2i')
uacr(rng, 84, [(24, -0.26, -0.005), (48, -0.17, -0.003)], 'uacr_sglt2i_finerenone')
uacr(rng, 96, [(24, -0.26, -0.005), (48, -0.17, -0.003), (72, -0.11, -0.002)], 'uacr_sglt2i_finerenone_glp1ra')

# ---- 4. pruebas de VIH en el último año por región (trimestral 2014-2022; confinamiento COVID-19 en 2020 T2); panel con denominadores
regions = ['Amazonas', 'Ancash', 'Arequipa', 'Callao', 'Cusco', 'Ica', 'Junin', 'La Libertad', 'Lima', 'Loreto', 'Piura', 'Puno']
rng = np.random.default_rng(202); N = 36; T = 25; t = np.arange(N); rows = []
for r in regions:
    base = rng.uniform(12, 24); sl = rng.uniform(0.15, 0.45); lev = -rng.uniform(4, 9); sl2 = rng.uniform(0, 0.25); e = ar1(rng, N, 0.4, 0.8)
    for i in range(N):
        p = np.clip(base + sl * i + e[i] + (lev + sl2 * (i - T) if i >= T else 0), 3, 60); n = int(rng.integers(900, 2400)); k = rng.binomial(n, p / 100)
        rows.append({'region': r, 'quarter': months('2014-01-01', N * 3)[3 * i].strftime('%Y-%m-%d'), 'tested': int(k), 'surveyed': n, 'pct_tested': round(100 * k / n, 2)})
pd.DataFrame(rows).to_csv('examples/hiv_testing_regions.csv', index=False)

# ---- 5. datos individuales: UACR repetida por paciente, tiempo en meses desde el inicio del SGLT2i (0); finerenona al mes 12 y AR-GLP-1 al mes 24
rng = np.random.default_rng(303); P = 300; rows = []
for pid in range(1, P + 1):
    b = rng.normal(0, 0.9); age = int(rng.normal(62, 11)); sex = rng.choice(['Female', 'Male'], p=[.42, .58]); base = np.log(430)
    visits = np.sort(rng.choice(np.arange(-24, 37), size=int(rng.integers(10, 18)), replace=False)); u = 0.0
    for m in visits:
        u = 0.4 * u + rng.normal(0, 0.18); ln = base + b + 0.004 * m + u
        if m >= 0: ln += -0.26 - 0.006 * m
        if m >= 12: ln += -0.17 - 0.004 * (m - 12)
        if m >= 24: ln += -0.11 - 0.003 * (m - 24)
        rows.append({'patient_id': pid, 'month': int(m), 'uacr': round(float(np.exp(ln)), 1), 'ln_uacr': round(float(ln), 4), 'age': age, 'sex': sex})
pd.DataFrame(rows).to_csv('examples/uacr_individual.csv', index=False)

# ---- 6. encuesta individual: persona-trimestre con prueba de VIH en el último año (0/1), 8 regiones; confinamiento en 2020 T2
rng = np.random.default_rng(404); rows = []; rid = 0
for r in regions[:8]:
    base = rng.uniform(-1.9, -1.1); sl = rng.uniform(0.01, 0.03); lev = -rng.uniform(0.3, 0.6)
    for i in range(N):
        p = 1 / (1 + np.exp(-(base + sl * i + (lev if i >= T else 0)))); y = rng.random(30) < p
        for v in y: rid += 1; rows.append({'respondent_id': rid, 'region': r, 'quarter': months('2014-01-01', N * 3)[3 * i].strftime('%Y-%m-%d'), 'tested': int(v)})
pd.DataFrame(rows).to_csv('examples/hiv_testing_survey.csv', index=False)
print('ok')
