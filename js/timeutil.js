// Fechas de las series: lectura de formatos habituales, detección de la unidad (día, semana, mes, trimestre, año o índice),
// índice de periodo, rótulos y marcas del eje. Sin dependencias del navegador.
const MS = 86400000, DAYS = { day: 1, week: 7, month: 30.44, quarter: 91.31, year: 365.25 };
export const UNITS = [['auto', 'Detectar'], ['month', 'Mes'], ['quarter', 'Trimestre'], ['year', 'Año'], ['week', 'Semana'], ['day', 'Día'], ['index', 'Número de periodo']];
export const PERIOD = { month: 12, quarter: 4, week: 52, day: 7, year: 0, index: 0 };
const pad = n => String(n).padStart(2, '0');

// → {day: días desde 1970-01-01 (UTC), y, m (1–12)} o null; admite ISO, d/m/a, a/m, año, «Qn aaaa», «aaaa-Qn», Date y serie de Excel
export function parseTime(v) {
  if (v instanceof Date && !isNaN(v)) return fromYMD(v.getUTCFullYear(), v.getUTCMonth() + 1, v.getUTCDate());
  if (typeof v === 'number') { if (v > 20000 && v < 80000) return parseTime(new Date(Math.round((v - 25569) * MS))); if (v >= 1800 && v <= 2200 && Number.isInteger(v)) return fromYMD(v, 1, 1); return null; }
  const s = String(v ?? '').trim(); if (!s) return null; let m;
  if ((m = /^(\d{4})[-/.](\d{1,2})(?:[-/.](\d{1,2}))?/.exec(s))) return fromYMD(+m[1], +m[2], m[3] ? +m[3] : 1);
  if ((m = /^(\d{1,2})[-/.](\d{1,2})[-/.](\d{4})$/.exec(s))) return fromYMD(+m[3], +m[2], +m[1]);
  if ((m = /^(\d{1,2})[-/.](\d{4})$/.exec(s))) return fromYMD(+m[2], +m[1], 1);
  if ((m = /^(?:Q([1-4])[\s-]*(\d{4})|(\d{4})[\s-]*Q([1-4]))$/i.exec(s))) return fromYMD(+(m[2] || m[3]), 3 * (+(m[1] || m[4]) - 1) + 1, 1);
  if ((m = /^(\d{4})$/.exec(s))) return fromYMD(+m[1], 1, 1);
  const d = new Date(s); return isNaN(d) ? null : fromYMD(d.getFullYear(), d.getMonth() + 1, d.getDate());
}
const fromYMD = (y, m, d) => ({ y, m, day: Math.round(Date.UTC(y, m - 1, d) / MS) });

// unidad por la mediana de las diferencias entre fechas distintas
export function detectUnit(parsed) {
  const days = [...new Set(parsed.filter(Boolean).map(p => p.day))].sort((a, b) => a - b); if (days.length < 3) return 'index';
  const d = days.slice(1).map((v, i) => v - days[i]).sort((a, b) => a - b), med = d[Math.floor(d.length / 2)];
  return med <= 1.5 ? 'day' : med >= 6 && med <= 8 ? 'week' : med >= 27 && med <= 32 ? 'month' : med >= 88 && med <= 93 ? 'quarter' : med >= 360 && med <= 370 ? 'year' : 'index';
}
const ord = (p, unit) => (unit === 'month' ? p.y * 12 + p.m - 1 : unit === 'quarter' ? p.y * 4 + Math.floor((p.m - 1) / 3) : unit === 'year' ? p.y : unit === 'week' ? Math.floor(p.day / 7) : p.day);
// índice entero del periodo desde el primer dato; periods = [{index, label, sort}] para todo el rango (los huecos quedan vacíos)
export function timeAxis(rawValues, unitIn = 'auto') {
  const parsed = rawValues.map(parseTime), ok = parsed.filter(Boolean).length; let unit = unitIn;
  if (unit === 'auto') unit = ok >= rawValues.length * 0.9 ? detectUnit(parsed) : 'index';
  if (unit === 'index' || ok < rawValues.length * 0.9) { const nums = rawValues.map(Number); if (nums.every(Number.isFinite)) { const lo = Math.min(...nums); return { unit: 'index', idx: nums.map(v => Math.round(v - lo)), N: Math.max(...nums) - lo + 1, label: i => String(Math.round(i + lo)), start: lo, step: 1 }; } throw new Error('No se pudo interpretar la columna de tiempo como fechas ni como números.'); }
  if (parsed.some(p => !p)) throw new Error('Hay valores de tiempo que no se pueden leer como fecha.');
  const o = parsed.map(p => ord(p, unit)), lo = Math.min(...o), idx = o.map(v => v - lo), N = Math.max(...idx) + 1;
  const label = i => { const k = lo + i; if (unit === 'month') return `${pad((k % 12) + 1)}-${Math.floor(k / 12)}`; if (unit === 'quarter') return `Q${(k % 4) + 1} ${Math.floor(k / 4)}`; if (unit === 'year') return String(k); const dt = new Date(unit === 'week' ? k * 7 * MS : k * MS); return dt.toISOString().slice(0, 10); };
  return { unit, idx, N, label, start: lo, step: 1 };
}
// marcas del eje: ≈ 8–10 rótulos, preferentemente en enero; minor = todos los periodos si son pocos
export function makeTicks(ax, maxTicks = 10) {
  const { unit, N, label, start } = ax, per = unit === 'month' ? 12 : unit === 'quarter' ? 4 : 1, ticks = [];
  if (unit === 'month' || unit === 'quarter' || unit === 'year') {
    const years = Math.ceil(N / per), step = Math.max(1, Math.ceil(years / maxTicks)), first = unit === 'year' ? 0 : (per - (((start % per) + per) % per)) % per;
    for (let i = unit === 'year' ? 0 : first; i < N; i += unit === 'year' ? step : step * per) ticks.push({ i, text: unit === 'year' ? label(i) : unit === 'month' ? `01-${label(i).slice(3)}` : `Q1 ${label(i).slice(3)}` });
    if (!ticks.length) ticks.push({ i: 0, text: label(0) });
  } else { const step = Math.max(1, Math.ceil(N / maxTicks)); for (let i = 0; i < N; i += step) ticks.push({ i, text: label(i) }); }
  return { ticks, minor: N <= 80 ? Array.from({ length: N }, (_, i) => i) : [] };
}
