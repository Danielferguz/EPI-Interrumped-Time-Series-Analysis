// js/timeutil.js: lectura de fechas, unidad, índice de periodo y marcas del eje
import { parseTime, detectUnit, timeAxis, makeTicks } from '../js/timeutil.js';
let fails = 0, checks = 0; const ok = (n, c, x = '') => { checks++; if (!c) { fails++; console.log('FAIL', n, x); } else console.log('ok  ', n); };
const ymd = p => p && `${p.y}-${p.m}`;
ok('ISO, d/m/a, a/m, trimestres y año', ymd(parseTime('2020-04-15')) === '2020-4' && ymd(parseTime('15/04/2020')) === '2020-4' && ymd(parseTime('2020-04')) === '2020-4' && ymd(parseTime('Q3 2019')) === '2019-7' && ymd(parseTime('2019-Q2')) === '2019-4' && ymd(parseTime(2018)) === '2018-1' && ymd(parseTime(43831)) === '2020-1');
const months = Array.from({ length: 30 }, (_, i) => `${2018 + Math.floor(i / 12)}-${String(i % 12 + 1).padStart(2, '0')}-01`);
ok('detecta meses', detectUnit(months.map(parseTime)) === 'month' && timeAxis(months).unit === 'month');
const ax = timeAxis(months); ok('índice de periodo y rótulos', ax.N === 30 && ax.idx[0] === 0 && ax.idx[29] === 29 && ax.label(0) === '01-2018' && ax.label(14) === '03-2019');
const q = ['2014-01-01', '2014-04-01', '2014-07-01', '2014-10-01', '2015-01-01', '2015-04-01']; const aq = timeAxis(q); ok('trimestres', aq.unit === 'quarter' && aq.N === 6 && aq.label(4) === 'Q1 2015');
const gap = ['2020-01-01', '2020-02-01', '2020-05-01', '2020-06-01']; ok('huecos: el índice respeta los meses faltantes', timeAxis(gap, 'month').idx.join() === '0,1,4,5' && timeAxis(gap, 'month').N === 6);
const ix = timeAxis(['-24', '-23', '0', '36']); ok('números: el índice parte de 0 y conserva el rótulo original', ix.unit === 'index' && ix.idx.join() === '0,1,24,60' && ix.label(24) === '0');
const t = makeTicks(timeAxis(months)); ok('marcas en enero', t.ticks.length >= 2 && t.ticks[0].text === '01-2018' && t.ticks[0].i === 0 && t.ticks[1].i === 12);
const t2 = makeTicks(timeAxis(Array.from({ length: 36 }, (_, i) => `${2014 + Math.floor(i / 4)}-${String((i % 4) * 3 + 1).padStart(2, '0')}-01`))); ok('marcas trimestrales en T1', t2.ticks.every(k => k.i % 4 === 0));
console.log(`\n${checks - fails}/${checks} comprobaciones correctas`); process.exit(fails ? 1 : 0);
