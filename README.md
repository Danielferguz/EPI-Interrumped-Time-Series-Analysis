# EPI-ITS

Módulo de **EPI-LAB** para **analizar series temporales interrumpidas (ITS)**: cuánto cambia un resultado medido repetidamente cuando ocurre una intervención o un evento, con una a cuatro interrupciones consecutivas (p. ej., iSGLT2 → finerenona → AR-GLP-1), datos agregados o poblacionales (con paneles por región o centro) y datos individuales, y con pronóstico del contrafactual por varios métodos.

Todo se calcula en el navegador: los datos no salen de tu computadora. Sitio estático (GitHub Pages), sin instalación.

## Funciones

| Paso | Contenido |
|---|---|
| 1 · Guía y ejemplos | Guía de decisión y seis ejemplos **simulados**: proteinuria (UACR) antes y después de un iSGLT2 (1 interrupción), + finerenona (2) y + AR-GLP-1 (3); pruebas de VIH por región con un confinamiento por COVID-19 (panel binomial); UACR repetida por paciente (modelo mixto); encuesta individual binaria (GEE) |
| 2 · Datos | CSV/Excel; fechas en ISO, d/m/a, trimestres («Q2 2020»), años o números; unidad detectada (día, semana, mes, trimestre, año); agregados (resultado, denominador, grupo) o individuales (resultado, identificador) |
| 3 · Interrupciones y modelo | Hasta 4 interrupciones con nombre, efecto de nivel, pendiente o ambos, periodo de transición; continuo (MCO), Poisson con desplazamiento log(población), binomial agrupada, modelo mixto o GEE; errores de Newey–West o AR(1) de Prais–Winsten; estacionalidad de Fourier |
| 4 · Regresión segmentada | Figura con ajuste por tramo, IC, contrafactual y puente con el efecto; efecto inmediato, promedio y al final de cada tramo (diferencia absoluta y relativa por el método delta); coeficientes; ACF/PACF, Durbin–Watson y Ljung–Box |
| 5 · Pronóstico | Contrafactual con tendencia lineal + Fourier, Holt-Winters aditivo (tendencia amortiguada, estados iniciales por mínimos cuadrados) y ARIMA automático (AICc, ADF); validación de origen móvil dentro del periodo previo; diferencia observado − pronóstico |
| 6 · Grupos y paneles | Paneles pequeños por grupo, forest del efecto por grupo y efecto agrupado (DerSimonian–Laird) con I² |
| 7 · Reporte | Métodos (ES/EN), **código Python** (statsmodels, sin figuras) y tablas en CSV; figuras en SVG/PNG/TIFF a 300–1200 dpi |

## Validación de los motores

```bash
node tests/validate_its.mjs        # regresión segmentada, HAC, Prais–Winsten, Poisson/binomial, LMM y GEE frente a statsmodels
node tests/validate_forecast.mjs   # tendencia, Holt-Winters y ARIMA frente a statsmodels; origen móvil y efectos
node tests/validate_time.mjs       # fechas, unidades y ejes
python3 tools/make_reference.py    # regenera tests/reference_its.json (requiere numpy, pandas, statsmodels)
python3 tools/make_examples.py     # regenera los ejemplos simulados
```

El script de Python que descarga la herramienta reproduce exactamente los efectos mostrados en pantalla (comprobado en los seis ejemplos).

## Referencias

Bernal JL, Cummins S, Gasparrini A. *Int J Epidemiol* 2017;46:348–355 · Kontopantelis E et al. *BMJ* 2015;350:h2750 · Wagner AK et al. *J Clin Pharm Ther* 2002;27:299–309 · Schaffer AL et al. *BMC Med Res Methodol* 2021;21:58 · Turner SL et al. *BMC Med Res Methodol* 2021;21:134 · Newey WK, West KD. *Econometrica* 1987;55:703–708 · Hyndman RJ, Athanasopoulos G. *Forecasting: Principles and Practice*, 3.ª ed.

Los datos de los ejemplos son simulados. © Daniel Fernández-Guzmán · CC BY-NC-ND 4.0.
