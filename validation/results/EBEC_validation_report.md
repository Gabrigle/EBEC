# EBEC / EBEC_FN / EBEC_FP - informe de validacion por simulacion

Generado: 2026-10-18 19:36 | analisis de EBEC_all_pvalues.csv
Replicas H0 por celda y test: 2000; alpha = 0.05

## 1. Error tipo I (H0 verdadera: y ~ Bernoulli(p), p fijo)

Criterio fijado de antemano: la celda pasa si el IC de Wilson al 95% de la tasa de rechazo contiene 0.05.

| test | q | tasa media | min celda | max celda | celdas OK |
|---|---|---|---|---|---|
| EBEC_FN | 0.05 | 0.049 | 0.044 | 0.057 | 12/12 |
| EBEC_FN | 0.1 | 0.049 | 0.043 | 0.054 | 12/12 |
| EBEC_FN | 0.2 | 0.049 | 0.042 | 0.056 | 12/12 |
| EBEC_FP | 0.05 | 0.048 | 0.038 | 0.056 | 10/12 |
| EBEC_FP | 0.1 | 0.050 | 0.040 | 0.061 | 11/12 |
| EBEC_FP | 0.2 | 0.050 | 0.040 | 0.059 | 11/12 |
| EBEC | 0.05 | 0.049 | 0.040 | 0.054 | 12/12 |
| EBEC | 0.1 | 0.049 | 0.043 | 0.056 | 12/12 |
| EBEC | 0.2 | 0.050 | 0.043 | 0.060 | 11/12 |
| Spiegelhalter | - | 0.054 | 0.046 | 0.061 | 11/12 |

Celdas que NO pasan: 6 de 120.

| test | n | prev | q | tasa | veredicto |
|---|---|---|---|---|---|
| EBEC_FP | 1000 | 0.3 | 0.05 | 0.040 | CONSERVADOR |
| EBEC_FP | 1000 | 0.5 | 0.10 | 0.061 | ANTICONSERVADOR |
| EBEC_FP |  500 | 0.1 | 0.05 | 0.038 | CONSERVADOR |
| EBEC_FP |  500 | 0.1 | 0.20 | 0.040 | CONSERVADOR |
| EBEC |  100 | 0.5 | 0.20 | 0.060 | ANTICONSERVADOR |
| Spiegelhalter |  250 | 0.1 | NA | 0.061 | ANTICONSERVADOR |

Nota: con n pequeno y q pequeno el subconjunto de cola tiene pocas observaciones y el p-valor es muy discreto (tests conservadores).

## 2. Potencia (media sobre n y prevalencia; por escenario y q)

Detalle completo en EBEC_power_by_cell.csv. Se informan TODOS los escenarios.

**FN_infraprediccion c=0.05 d=1.5**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.027 |
| EBEC_FN | 0.1 | 0.040 |
| EBEC_FN | 0.2 | 0.058 |
| EBEC_FP | 0.05 | 0.034 |
| EBEC_FP | 0.1 | 0.029 |
| EBEC_FP | 0.2 | 0.028 |
| EBEC | 0.05 | 0.026 |
| EBEC | 0.1 | 0.039 |
| EBEC | 0.2 | 0.046 |
| Spiegelhalter | - | 0.233 |

**FN_infraprediccion c=0.05 d=3.0**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.014 |
| EBEC_FN | 0.1 | 0.022 |
| EBEC_FN | 0.2 | 0.058 |
| EBEC_FP | 0.05 | 0.020 |
| EBEC_FP | 0.1 | 0.017 |
| EBEC_FP | 0.2 | 0.013 |
| EBEC | 0.05 | 0.009 |
| EBEC | 0.1 | 0.039 |
| EBEC | 0.2 | 0.074 |
| Spiegelhalter | - | 0.558 |

**FN_infraprediccion c=0.10 d=1.5**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.016 |
| EBEC_FN | 0.1 | 0.023 |
| EBEC_FN | 0.2 | 0.046 |
| EBEC_FP | 0.05 | 0.023 |
| EBEC_FP | 0.1 | 0.019 |
| EBEC_FP | 0.2 | 0.020 |
| EBEC | 0.05 | 0.017 |
| EBEC | 0.1 | 0.033 |
| EBEC | 0.2 | 0.047 |
| Spiegelhalter | - | 0.494 |

**FN_infraprediccion c=0.10 d=3.0**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.002 |
| EBEC_FN | 0.1 | 0.004 |
| EBEC_FN | 0.2 | 0.024 |
| EBEC_FP | 0.05 | 0.013 |
| EBEC_FP | 0.1 | 0.008 |
| EBEC_FP | 0.2 | 0.003 |
| EBEC | 0.05 | 0.001 |
| EBEC | 0.1 | 0.011 |
| EBEC | 0.2 | 0.093 |
| Spiegelhalter | - | 0.866 |

**FP_sobreprediccion c=0.05 d=1.5**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.040 |
| EBEC_FN | 0.1 | 0.040 |
| EBEC_FN | 0.2 | 0.036 |
| EBEC_FP | 0.05 | 0.023 |
| EBEC_FP | 0.1 | 0.025 |
| EBEC_FP | 0.2 | 0.040 |
| EBEC | 0.05 | 0.080 |
| EBEC | 0.1 | 0.049 |
| EBEC | 0.2 | 0.037 |
| Spiegelhalter | - | 0.124 |

**FP_sobreprediccion c=0.05 d=3.0**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.034 |
| EBEC_FN | 0.1 | 0.031 |
| EBEC_FN | 0.2 | 0.026 |
| EBEC_FP | 0.05 | 0.023 |
| EBEC_FP | 0.1 | 0.027 |
| EBEC_FP | 0.2 | 0.050 |
| EBEC | 0.05 | 0.096 |
| EBEC | 0.1 | 0.049 |
| EBEC | 0.2 | 0.039 |
| Spiegelhalter | - | 0.217 |

**FP_sobreprediccion c=0.10 d=1.5**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.031 |
| EBEC_FN | 0.1 | 0.023 |
| EBEC_FN | 0.2 | 0.023 |
| EBEC_FP | 0.05 | 0.011 |
| EBEC_FP | 0.1 | 0.009 |
| EBEC_FP | 0.2 | 0.026 |
| EBEC | 0.05 | 0.133 |
| EBEC | 0.1 | 0.045 |
| EBEC | 0.2 | 0.042 |
| Spiegelhalter | - | 0.241 |

**FP_sobreprediccion c=0.10 d=3.0**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.027 |
| EBEC_FN | 0.1 | 0.026 |
| EBEC_FN | 0.2 | 0.023 |
| EBEC_FP | 0.05 | 0.004 |
| EBEC_FP | 0.1 | 0.009 |
| EBEC_FP | 0.2 | 0.025 |
| EBEC | 0.05 | 0.186 |
| EBEC | 0.1 | 0.035 |
| EBEC | 0.2 | 0.037 |
| Spiegelhalter | - | 0.330 |

**inversion c=0.05 s=0.5**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.026 |
| EBEC_FN | 0.1 | 0.024 |
| EBEC_FN | 0.2 | 0.029 |
| EBEC_FP | 0.05 | 0.043 |
| EBEC_FP | 0.1 | 0.046 |
| EBEC_FP | 0.2 | 0.041 |
| EBEC | 0.05 | 0.026 |
| EBEC | 0.1 | 0.031 |
| EBEC | 0.2 | 0.124 |
| Spiegelhalter | - | 0.320 |

**inversion c=0.05 s=1.0**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.021 |
| EBEC_FN | 0.1 | 0.023 |
| EBEC_FN | 0.2 | 0.023 |
| EBEC_FP | 0.05 | 0.046 |
| EBEC_FP | 0.1 | 0.050 |
| EBEC_FP | 0.2 | 0.056 |
| EBEC | 0.05 | 0.020 |
| EBEC | 0.1 | 0.024 |
| EBEC | 0.2 | 0.202 |
| Spiegelhalter | - | 0.557 |

**inversion c=0.10 s=0.5**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.023 |
| EBEC_FN | 0.1 | 0.027 |
| EBEC_FN | 0.2 | 0.023 |
| EBEC_FP | 0.05 | 0.049 |
| EBEC_FP | 0.1 | 0.057 |
| EBEC_FP | 0.2 | 0.063 |
| EBEC | 0.05 | 0.019 |
| EBEC | 0.1 | 0.024 |
| EBEC | 0.2 | 0.211 |
| Spiegelhalter | - | 0.534 |

**inversion c=0.10 s=1.0**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.010 |
| EBEC_FN | 0.1 | 0.015 |
| EBEC_FN | 0.2 | 0.020 |
| EBEC_FP | 0.05 | 0.045 |
| EBEC_FP | 0.1 | 0.053 |
| EBEC_FP | 0.2 | 0.062 |
| EBEC | 0.05 | 0.013 |
| EBEC | 0.1 | 0.013 |
| EBEC | 0.2 | 0.057 |
| Spiegelhalter | - | 0.811 |

**sobreconfianza global k=0.50**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.002 |
| EBEC_FN | 0.1 | 0.002 |
| EBEC_FN | 0.2 | 0.004 |
| EBEC_FP | 0.05 | 0.105 |
| EBEC_FP | 0.1 | 0.146 |
| EBEC_FP | 0.2 | 0.198 |
| EBEC | 0.05 | 0.003 |
| EBEC | 0.1 | 0.003 |
| EBEC | 0.2 | 0.006 |
| Spiegelhalter | - | 0.979 |

**sobreconfianza global k=0.75**

| test | q | potencia media |
|---|---|---|
| EBEC_FN | 0.05 | 0.007 |
| EBEC_FN | 0.1 | 0.007 |
| EBEC_FN | 0.2 | 0.007 |
| EBEC_FP | 0.05 | 0.079 |
| EBEC_FP | 0.1 | 0.097 |
| EBEC_FP | 0.2 | 0.121 |
| EBEC | 0.05 | 0.003 |
| EBEC | 0.1 | 0.008 |
| EBEC | 0.2 | 0.210 |
| Spiegelhalter | - | 0.780 |

## 3. Cobertura del IC bootstrap (la verdad es exceso = 0 bajo H0)

| test | n | q | cobertura de 0 | anchura media |
|---|---|---|---|---|
| no calculado | NA | NA |   NA |   NA |

Cobertura esperada: 0.95. Si es claramente menor, el IC no esta validado y no debe usarse en el paper.
