# EBEC-scan: informe de validacion v2 (modo full)

## C1. Error tipo I bajo H0

| test | tasa agregada | min celda | max celda | celdas | liberales | conservadoras |
|---|---|---|---|---|---|---|
| scan_all | 0.046 | 0.037 | 0.053 | 12 | 0 | 1 |
| scan_tot | 0.046 | 0.041 | 0.051 | 12 | 0 | 0 |
| spiegelhalter | 0.054 | 0.045 | 0.061 | 12 | 2 | 0 |
| hl | 0.057 | 0.038 | 0.096 | 12 | 3 | 1 |
| cox | 0.052 | 0.042 | 0.060 | 12 | 1 | 0 |

**C1 (tamano) = PASA**  (tasa 0.0457; limite superior 0.0542; celdas liberales 0 de 12)

## Potencia por escenario (media de las celdas ejecutadas)

| escenario | grupo | scan_all | scan_tot | Spiegelhalter | HL | Cox |
|---|---|---|---|---|---|---|
| FN under-prediction c=0.05 d=1.5 | fn_under | 0.221 | 0.238 | 0.233 | 0.173 | 0.164 |
| FN under-prediction c=0.05 d=3.0 | fn_under | 0.613 | 0.631 | 0.555 | 0.603 | 0.525 |
| FN under-prediction c=0.10 d=1.5 | fn_under | 0.518 | 0.520 | 0.505 | 0.481 | 0.450 |
| FN under-prediction c=0.10 d=3.0 | fn_under | 0.915 | 0.918 | 0.854 | 0.951 | 0.901 |
| FP over-prediction c=0.05 d=1.5 | fp_over | 0.154 | 0.093 | 0.123 | 0.123 | 0.150 |
| FP over-prediction c=0.05 d=3.0 | fp_over | 0.345 | 0.180 | 0.213 | 0.251 | 0.324 |
| FP over-prediction c=0.10 d=1.5 | fp_over | 0.460 | 0.215 | 0.252 | 0.376 | 0.501 |
| FP over-prediction c=0.10 d=3.0 | fp_over | 0.689 | 0.315 | 0.332 | 0.635 | 0.747 |
| Global overconfidence k=0.50 | global | 0.956 | 0.969 | 0.979 | 0.890 | 0.947 |
| Global overconfidence k=0.75 | global | 0.690 | 0.730 | 0.780 | 0.552 | 0.640 |
| Inversion c=0.05 s=0.5 | inv | 0.304 | 0.346 | 0.320 | 0.292 | 0.266 |
| Inversion c=0.05 s=1.0 | inv | 0.496 | 0.527 | 0.537 | 0.425 | 0.450 |
| Inversion c=0.10 s=0.5 | inv | 0.487 | 0.528 | 0.537 | 0.423 | 0.441 |
| Inversion c=0.10 s=1.0 | inv | 0.751 | 0.776 | 0.810 | 0.612 | 0.697 |
| Rare confident FN c=0.01 m=0.02 | rare | 0.411 | 0.438 | 0.136 | 0.094 | 0.138 |
| Rare confident FN c=0.01 m=0.10 | rare | 0.154 | 0.165 | 0.079 | 0.072 | 0.064 |
| Rare confident FN c=0.02 m=0.02 | rare | 0.559 | 0.586 | 0.243 | 0.196 | 0.308 |
| Rare confident FN c=0.02 m=0.10 | rare | 0.279 | 0.288 | 0.156 | 0.110 | 0.113 |
| Rare confident FP c=0.01 m=0.02 | rare | 0.909 | 0.728 | 0.229 | 0.097 | 0.556 |
| Rare confident FP c=0.01 m=0.10 | rare | 0.746 | 0.261 | 0.173 | 0.080 | 0.303 |
| Rare confident FP c=0.02 m=0.02 | rare | 0.980 | 0.913 | 0.522 | 0.293 | 0.825 |
| Rare confident FP c=0.02 m=0.10 | rare | 0.871 | 0.518 | 0.354 | 0.230 | 0.573 |
| Rare confident mixed c=0.01 m=0.02 | rare | 0.846 | 0.582 | 0.179 | 0.084 | 0.338 |
| Rare confident mixed c=0.01 m=0.10 | rare | 0.563 | 0.207 | 0.116 | 0.074 | 0.148 |
| Rare confident mixed c=0.02 m=0.02 | rare | 0.955 | 0.793 | 0.378 | 0.149 | 0.694 |
| Rare confident mixed c=0.02 m=0.10 | rare | 0.793 | 0.380 | 0.226 | 0.103 | 0.340 |

**C2 (errores raros confiados) = PASA**  (10 de 12 escenarios con ventaja >= 0.15 sobre Spiegelhalter; se exigen 8)

**C3 (componentes por clase) = FALLA**  (FN: 0.000 de escenarios con ventaja >= 0.05; FP: 1.000)

**C4 (informativo) global:** scan_all - Spiegelhalter = Global overconfidence k=0.50: -0.023; Global overconfidence k=0.75: -0.090

## C5. Predicciones de un modelo ajustado (tipo I)

| n | ajustadas en la misma muestra | validacion cruzada 5-fold | repeticiones |
|---|---|---|---|
| 250 | 0.010 | 0.067 | 300 |
| 500 | 0.003 | 0.023 | 300 |

## Veredicto segun PREREGISTRATION.md
C1 y C2 PASAN, C3 FALLA: ventaja en errores raros; NO afirmar atribucion FN/FP.
