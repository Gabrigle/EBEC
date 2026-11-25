# EBEC-scan: informe de validacion v3 (scan_fp) (modo full)

## C1. Error tipo I bajo H0

| test | tasa agregada | min celda | max celda | celdas | liberales | conservadoras |
|---|---|---|---|---|---|---|
| scan_fp | 0.049 | 0.044 | 0.053 | 12 | 0 | 0 |
| scan_both | 0.048 | 0.043 | 0.053 | 12 | 0 | 0 |
| scan_tot | 0.047 | 0.042 | 0.052 | 12 | 0 | 0 |
| spiegelhalter | 0.055 | 0.049 | 0.061 | 12 | 3 | 0 |
| hl | 0.058 | 0.042 | 0.088 | 12 | 4 | 0 |
| cox | 0.052 | 0.045 | 0.060 | 12 | 1 | 0 |

**C1 (tamano) = PASA**  (tasa 0.0490; limite superior 0.0542; celdas liberales 0 de 12)

## Potencia por escenario (media de las celdas ejecutadas)

| escenario | grupo | scan_fp | scan_both | scan_tot | Spiegelhalter | HL | Cox |
|---|---|---|---|---|---|---|---|
| FN under-prediction c=0.05 d=1.5 | fn_under | 0.189 | 0.214 | 0.230 | 0.232 | 0.175 | 0.169 |
| FN under-prediction c=0.05 d=3.0 | fn_under | 0.572 | 0.624 | 0.642 | 0.563 | 0.614 | 0.553 |
| FN under-prediction c=0.10 d=1.5 | fn_under | 0.444 | 0.513 | 0.506 | 0.496 | 0.490 | 0.457 |
| FN under-prediction c=0.10 d=3.0 | fn_under | 0.885 | 0.920 | 0.916 | 0.858 | 0.953 | 0.906 |
| FP over-prediction c=0.05 d=1.5 | fp_over | 0.186 | 0.174 | 0.101 | 0.135 | 0.117 | 0.173 |
| FP over-prediction c=0.05 d=3.0 | fp_over | 0.353 | 0.338 | 0.185 | 0.211 | 0.248 | 0.333 |
| FP over-prediction c=0.10 d=1.5 | fp_over | 0.481 | 0.466 | 0.219 | 0.251 | 0.365 | 0.493 |
| FP over-prediction c=0.10 d=3.0 | fp_over | 0.695 | 0.692 | 0.321 | 0.336 | 0.613 | 0.752 |
| Global overconfidence k=0.50 | global | 0.967 | 0.963 | 0.973 | 0.980 | 0.902 | 0.950 |
| Global overconfidence k=0.75 | global | 0.686 | 0.686 | 0.723 | 0.774 | 0.553 | 0.630 |
| Inversion c=0.05 s=0.5 | inv | 0.298 | 0.297 | 0.352 | 0.319 | 0.297 | 0.257 |
| Inversion c=0.05 s=1.0 | inv | 0.501 | 0.504 | 0.538 | 0.547 | 0.428 | 0.458 |
| Inversion c=0.10 s=0.5 | inv | 0.501 | 0.499 | 0.544 | 0.540 | 0.427 | 0.447 |
| Inversion c=0.10 s=1.0 | inv | 0.761 | 0.752 | 0.790 | 0.801 | 0.617 | 0.703 |
| Rare confident FN c=0.01 m=0.02 | rare | 0.406 | 0.414 | 0.430 | 0.117 | 0.082 | 0.125 |
| Rare confident FN c=0.01 m=0.10 | rare | 0.153 | 0.181 | 0.184 | 0.093 | 0.063 | 0.070 |
| Rare confident FN c=0.02 m=0.02 | rare | 0.556 | 0.557 | 0.583 | 0.229 | 0.184 | 0.297 |
| Rare confident FN c=0.02 m=0.10 | rare | 0.246 | 0.269 | 0.289 | 0.148 | 0.112 | 0.105 |
| Rare confident FP c=0.01 m=0.02 | rare | 0.929 | 0.926 | 0.733 | 0.229 | 0.098 | 0.567 |
| Rare confident FP c=0.01 m=0.10 | rare | 0.774 | 0.749 | 0.267 | 0.177 | 0.098 | 0.312 |
| Rare confident FP c=0.02 m=0.02 | rare | 0.981 | 0.981 | 0.917 | 0.516 | 0.296 | 0.827 |
| Rare confident FP c=0.02 m=0.10 | rare | 0.889 | 0.865 | 0.512 | 0.362 | 0.243 | 0.567 |
| Rare confident mixed c=0.01 m=0.02 | rare | 0.854 | 0.851 | 0.579 | 0.163 | 0.072 | 0.324 |
| Rare confident mixed c=0.01 m=0.10 | rare | 0.558 | 0.558 | 0.207 | 0.124 | 0.060 | 0.133 |
| Rare confident mixed c=0.02 m=0.02 | rare | 0.955 | 0.956 | 0.786 | 0.377 | 0.141 | 0.698 |
| Rare confident mixed c=0.02 m=0.10 | rare | 0.789 | 0.782 | 0.367 | 0.231 | 0.114 | 0.345 |

**C2 (errores raros confiados) = PASA**  (10 de 12 escenarios con ventaja >= 0.15 sobre Spiegelhalter; se exigen 8)

**C3 (componente FP) = PASA**  (escenarios FP con scan_fp - scan_tot >= 0.05: 1.000)
(Informativo, sin criterio: escenarios FN con scan_fp - scan_tot >= 0.05: 0.000)

**C6 (coste fuera de FP, informativo) = PASA**  (media scan_fp - scan_tot en escenarios no FP = 0.038; se tolera >= -0.03)
Referencia: media scan_fp - scan_both en todos los escenarios = -0.005

**C4 (informativo) global:** scan_fp - Spiegelhalter = Global overconfidence k=0.50: -0.013; Global overconfidence k=0.75: -0.088

## C5. Predicciones de un modelo ajustado (tipo I)

| n | ajustadas en la misma muestra | validacion cruzada 5-fold | repeticiones |
|---|---|---|---|
| 250 | 0.017 | 0.057 | 300 |
| 500 | 0.013 | 0.033 | 300 |

## Veredicto segun PREREGISTRATION.md
C1, C2 y C3 PASAN: test valido con ventaja demostrada en errores raros y con componente de falsos positivos. FN: sin afirmacion.
