# EBEC-scan: informe de validacion v2 (modo pilot)

## C1. Error tipo I bajo H0

| test | tasa agregada | min celda | max celda | celdas | liberales | conservadoras |
|---|---|---|---|---|---|---|
| scan_all | 0.046 | 0.030 | 0.065 | 4 | 0 | 0 |
| scan_tot | 0.050 | 0.030 | 0.060 | 4 | 0 | 0 |
| spiegelhalter | 0.049 | 0.035 | 0.065 | 4 | 0 | 0 |
| hl | 0.041 | 0.030 | 0.055 | 4 | 0 | 0 |
| cox | 0.048 | 0.030 | 0.060 | 4 | 0 | 0 |

**C1 (tamano) = PASA**  (tasa 0.0462; limite superior 0.0731; celdas liberales 0 de 4)

## Potencia por escenario (media de las celdas ejecutadas)

| escenario | grupo | scan_all | scan_tot | Spiegelhalter | HL | Cox |
|---|---|---|---|---|---|---|
| FN under-prediction c=0.05 d=1.5 | fn_under | 0.210 | 0.270 | 0.200 | 0.130 | 0.120 |
| FN under-prediction c=0.05 d=3.0 | fn_under | 0.860 | 0.930 | 0.860 | 0.740 | 0.760 |
| FN under-prediction c=0.10 d=1.5 | fn_under | 0.690 | 0.790 | 0.710 | 0.590 | 0.550 |
| FN under-prediction c=0.10 d=3.0 | fn_under | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| FP over-prediction c=0.05 d=1.5 | fp_over | 0.130 | 0.040 | 0.090 | 0.080 | 0.130 |
| FP over-prediction c=0.05 d=3.0 | fp_over | 0.320 | 0.040 | 0.030 | 0.140 | 0.290 |
| FP over-prediction c=0.10 d=1.5 | fp_over | 0.330 | 0.050 | 0.040 | 0.370 | 0.520 |
| FP over-prediction c=0.10 d=3.0 | fp_over | 0.810 | 0.050 | 0.070 | 0.690 | 0.820 |
| Global overconfidence k=0.50 | global | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Global overconfidence k=0.75 | global | 0.660 | 0.770 | 0.830 | 0.480 | 0.690 |
| Inversion c=0.05 s=0.5 | inv | 0.150 | 0.200 | 0.230 | 0.080 | 0.100 |
| Inversion c=0.05 s=1.0 | inv | 0.360 | 0.370 | 0.400 | 0.240 | 0.270 |
| Inversion c=0.10 s=0.5 | inv | 0.370 | 0.470 | 0.440 | 0.170 | 0.230 |
| Inversion c=0.10 s=1.0 | inv | 0.810 | 0.900 | 0.920 | 0.610 | 0.780 |
| Rare confident FN c=0.01 m=0.02 | rare | 0.260 | 0.340 | 0.050 | 0.040 | 0.070 |
| Rare confident FN c=0.01 m=0.10 | rare | 0.080 | 0.070 | 0.090 | 0.120 | 0.100 |
| Rare confident FN c=0.02 m=0.02 | rare | 0.590 | 0.640 | 0.160 | 0.140 | 0.170 |
| Rare confident FN c=0.02 m=0.10 | rare | 0.060 | 0.140 | 0.130 | 0.120 | 0.060 |
| Rare confident FP c=0.01 m=0.02 | rare | 0.960 | 0.900 | 0.140 | 0.030 | 0.520 |
| Rare confident FP c=0.01 m=0.10 | rare | 0.830 | 0.160 | 0.150 | 0.120 | 0.180 |
| Rare confident FP c=0.02 m=0.02 | rare | 1.000 | 1.000 | 0.480 | 0.290 | 0.950 |
| Rare confident FP c=0.02 m=0.10 | rare | 0.990 | 0.480 | 0.380 | 0.180 | 0.610 |
| Rare confident mixed c=0.01 m=0.02 | rare | 0.830 | 0.670 | 0.170 | 0.070 | 0.270 |
| Rare confident mixed c=0.01 m=0.10 | rare | 0.510 | 0.110 | 0.120 | 0.070 | 0.080 |
| Rare confident mixed c=0.02 m=0.02 | rare | 0.980 | 0.940 | 0.350 | 0.110 | 0.670 |
| Rare confident mixed c=0.02 m=0.10 | rare | 0.820 | 0.310 | 0.230 | 0.070 | 0.230 |

**C2 (errores raros confiados) = PASA**  (10 de 12 escenarios con ventaja >= 0.15 sobre Spiegelhalter; se exigen 8)

**C3 (componentes por clase) = FALLA**  (FN: 0.000 de escenarios con ventaja >= 0.05; FP: 0.875)

**C4 (informativo) global:** scan_all - Spiegelhalter = Global overconfidence k=0.50: 0.000; Global overconfidence k=0.75: -0.170

## C5. Predicciones de un modelo ajustado (tipo I)

| n | ajustadas en la misma muestra | validacion cruzada 5-fold | repeticiones |
|---|---|---|---|
| 250 | 0.000 | 0.050 | 100 |

## Veredicto segun PREREGISTRATION.md
C1 y C2 PASAN, C3 FALLA: ventaja en errores raros; NO afirmar atribucion FN/FP.
(Modo pilot: orientativo, NO usar en el paper.)
