# EBEC-scan: informe de validacion v3 (scan_fp) (modo lowprev)

## C1. Error tipo I bajo H0

| test | tasa agregada | min celda | max celda | celdas | liberales | conservadoras |
|---|---|---|---|---|---|---|
| scan_fp | 0.043 | 0.035 | 0.052 | 3 | 0 | 1 |
| scan_both | 0.043 | 0.034 | 0.052 | 3 | 0 | 1 |
| scan_tot | 0.047 | 0.042 | 0.054 | 3 | 0 | 0 |
| spiegelhalter | 0.053 | 0.050 | 0.055 | 3 | 0 | 0 |
| hl | 0.129 | 0.094 | 0.189 | 3 | 3 | 0 |
| cox | 0.053 | 0.050 | 0.055 | 3 | 0 | 0 |

**C1 (tamano) = PASA**  (tasa 0.0427; limite superior 0.0619; celdas liberales 0 de 3)

## Potencia por escenario (media de las celdas ejecutadas)

| escenario | grupo | scan_fp | scan_both | scan_tot | Spiegelhalter | HL | Cox |
|---|---|---|---|---|---|---|---|
| FN under-prediction c=0.10 d=3.0 | fn_under | 0.157 | 0.155 | 0.242 | 0.143 | 0.730 | 0.237 |
| FP over-prediction c=0.10 d=3.0 | fp_over | 0.217 | 0.212 | 0.007 | 0.000 | 0.198 | 0.413 |
| Global overconfidence k=0.50 | global | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Inversion c=0.10 s=1.0 | inv | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Rare confident FN c=0.01 m=0.02 | rare | 0.068 | 0.068 | 0.075 | 0.078 | 0.128 | 0.050 |
| Rare confident FN c=0.01 m=0.10 | rare | 0.050 | 0.050 | 0.045 | 0.052 | 0.123 | 0.070 |
| Rare confident FN c=0.02 m=0.02 | rare | 0.050 | 0.053 | 0.058 | 0.062 | 0.107 | 0.055 |
| Rare confident FN c=0.02 m=0.10 | rare | 0.027 | 0.028 | 0.028 | 0.028 | 0.115 | 0.062 |
| Rare confident FP c=0.01 m=0.02 | rare | 1.000 | 1.000 | 0.685 | 0.678 | 0.295 | 1.000 |
| Rare confident FP c=0.01 m=0.10 | rare | 1.000 | 1.000 | 0.472 | 0.515 | 0.287 | 0.970 |
| Rare confident FP c=0.02 m=0.02 | rare | 1.000 | 1.000 | 0.978 | 0.983 | 0.730 | 1.000 |
| Rare confident FP c=0.02 m=0.10 | rare | 1.000 | 1.000 | 0.900 | 0.918 | 0.643 | 1.000 |
| Rare confident mixed c=0.01 m=0.02 | rare | 0.988 | 0.988 | 0.372 | 0.343 | 0.152 | 0.855 |
| Rare confident mixed c=0.01 m=0.10 | rare | 0.948 | 0.948 | 0.182 | 0.220 | 0.157 | 0.627 |
| Rare confident mixed c=0.02 m=0.02 | rare | 1.000 | 1.000 | 0.660 | 0.678 | 0.327 | 0.992 |
| Rare confident mixed c=0.02 m=0.10 | rare | 0.998 | 0.998 | 0.418 | 0.467 | 0.298 | 0.910 |

**C2 (errores raros confiados) = FALLA**  (6 de 12 escenarios con ventaja >= 0.15 sobre Spiegelhalter; se exigen 8)

**C3 (componente FP) = PASA**  (escenarios FP con scan_fp - scan_tot >= 0.05: 0.800)
(Informativo, sin criterio: escenarios FN con scan_fp - scan_tot >= 0.05: 0.000)

**C6 (coste fuera de FP, informativo) = PASA**  (media scan_fp - scan_tot en escenarios no FP = 0.201; se tolera >= -0.03)
Referencia: media scan_fp - scan_both en todos los escenarios = 0.000

**C4 (informativo) global:** scan_fp - Spiegelhalter = Global overconfidence k=0.50: 0.000

## Veredicto segun PREREGISTRATION.md
C1 PASA, C2 FALLA: test valido pero sin ventaja de potencia demostrada; presentar como diagnostico descriptivo.
