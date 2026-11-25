# EBEC-scan: informe de validacion v2 (modo dry)

## C1. Error tipo I bajo H0

| test | tasa agregada | min celda | max celda | celdas | liberales | conservadoras |
|---|---|---|---|---|---|---|
| scan_all | 0.000 | 0.000 | 0.000 | 1 | 0 | 0 |
| scan_tot | 0.000 | 0.000 | 0.000 | 1 | 0 | 0 |
| spiegelhalter | 0.000 | 0.000 | 0.000 | 1 | 0 | 0 |
| hl | 0.050 | 0.050 | 0.050 | 1 | 0 | 0 |
| cox | 0.000 | 0.000 | 0.000 | 1 | 0 | 0 |

**C1 (tamano) = PASA**  (tasa 0.0000; limite superior 0.1962; celdas liberales 0 de 1)

## Potencia por escenario (media de las celdas ejecutadas)

| escenario | grupo | scan_all | scan_tot | Spiegelhalter | HL | Cox |
|---|---|---|---|---|---|---|
| FN under-prediction c=0.05 d=1.5 | fn_under | 0.000 | 0.000 | 0.200 | 0.200 | 0.100 |
| Inversion c=0.05 s=0.5 | inv | 0.100 | 0.100 | 0.200 | 0.100 | 0.000 |
| Rare confident FN c=0.02 m=0.10 | rare | 0.000 | 0.000 | 0.000 | 0.100 | 0.200 |

**C2 (errores raros confiados) = FALLA**  (0 de 1 escenarios con ventaja >= 0.15 sobre Spiegelhalter; se exigen 1)

**C3 (componentes por clase) = NO EVALUABLE**  (FN: 0.000 de escenarios con ventaja >= 0.05; FP: NA)

## C5. Predicciones de un modelo ajustado (tipo I)

| n | ajustadas en la misma muestra | validacion cruzada 5-fold | repeticiones |
|---|---|---|---|
| 250 | 0.000 | 0.000 | 10 |

## Veredicto segun PREREGISTRATION.md
C1 PASA, C2 FALLA: test valido pero sin ventaja de potencia demostrada; presentar como diagnostico descriptivo.
(Modo dry: PRUEBA TECNICA del codigo, sin significado estadistico.)
