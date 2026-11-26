# EBEC: Excess Brier Error Concentration

Test for whether the worst-predicted cases of a binary model accumulate
more Brier loss than expected if the predicted probabilities were true.

## Install

```r
# install.packages('remotes')
remotes::install_github('Gabrigle/EBEC')
```

## Example

```r
library(EBEC)
set.seed(1)
probs  <- runif(200, 0.05, 0.60)   # predicted probabilities
y_real <- rbinom(200, 1, probs)    # observed outcomes (0/1)

ebec_result <- EBEC::ebec(y_real, probs, B = 9999, classes = "none", seed = 42)
print(ebec_result, explain = FALSE, lang = "es")
```

## Citation

Iglesias-Caride, G. (2026). EBEC: Excess Brier Error Concentration (R package version 0.1.0). Zenodo. [DOI pending]
