# EBEC: Excess Brier Error Concentration

Test for whether the worst-predicted cases of a binary model accumulate
more Brier loss than expected if the predicted probabilities were true.

## Install

```r
# install.packages('remotes')
remotes::install_github('[tu_usuario]/EBEC')
```

## Example

```r
library(EBEC)
set.seed(1); p <- runif(200, 0.2, 0.8); y <- rbinom(200, 1, p)
ebec(y, p, B = 199, ci = FALSE)
```

## Citation

Iglesias-Caride, G. (2026). EBEC: Excess Brier Error Concentration. Zenodo. [DOI]
