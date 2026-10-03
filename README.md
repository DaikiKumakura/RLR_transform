# RLR_transform

Research scripts for compositional time series and CCM (Convergent Cross Mapping). These are figure scripts, not an R package. A small executable check is included; it does not reproduce all published figures or establish causal findings.

## Compatible environment

The figure scripts call the historical `simplex()`, `ccm()` and `ccm_means()` interfaces. The verified small example uses **R 4.5.1 and rEDM 0.7.3**. Current rEDM 2.x uses a different API and cannot be substituted into these scripts without a reviewed migration. The original analysis environment was not recorded; the verified environment below is not claimed to be that original environment.

Use a separate library. Installing a source package requires Rtools on Windows, or an R-compatible compiler on Linux/macOS.

```r
dir.create(".r-library", showWarnings = FALSE)
.libPaths(c(normalizePath(".r-library"), .libPaths()))
install.packages(c("Rcpp", "RcppEigen"), repos = "https://cloud.r-project.org")
install.packages("https://cran.r-project.org/src/contrib/Archive/rEDM/rEDM_0.7.3.tar.gz", repos = NULL, type = "source")
```

From this repository directory, run:

```sh
R_LIBS_USER="$PWD/.r-library" Rscript smoke.R
```

On PowerShell, set `$env:R_LIBS_USER=(Resolve-Path .r-library).Path` before `Rscript smoke.R`.

The check parses all six figure scripts, extracts the original transformations and embedding/CCM functions from Fig1.R, and runs 120 retained synthetic observations with three library sizes and ten samples each. An explicit seed is supplied to CCM. It prints results for raw, compositional and RLR series. Set `RLR_SMOKE_OUTPUT` to save the table. Comparing these small runs is a runtime check, not a new scientific result.

## Figure inputs and dependencies

| Script | Input | Additional packages |
| --- | --- | --- |
| Fig1.R | Coupled synthetic dynamics generated in the script | deSolve, rgr, tidyverse, tseries |
| Fig2.R | Synthetic dynamics generated in the script | tidyverse |
| Fig3.R | Synthetic dynamics generated in the script | deSolve, tidyverse, tseries, DT, ggsci |
| Fig4.R | rEDM `paramecium_didinium` data | tidyverse, plotly, tseries |
| Fig5.R | Included Fig5.csv | tidyverse, plotly, tseries |
| FigS2.R | rEDM `paramecium_didinium` data | tseries, DT, tidyverse |

Run from the repository directory so Fig5.csv resolves. The scripts can be run individually, normally Fig1.R through Fig5.R, then FigS2.R. Full runs require the additional packages above. Run those scripts in a separate copy because Fig3 writes named CSV/EPS outputs into the working directory. Do not overwrite prior results. Several scripts perform repeated surrogate analyses and can be substantially slower than smoke.R. Their stochastic figure results and numerical agreement with previous outputs have not been checked.

RLR uses logarithms of positive ratios and has singularities where its denominator is zero. Zero/nonpositive counts or near-singular ratios require a deliberate analytical treatment; no pseudocount policy is imposed here.

## Validation and maintenance

On 2026-10-03 the original-API small CCM check passed on Windows. GitHub Actions runs the same check in a clean Linux environment, with the legacy rEDM version pinned. Full figure dependencies are deliberately outside that lightweight job. Changes to the historical scientific scripts should be assessed separately from maintenance of smoke.R.

References: [CRAN archive](https://cran.r-project.org/src/contrib/Archive/rEDM/), [original ccm_means documentation](https://ha0ye.github.io/rEDM/reference/ccm_means.html), [current upstream rEDM](https://github.com/SugiharaLab/rEDM).
