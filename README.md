# RLR_transform

**Research scripts on how transforming compositional time series changes
the causal links that Convergent Cross Mapping (CCM) infers.**

Relative abundances, such as microbial community data, are compositional: the
parts of each sample sum to a constant. These R scripts compare CCM on raw,
compositional and RLR-transformed (log ratio) series, using simulated
dynamics and public predator–prey data, and draw the figures of that
analysis.

> [!NOTE]
> These are figure scripts, not an R package. A small executable check runs
> in CI, but the full figure analyses have not been rerun, and the check does
> not establish causal findings.

Start with [Environment](#environment) and [Quick check](#quick-check).
Use [Running the figures](#running-the-figures) for the full analysis.

## Contents

| Script | Input | Additional packages |
| --- | --- | --- |
| `Fig1.R` | coupled synthetic dynamics generated in the script | deSolve, rgr, tidyverse, tseries |
| `Fig2.R` | synthetic dynamics generated in the script | tidyverse |
| `Fig3.R` | synthetic dynamics generated in the script | deSolve, tidyverse, tseries, DT, ggsci |
| `Fig4.R` | rEDM `paramecium_didinium` data | tidyverse, plotly, tseries |
| `Fig5.R` | included `Fig5.csv` | tidyverse, plotly, tseries |
| `FigS2.R` | rEDM `paramecium_didinium` data | tseries, DT, tidyverse |
| `smoke.R` | synthetic data | rEDM 0.7.3 only |

## Environment

The scripts use the historical rEDM interface (`simplex()`, `ccm()`,
`ccm_means()`). They run with **R 4.5.1 and rEDM 0.7.3**; current rEDM 2.x
has a different interface and cannot be used without a reviewed migration.
The original analysis environment was not recorded.

Install rEDM 0.7.3 into a separate library (a source package: Rtools on
Windows, a compiler on Linux and macOS):

```r
dir.create(".r-library", showWarnings = FALSE)
.libPaths(c(normalizePath(".r-library"), .libPaths()))
install.packages(c("Rcpp", "RcppEigen"), repos = "https://cloud.r-project.org")
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/rEDM/rEDM_0.7.3.tar.gz",
  repos = NULL, type = "source"
)
```

## Quick check

From the repository folder:

```sh
R_LIBS_USER="$PWD/.r-library" Rscript smoke.R
```

On PowerShell, set `$env:R_LIBS_USER=(Resolve-Path .r-library).Path` first.
The check parses all six figure scripts, takes the transformations and the
embedding and CCM functions from `Fig1.R`, and runs CCM on 120 synthetic
observations (three library sizes, ten samples each, fixed seed) for raw,
compositional and RLR series. Set `RLR_SMOKE_OUTPUT` to save the table. It
checks that the code runs; it is not a new scientific result.

## Running the figures

Run from the repository folder so that `Fig5.csv` is found, normally `Fig1.R`
to `Fig5.R` and then `FigS2.R`, after installing the packages in the table.
Use a separate copy of the folder: `Fig3.R` writes CSV and EPS files into the
working folder. Several scripts repeat surrogate analyses and take much longer
than `smoke.R`. Their stochastic results have not been compared with earlier
outputs.

RLR takes logarithms of ratios, so zero or negative counts and near-zero
denominators need a deliberate treatment; no pseudocount is imposed here.

## Validation

On 2026-10-03 the small CCM check passed on Windows. GitHub Actions runs the
same check on Linux with rEDM 0.7.3 pinned; the full figure dependencies are
outside that job.

## References

- [rEDM archive on CRAN](https://cran.r-project.org/src/contrib/Archive/rEDM/)
- [Original `ccm_means` documentation](https://ha0ye.github.io/rEDM/reference/ccm_means.html)
- [Current rEDM](https://github.com/SugiharaLab/rEDM)

## License

MIT
