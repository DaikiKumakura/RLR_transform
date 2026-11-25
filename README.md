# RLR_transform

**Research scripts on how transforming compositional time series changes
the causal links that Convergent Cross Mapping (CCM) infers.**

## Installation

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

## Quick start

From the repository folder:

```sh
R_LIBS_USER="$PWD/.r-library" Rscript smoke.R
```

## Documentation

[Usage guide](https://github.com/DaikiKumakura/RLR_transform/blob/main/README-GUIDE.md)

## Validation and limitations

- **Validation:** Ran the small CCM check on Windows and configured the same check on Linux with rEDM 0.7.3.
- **Results:** The Windows smoke check passed.
- **Limitations:** The full figure workflow and original analysis environment have not been reproduced.
- **Reproduction:** [smoke.R](https://github.com/DaikiKumakura/RLR_transform/blob/main/smoke.R).

## Citation

Kumakura, D. *RLR_transform* [Software]. [GitHub repository](https://github.com/DaikiKumakura/RLR_transform).

Related references:

- [rEDM archive on CRAN](https://cran.r-project.org/src/contrib/Archive/rEDM/)
- [Original `ccm_means` documentation](https://ha0ye.github.io/rEDM/reference/ccm_means.html)
- [Current rEDM](https://github.com/SugiharaLab/rEDM)

## License

This project is licensed under MIT. See [LICENSE](LICENSE) for the license terms.
