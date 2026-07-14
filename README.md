# Assessing Covariate-Adjusted Risk Differences in Clinical Trials with Small Samples

This repository contains all code required to reproduce the simulation study and figures presented in the manuscript.

The repository provides a fully reproducible workflow, including data generation, implementation of statistical methods, simulation execution, aggregation of performance measures, and generation of manuscript figures.

All analyses are simulation-based and do not require any external or proprietary datasets.

---

# Repository Structure

```text
binary-endpoints/
│
├── reproducible_simulation.qmd
├── reproducible_simulation.html
│
├── figures/
│   └── Output directory containing all manuscript figures
│
├── R/
│   └── Supporting R functions for
│       - data generation
│       - implementation of statistical methods
│       - simulation execution
│       - analysis and summarization of simulation results
│
└── src/
    └── Supporting C++ code implementing an efficient
        bisection-search algorithm used within the simulation framework
```

---

# Master File

## `reproducible_simulation.qmd`

This document contains the complete simulation workflow used to generate the results reported in the manuscript.

Specifically, it:

1. loads all required packages and supporting source files,
2. defines simulation scenarios and data-generating mechanisms,
3. executes the simulation study,
4. summarizes simulation results,
5. computes performance measures,
6. generates all manuscript figures,
7. documents the computational environment used for the analyses.

---

# Reproducing the Results

From the repository root directory, render

```bash
quarto render reproducible_simulation.qmd
```

This will:

- execute the simulation study,
- generate the manuscript figures,
- save figure files to the `figures/` directory,
- render the supplementary computational document.

Depending on the available computing resources, execution may require substantial computation time.

---

# Manuscript Figures

All figures reported in the manuscript are generated automatically during execution of

```text
reproducible_simulation.qmd
```

and are saved to

```text
figures/
```

in PDF format.

No manual post-processing is required.

---

# Data Availability

The analyses are entirely simulation-based.

No patient-level data, proprietary datasets, or external data sources are required to reproduce the results.

All data used in the analyses are generated directly from code contained in this repository.

---

# Random Number Generation

Simulation results rely on pseudorandom number generation.

Random-number generator seeds are explicitly specified within the simulation workflow to facilitate reproducibility.

---

# Software Requirements

The simulation study was developed and executed using

- R version 4.5.2 (2025-10-31)
- Ubuntu 24.04.4 LTS (x86_64-pc-linux-gnu)

Supporting source files are loaded automatically from:

```text
R/
src/
```

and therefore must be present in the repository root directory.


## Computational Requirements

The simulation workflow was developed for execution on a high-performance computing (HPC) cluster using the `clustermq` package.

The default configuration assumes access to a large number of compute cores (more than 1000 concurrent jobs were available for the original analyses).

Users reproducing the analyses in a different environment will need to modify the cluster configuration and adjust the number of parallel jobs accordingly. Depending on available computational resources, this may substantially increase execution time.

In particular, users without access to an HPC environment should review the settings in the simulation section governing:

- cluster registration,
- number of jobs,
- number of parallel workers,
- memory allocation.

The simulation code can be executed with fewer cores, but runtime may increase considerably.


# Session Information

```text
R version 4.5.2 (2025-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.4 LTS

Matrix products: default

BLAS:
  /opt/R/openval/2026.03.00/4.5.2/lib/R/lib/libRblas.so

LAPACK:
  /opt/R/openval/2026.03.00/4.5.2/lib/R/lib/libRlapack.so

LAPACK version 3.12.1
```

## Random Number Generation

```text
RNG:     L'Ecuyer-CMRG
Normal:  Inversion
Sample:  Rejection
```

Random seeds are explicitly specified within the simulation workflow.

## Locale

```text
LC_CTYPE=en_US.UTF-8
LC_NUMERIC=C
LC_TIME=en_US.UTF-8
LC_COLLATE=en_US.UTF-8
LC_MONETARY=en_US.UTF-8
LC_MESSAGES=en_US.UTF-8
LC_PAPER=en_US.UTF-8
LC_NAME=C
LC_ADDRESS=C
LC_TELEPHONE=C
LC_MEASUREMENT=en_US.UTF-8
LC_IDENTIFICATION=C
```

```text
time zone: Etc/UTC
tzcode source: system (glibc)
```

## Attached Base Packages

```text
grid
parallel
stats
graphics
grDevices
utils
datasets
methods
base
```

## Attached Packages

```text
logistf_1.26.1
detectseparation_0.3
marginaleffects_0.25.1
margins_0.3.28
sandwich_3.1-1
Exact_3.3
RobinCar_1.0.0
kableExtra_1.4.0
doRNG_1.8.6.2
rngtools_1.5.2
foreach_1.5.2
clustermq_0.9.9
svglite_2.2.2
patchwork_1.3.2
ggplot2_4.0.1
DT_0.34.0
Rcpp_1.1.0
tibble_3.3.0
furrr_0.3.1
future_1.67.0-9017
purrr_1.0.4
tidyr_1.3.1
dplyr_1.1.4
```

## Additional Loaded Namespaces

```text
Rdpack_2.6.4
lpSolveAPI_5.5.2.0-17.14
rlang_1.1.6
magrittr_2.0.4
prediction_0.3.18
compiler_4.5.2
mgcv_1.9-4
systemfonts_1.3.1
vctrs_0.6.5
stringr_1.6.0
crayon_1.5.3
pkgconfig_2.0.3
shape_1.4.6.1
fastmap_1.2.0
backports_1.5.0
labeling_0.4.3
ROI_1.0-1
rmarkdown_2.30
nloptr_2.2.1
ragg_1.5.0
xfun_0.54
glmnet_4.1-8
jomo_2.7-6
AIPW_0.6.9.2
ROI.plugin.lpsolve_1.0-2
jsonlite_2.0.0
progress_1.2.3
pan_1.9
prettyunits_1.2.0
broom_1.0.10
R6_2.6.1
stringi_1.8.7
RColorBrewer_1.1-3
Rsolnp_1.16
rpart_4.1.24
boot_1.3-32
parallelly_1.44.0
estimability_1.5.1
numDeriv_2016.8-1.1
iterators_1.0.14
knitr_1.50
future.apply_1.11.3
zoo_1.8-14
nnet_7.3-20
Matrix_1.7-4
nnls_1.6
splines_4.5.2
tidyselect_1.2.1
emulator_1.2-24
rstudioapi_0.17.1
yaml_2.3.11
codetools_0.2-20
listenv_0.9.1
lattice_0.22-7
withr_3.0.2
S7_0.2.1
coda_0.19-4.1
evaluate_1.0.5
SuperLearner_2.0-29
fastDummies_1.7.5
survival_3.8-3
gam_1.22-5
xml2_1.5.1
pillar_1.11.1
mice_3.17.0
checkmate_2.3.3
reformulas_0.4.2
generics_0.1.4
truncnorm_1.0-9
hms_1.1.4
scales_1.4.0
rootSolve_1.8.2.4
minqa_1.2.8
xtable_1.8-4
globals_0.18.0
glue_1.8.0
slam_0.1-55
emmeans_1.11.1
tools_4.5.2
data.table_1.17.8
lme4_1.1-37
registry_0.5-1
mvtnorm_1.3-3
rbibutils_2.4
tidyverse_2.0.0
nlme_3.1-168
formula.tools_1.7.1
cli_3.6.5
textshaping_1.0.4
viridisLite_0.4.2
gtable_0.3.6
digest_0.6.39
operator.tools_1.6.3
progressr_0.15.1
htmlwidgets_1.6.4
farver_2.1.2
htmltools_0.5.8.1
lifecycle_1.0.4
mitml_0.4-5
MASS_7.3-65
```
