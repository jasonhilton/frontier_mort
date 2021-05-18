README
================

# Modelling Frontier Mortality using Bayesian Generalised Additive Models

This repository contains the code need to reproduce the results in the
paper “Modelling Frontier Mortality using Bayesian Generalised Additive
Models”. The paper is available from this repo at:
<https://github.com/jasonhilton/frontier_mort/raw/master/paper/mort_pacer.pdf>.
The appendix is similarly at:
<https://github.com/jasonhilton/frontier_mort/raw/master/paper/appendix.pdf>

**Link to Presentation:**
<https://jasonhilton.github.io/frontier_mort/presentation/LSHTM_mortality.html>

# Requirements

Running the code requires a recent version of `R` (\>3.5 is probably
best) and of `rstan`, the `R` interface to stan ([see
here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

Although it is not strictly necessary, a High Performance Computing
platform was used to minimise runtime. The portable batch system (PBS)
scripts needed to conduct the runs on this platform are included, and it
should be possible to run these on similar systems with minimal
adaptation.

Processing the results is memory intensive. The systems used have a
minimum of 16gb of RAM. The scripts may run on systems with less RAM,
but this is untested. It is possible the process may hang as the
available ram is exhausted and swap space begins to be used. If you wish
to replicate the result but do not have access to a system with enough
RAM, let me know - I may be able to optimise the code to be less memory
intensive.

To produce the final manuscript, Rmarkdown, knitr, latex and pandoc are
used.

## Storage

The complete set of MCMC results take up less than 1GB.

-----

# Reproducing the results.

The results presented in the paper can be reproduced by following the
steps below. All code chucks in this document should be run from the
terminal from the base directory of this repository (that is, the one
that contains this file).

## Installing R packages

  - The first step is to install a number of R packages. The command
    below does so programmatically, and can be run from the command
    line.

<!-- end list -->

``` bash
Rscript install_packages.R
```

## Getting the data

UK mortality data from the human mortality database is required. This
data can also be obtained using an script. Running this script will
create a `data` subfolder within this repository and download files from
the [Human Mortality Database](http://www.mortality.org/) (HMD) to it.
However, an account and password is needed to access this data. For this
script to work, you need to set environment variables containing your
username and password. You can do this by adding the following lines to
your `.Renviron.site` file (create one if you do not already have one),
which will be in the directory specified by running the command
`R.home()` in the R console, in the `etc` sub-folder.

``` bash
HFD_user=user
HFD_pass=pass
```

The ‘user’ and ‘pass’ text after the equals signs should be replaced
with your account username and password. Then, run the command below
from the command line to download the data.

``` bash
Rscript download_mortality_data.R
```

This script uses the third-party `HMDHFDplus` package to obtain the
data. Please ensure you use a unique password for HMD and don’t reuse
one from something else, and take appropriate precautions with your
credential information. I am not a security expert, and I cannot
guarantee the security of these credentials.

## Running the model

The script `scripts/run_model.R` produces samples for a single model.
Exactly what data and which model is used is determined by the
configuration file you specify. These are located in the `config`
folder, while the models themselves are located in the `stan` folder. To
run the models required the results presented in the paper, the
following commands are needed.

``` bash
Rscript run_model.R linear
Rscript run_model.R quadratic
Rscript run_model.R independent
```

You should allow two days for each model to reach the required 8000
samples per chain. Outputs are saved to the `results` directory, in
subfolders named according to the date and time the models were run.

## Running with Portable Batch System

If you wish to run the model on a cluster (assuming you have one
available), the following commands will submit jobs using PBS.

``` bash

qsub PBS/run_stan.pbs -v config=linear
qsub PBS/run_stan.pbs -v config=quadratic
qsub PBS/run_stan.pbs -v config=independent
```

This is convenient, because you can run all three models at the same
time (if you have enough nodes available), and you avoid tying up your
own machine for the duration\!

## Assessment

To assess the forecast root-mean-squared error and coverage of each
mode, another script must be run. The argument to this script should be
the name of the folder where the results for the model to be assessed
are saved. This folder will have been created by the script run in the
previous step, and will have a `date_time` type format:

``` bash

Rscript scripts/run_assess.R YYYYMMDD_HHMMSS
```

Again, this can be run on a cluster if necessary:

``` bash

qsub PBS/run_assess.pbs -v run_stamp=YYYYMMDD_HHMMSS
```

This script will save the assessment outputs in the same folder as the
original stan model output.

## Compile the paper

To produce the `tex` files needed to produce the manuscript, the below
commands run knitr and pandoc to process the R markdown files
`mort_pacer.Rmd` and `appendix.Rmd`. However, you need to specify the
folders where the relevant results are located by editing the two
`process_` scripts to include the relevant timestamps corresponding to
the linear, quadratic and independent models. These are included as
`params` to the `rmarkdown::render` call.

``` bash

Rscript scripts/process_rmd.R
Rscript scripts/process_appendix.R
```

## R session information

### Local Machine

``` r
library(tibble)
library(magrittr)
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     extract

``` r
library(purrr)
```

    ## 
    ## Attaching package: 'purrr'

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     set_names

``` r
library(lazyeval)
```

    ## 
    ## Attaching package: 'lazyeval'

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     is_atomic, is_formula

``` r
library(tidyselect)
library(rstan)
```

    ## Loading required package: StanHeaders

    ## Loading required package: ggplot2

    ## rstan (Version 2.21.2, GitRev: 2e1f913d3ca3)

    ## For execution on a local, multicore CPU with excess RAM we recommend calling
    ## options(mc.cores = parallel::detectCores()).
    ## To avoid recompilation of unchanged Stan programs, we recommend calling
    ## rstan_options(auto_write = TRUE)

    ## 
    ## Attaching package: 'rstan'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     extract

``` r
library(yaml)
library(here)
```

    ## here() starts at /storage/frontier_mort

``` r
library(git2r)
```

    ## 
    ## Attaching package: 'git2r'

    ## The following object is masked from 'package:rstan':
    ## 
    ##     lookup

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     is_empty, when

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     add

``` r
library(HMDHFDplus)
library(curl)
library(ggplot2)
library(ggfan)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:git2r':
    ## 
    ##     pull

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 3.6.3 (2020-02-29)
    ##  os       Ubuntu 18.04.5 LTS          
    ##  system   x86_64, linux-gnu           
    ##  ui       X11                         
    ##  language en_GB:en                    
    ##  collate  en_GB.UTF-8                 
    ##  ctype    en_GB.UTF-8                 
    ##  tz       Europe/London               
    ##  date     2020-10-28                  
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package      * version  date       lib source        
    ##  assertthat     0.2.1    2019-03-21 [1] CRAN (R 3.6.0)
    ##  backports      1.1.9    2020-08-24 [1] CRAN (R 3.6.3)
    ##  bitops         1.0-6    2013-08-17 [1] CRAN (R 3.6.0)
    ##  callr          3.4.3    2020-03-28 [1] CRAN (R 3.6.3)
    ##  cli            2.0.2    2020-02-28 [1] CRAN (R 3.6.3)
    ##  codetools      0.2-16   2018-12-24 [4] CRAN (R 3.6.3)
    ##  colorspace     1.4-1    2019-03-18 [1] CRAN (R 3.6.0)
    ##  crayon         1.3.4    2017-09-16 [1] CRAN (R 3.6.0)
    ##  curl         * 4.3      2019-12-02 [1] CRAN (R 3.6.2)
    ##  desc           1.2.0    2018-05-01 [1] CRAN (R 3.6.0)
    ##  devtools     * 2.3.1    2020-07-21 [1] CRAN (R 3.6.3)
    ##  digest         0.6.25   2020-02-23 [1] CRAN (R 3.6.3)
    ##  dplyr        * 1.0.2    2020-08-18 [1] CRAN (R 3.6.3)
    ##  ellipsis       0.3.1    2020-05-15 [1] CRAN (R 3.6.3)
    ##  evaluate       0.14     2019-05-28 [1] CRAN (R 3.6.1)
    ##  fansi          0.4.1    2020-01-08 [1] CRAN (R 3.6.2)
    ##  fs             1.5.0    2020-07-31 [1] CRAN (R 3.6.3)
    ##  generics       0.0.2    2018-11-29 [1] CRAN (R 3.6.3)
    ##  ggfan        * 0.1.3    2019-03-07 [1] CRAN (R 3.6.1)
    ##  ggplot2      * 3.3.2    2020-06-19 [1] CRAN (R 3.6.3)
    ##  git2r        * 0.27.1   2020-05-03 [1] CRAN (R 3.6.3)
    ##  glue           1.4.2    2020-08-27 [1] CRAN (R 3.6.3)
    ##  gridExtra      2.3      2017-09-09 [1] CRAN (R 3.6.0)
    ##  gtable         0.3.0    2019-03-25 [1] CRAN (R 3.6.0)
    ##  here         * 0.1      2017-05-28 [1] CRAN (R 3.6.1)
    ##  HMDHFDplus   * 1.9.13   2020-02-20 [1] CRAN (R 3.6.3)
    ##  htmltools      0.5.0    2020-06-16 [1] CRAN (R 3.6.3)
    ##  httr           1.4.2    2020-07-20 [1] CRAN (R 3.6.3)
    ##  inline         0.3.15   2018-05-18 [1] CRAN (R 3.6.0)
    ##  jsonlite       1.7.0    2020-06-25 [1] CRAN (R 3.6.3)
    ##  knitr          1.29     2020-06-23 [1] CRAN (R 3.6.3)
    ##  lazyeval     * 0.2.2    2019-03-15 [1] CRAN (R 3.6.1)
    ##  lifecycle      0.2.0    2020-03-06 [1] CRAN (R 3.6.3)
    ##  loo            2.3.1    2020-07-14 [1] CRAN (R 3.6.3)
    ##  magrittr     * 1.5      2014-11-22 [1] CRAN (R 3.6.1)
    ##  matrixStats    0.56.0   2020-03-13 [1] CRAN (R 3.6.3)
    ##  memoise        1.1.0    2017-04-21 [1] CRAN (R 3.6.0)
    ##  munsell        0.5.0    2018-06-12 [1] CRAN (R 3.6.0)
    ##  pillar         1.4.6    2020-07-10 [1] CRAN (R 3.6.3)
    ##  pkgbuild       1.1.0    2020-07-13 [1] CRAN (R 3.6.3)
    ##  pkgconfig      2.0.3    2019-09-22 [1] CRAN (R 3.6.2)
    ##  pkgload        1.1.0    2020-05-29 [1] CRAN (R 3.6.3)
    ##  prettyunits    1.1.1    2020-01-24 [1] CRAN (R 3.6.2)
    ##  processx       3.4.3    2020-07-05 [1] CRAN (R 3.6.3)
    ##  ps             1.3.4    2020-08-11 [1] CRAN (R 3.6.3)
    ##  purrr        * 0.3.4    2020-04-17 [1] CRAN (R 3.6.3)
    ##  R6             2.4.1    2019-11-12 [1] CRAN (R 3.6.2)
    ##  Rcpp           1.0.5    2020-07-06 [1] CRAN (R 3.6.3)
    ##  RcppParallel   5.0.1    2020-05-06 [1] CRAN (R 3.6.3)
    ##  RCurl          1.98-1.2 2020-04-18 [1] CRAN (R 3.6.3)
    ##  remotes        2.2.0    2020-07-21 [1] CRAN (R 3.6.3)
    ##  rlang          0.4.7    2020-07-09 [1] CRAN (R 3.6.3)
    ##  rmarkdown      2.3      2020-06-18 [1] CRAN (R 3.6.3)
    ##  rprojroot      1.3-2    2018-01-03 [1] CRAN (R 3.6.0)
    ##  rstan        * 2.21.2   2020-07-27 [1] CRAN (R 3.6.3)
    ##  scales         1.1.1    2020-05-11 [1] CRAN (R 3.6.3)
    ##  sessioninfo    1.1.1    2018-11-05 [1] CRAN (R 3.6.0)
    ##  StanHeaders  * 2.21.0-6 2020-08-16 [1] CRAN (R 3.6.3)
    ##  stringi        1.4.6    2020-02-17 [1] CRAN (R 3.6.3)
    ##  stringr        1.4.0    2019-02-10 [1] CRAN (R 3.6.0)
    ##  testthat       2.3.2    2020-03-02 [1] CRAN (R 3.6.3)
    ##  tibble       * 3.0.3    2020-07-10 [1] CRAN (R 3.6.3)
    ##  tidyr        * 1.1.2    2020-08-27 [1] CRAN (R 3.6.3)
    ##  tidyselect   * 1.1.0    2020-05-11 [1] CRAN (R 3.6.3)
    ##  usethis      * 1.6.1    2020-04-29 [1] CRAN (R 3.6.3)
    ##  V8             3.2.0    2020-06-19 [1] CRAN (R 3.6.3)
    ##  vctrs          0.3.4    2020-08-29 [1] CRAN (R 3.6.3)
    ##  withr          2.2.0    2020-04-20 [1] CRAN (R 3.6.3)
    ##  xfun           0.16     2020-07-24 [1] CRAN (R 3.6.3)
    ##  XML            3.99-0.3 2020-01-20 [1] CRAN (R 3.6.2)
    ##  yaml         * 2.2.1    2020-02-01 [1] CRAN (R 3.6.2)
    ## 
    ## [1] /home/jason/R/x86_64-pc-linux-gnu-library/3.6
    ## [2] /usr/local/lib/R/site-library
    ## [3] /usr/lib/R/site-library
    ## [4] /usr/lib/R/library

-----

### Iridis

``` r
print(readRDS(file.path("results", "hpc_sesh.Rds")))
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux Server release 6.10 (Santiago)
    ## 
    ## Matrix products: default
    ## BLAS:   /home/local/software/R/3.5.1/lib64/R/lib/libRblas.so
    ## LAPACK: /home/local/software/R/3.5.1/lib64/R/lib/libRlapack.so
    ## 
    ## Random number generation:
    ##  RNG:     
    ##  Normal:  
    ##  Sample:  
    ##  
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ## [1] splines   stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] rstan_2.19.2          ggplot2_3.2.1         StanHeaders_2.18.1-10
    ##  [4] git2r_0.26.1          yaml_2.2.0            stringr_1.4.0        
    ##  [7] stringi_1.4.3         readr_1.3.1           tibble_2.1.3         
    ## [10] tidyr_0.8.3           dplyr_0.8.3           purrr_0.3.2          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.2         pillar_1.4.2       compiler_3.5.1     prettyunits_1.0.2 
    ##  [5] tools_3.5.1        zeallot_0.1.0      pkgbuild_1.0.5     gtable_0.3.0      
    ##  [9] pkgconfig_2.0.2    rlang_0.4.0        cli_1.1.0          parallel_3.5.1    
    ## [13] loo_2.1.0          gridExtra_2.3      withr_2.1.2        vctrs_0.2.0       
    ## [17] hms_0.5.1          stats4_3.5.1       grid_3.5.1         tidyselect_0.2.5  
    ## [21] glue_1.3.1         inline_0.3.15      R6_2.4.0           processx_3.4.1    
    ## [25] callr_3.3.1        magrittr_1.5       matrixStats_0.54.0 backports_1.1.4   
    ## [29] scales_1.0.0       ps_1.3.0           assertthat_0.2.1   colorspace_1.4-1  
    ## [33] lazyeval_0.2.2     munsell_0.5.0      crayon_1.3.4
