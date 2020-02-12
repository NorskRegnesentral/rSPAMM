<img src="vignettes/figures/NR-logo_utvidet_r32g60b136_small.png" align="right" height="50px"/>

# The rSPAMM package

#### T. A. Øigård and M. Biuw

# Overview
The *rSPAMM* package is a R-package for performing assessment of harp and hooded seals in the West Ice (along the coast of Greenland) and the East Ice. Assessment of other populations is possible as long as similar input data is available. The package fits available pup production data and models the total abundance of the population of interest. Various catch options and the impact of the future abundance can be explored and graphical tools to visualize the results is available. 


Instructions on how to use *rSPAMM* is found in the vignette.


# Installation

The rSPAMM package depends on the following R-packages

  - *TMB*
  - *ggplot2*
  - *RcppEigen*
  
These dependencies will be automatically installed when installing the rSPAMM package.

The most recent version of *rSPAMM* is hosted on a git repository at
<https://github.com/NorskRegnesentral/rSPAMM.git>.


To install the R-package directly from the repository use the following command (note: the R-package devtools has to be installed first)
``` r
devtools::install_github("https://github.com/NorskRegnesentral/rSPAMM.git")
``` 

In order to load the package type:
```{r}
library(rSPAMM)
```

Instructions on how to use *rSPAMM* is found in the vignette.

To load the vignette type:
```{r}
vignette("howToUse_rSPAMM",package = "rSPAMM")
```
