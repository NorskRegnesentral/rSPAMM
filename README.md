<img src="R-package/man/figures/NR-logo_utvidet_r32g60b136_small.png" align="right" height="50px"/>

# The rSPAMM package

#### T. A. Øigård and M. Biuw

# Overview
The *rSPAMM* package is a R-package for performing assessment of harp and hooded seals in the West Ice (along the coast of Greenland) and the East Ice. Assessment of other populations is possible as long as similar input data is available. The package fits available pup production data and models the total abundance of the population of interest. Various catch options and the impact of the future abundance can be explored and graphical tools to visualize the results is available. 

This repository consists of two parts. One is the available data sets in the *Seal-data* folder and the package is found in the *R-package* folder.

Instructions on how to use *rSPAMM* is found in the *README.md* file in the *R-package* folder.


# Installation

The most recent version of *rSPAMM* is hosted on a git repository at
<https://github.com/NorskRegnesentral/rSPAMM.git>.

``` r
devtools::install_git(url = "https://github.com/NorskRegnesentral/rSPAMM.git", subdir = "R-package")
``` 
