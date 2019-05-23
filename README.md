# fdslrmAllinOne
All functions from R package fdslrm in one file. 
Functions are designed for modeling and prediction of time series using linear mixed models.

Installation and requirements
============================= 

Functions of fdslrm R package were developed with R version 3.4.2 (2017-09-28).

## Option 1: Clone with Git

To clone the fdslrm repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/gajdosandrej/fdslrmAllinOne.git
```
and you need to open fdslrmAllinOne.Rproj. 
Finally you load the functions from `fdslrmAllinOne.R` by running the following command in R console 
```
source("fdslrmAllinOne.R")
```
and the functions should be ready to use. 

## Option 2: Download R project from GitHub

On the GitHub page click on Clone or download and choose Donwload ZIP. 
Unzip the R project and open fdslrmAllinOne.Rproj in your RStudio. 
Consequently run the following command to load all R functions into your Environment:
```
source("fdslrmAllinOne.R")
```

## Option 3: Load R functions directly from GitHub

Just run the next command in R:
```
devtools::source_url("https://github.com/gajdosandrej/fdslrmAllinOne/blob/master/fdslrmAllinOne.R?raw=TRUE")
```

Acknowledgements
=======

This research is partially supported by the projects APVV-17-0568, VEGA 1/0311/18 and VVGS-PF-2018-792.

License
=======

See `LICENSE` for fdslrm licensing information.
