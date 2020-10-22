# script for handling R packages
# J. Pastorek, 17.11.2015


#  devtools::install_github("hadley/devtools")   # accesses latest devtools functions


# loading a package

if ( !require(devtools) )    { install.packages("devtools");    library(devtools) } 
setwd("D:/packs/023_bslQS/023_bslQS_01_pilot/bslQS")
devtools::load_all()


# installing a package to library

setwd("C:/Users/xy/Documents/uni/000_eawag/normal.pack")  # specify package folder
if ( !require(devtools) )    { install.packages("devtools");    library(devtools) } 
devtools::install()


# creating a package (".gz" archive)

setwd("C:/Users/xy/Documents/uni/000_eawag/normal.pack/")  # specify package folder
if ( !require(devtools) )    { install.packages("devtools");    library(devtools) }
if ( !require(roxygen2) )    { install.packages("roxygen2");    library(roxygen2) } 
if ( !require(testthat) )    { install.packages("testthat");    library(testthat) } 
if ( !require(knitr   ) )    { install.packages("knitr"   );    library(knitr) } 

devtools::build()