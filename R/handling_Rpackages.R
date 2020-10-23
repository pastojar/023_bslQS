# script for handling R packages
# J. Pastorek


#  devtools::install_github("hadley/devtools")   # accesses latest devtools functions


# loading a package
myPack_load <- function( path = getwd() ) {
  
  if ( path != getwd() )   { setwd( path ) }
  
  if ( !require(devtools) )    { install.packages("devtools");    library(devtools) } 
  
  pack <- devtools::load_all()
  
  return( environmentName( pack$env ) )  
}



# installing a package to library
myPack_inst <- function( path = getwd() ) {
  
  if ( path != getwd() )   { setwd( path ) } # specify package folder
  
  if ( !require(devtools) )    { install.packages("devtools");    library(devtools) }
  
  devtools::install()
}



# creating a package (".gz" archive)
myPack_make_gz <- function( path = getwd() ) {
  
  if ( path != getwd() )   { setwd( path ) } # specify package folder
  
  if ( !require(devtools) )    { install.packages("devtools");    library(devtools) }
  if ( !require(roxygen2) )    { install.packages("roxygen2");    library(roxygen2) }
  if ( !require(testthat) )    { install.packages("testthat");    library(testthat) }
  if ( !require(knitr   ) )    { install.packages("knitr"   );    library(knitr) }
  
  devtools::build()
}