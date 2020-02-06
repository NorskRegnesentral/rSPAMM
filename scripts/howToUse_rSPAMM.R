## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----setup, include=FALSE-----------------------------------------------------
library(rSPAMM)
library(TMB)

## ----eval = FALSE-------------------------------------------------------------
#  source("install_rSPAMM.R")

## ----eval = FALSE-------------------------------------------------------------
#  library(rSPAMM)
#  library(TMB)

## ----echo = FALSE-------------------------------------------------------------
catch_data <- read.table(paste("Data/harpeast/catch_data.dat",sep = ""),header = FALSE)
head(catch_data)

## ----echo = FALSE-------------------------------------------------------------
pup_production <- read.table(paste("Data/harpeast/pup_production.dat",sep = ""),header = FALSE)
pup_production

## ----echo = FALSE-------------------------------------------------------------
fecundity <- read.table(paste("Data/harpeast/fecundity.dat",sep = ""),header = FALSE)
fecundity

## ----echo = FALSE-------------------------------------------------------------
# Birth ogives for different periods
Pdat <- read.table(paste("Data/harpeast/wgharp.ogi",sep = ""),sep = "",header = TRUE)
# Which periods the various birth ogives applies to
Pper <- read.table(paste("Data/harpeast/wgharp.ogp",sep = ""),header = TRUE)
Pdat

## ----echo = FALSE-------------------------------------------------------------
Pper

## ----echo = FALSE-------------------------------------------------------------
priors <- read.table(paste("Data/harpeast/priors.dat",sep = ""),header = FALSE)
priors

## ----eval = FALSE-------------------------------------------------------------
#  data <- load.data(population = "harpeast")

## ----eval = TRUE,echo=FALSE---------------------------------------------------
library(rSPAMM)

#Set population
population = 'harpeast'
#Load the data
data <- load.data(population = population)

## ----eval = TRUE--------------------------------------------------------------
names(data)

## -----------------------------------------------------------------------------
data$Amax

data$pupProductionData


## ----eval = TRUE--------------------------------------------------------------
parameters <- load.initial.values(population = population)

## ----eval = TRUE--------------------------------------------------------------
parameters

## ----eval = FALSE-------------------------------------------------------------
#  obj = load.model.object(dat = data,par = parameters)

