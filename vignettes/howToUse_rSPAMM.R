## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----eval = FALSE--------------------------------------------------------
## source("install_rSPAMM.R")


## ----eval = FALSE--------------------------------------------------------
## library(rSPAMM)


## ----eval = FALSE--------------------------------------------------------
## data <- load.data(population = "harpeast")


## ----eval = TRUE,echo=FALSE----------------------------------------------
library(rSPAMM)

#Set population
population = 'harpeast'
#Load the data
data <- load.data(population = population)


## ----eval = TRUE---------------------------------------------------------
names(data)


## ------------------------------------------------------------------------
data$Amax

data$pupProductionData



## ----eval = TRUE---------------------------------------------------------
parameters <- load.initial.values(population = population)


## ----eval = TRUE---------------------------------------------------------
parameters


## ----eval = FALSE--------------------------------------------------------
## obj = load.model.object(dat = data,par = parameters)

## ----eval = TRUE, echo = FALSE,include = FALSE---------------------------
obj = load.model.object(dat = data,par = parameters,template='harps_and_hoods_population_model2')


## ------------------------------------------------------------------------
names(obj)


## ------------------------------------------------------------------------
obj$fn()

obj$gr()


## ----eval = FALSE--------------------------------------------------------
## opt <- run.model(object = obj)


## ----echo = FALSE--------------------------------------------------------
opt <- run.model(object = obj)


## ----eval = FALSE--------------------------------------------------------
## res <- model.results(dat = data,object = obj,optimized = opt)

## ----echo = FALSE,include = FALSE----------------------------------------
res <- model.results(dat = data,object = obj,optimized = opt)


## ------------------------------------------------------------------------
names(res)


## ----eval = FALSE--------------------------------------------------------
## partab <- par.table(results=res, dat=data, tab2flex=FALSE)

## ----include = FALSE-----------------------------------------------------
partab <- par.table(results=res, dat=data, tab2flex=FALSE) 

## ------------------------------------------------------------------------
partab 


## ----cars----------------------------------------------------------------
summary(cars)


## ----pressure, echo=FALSE------------------------------------------------
plot(pressure)

