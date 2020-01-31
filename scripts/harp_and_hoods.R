library(TMB)
source("functions.R")

# Lag funksjoner
# * Les inn data OK
# * Kjør modellen
# * Les inn data og kjør modellen
# * Finn likevektsfangsten
# * Lag figurer
# * Skriv ut resultater
# * Fortsett å lage flere funksjoner.

#compile("harp_and_hoods.cpp")
compile("harps_and_hoods_population_model.cpp","-O1 -g",DLLFLAGS="")
dyn.load(dynlib("harps_and_hoods_population_model"))

# Parameters -----------------
# Specify which population
which_population <- "harpeast"
# Number of years to run projections
years_of_prediction = 15
# Maximum age group
Amax = 20
# Annual catches in the projection
catch_quota = c(0,0)
# Set fecundity value for projections. Either set a fixed fecundity rate, "mean" which is the
# average fecundity rate observed or "NA" to use the last observed fecundity rate
Fproj = "mean"

# Read in data ---------------

# Catch data
catch_data <- read.table(paste("../Data/",which_population,"/catch_data.dat",sep = ""),header = FALSE)
# Pup production estimates
pup_production <- read.table(paste("../Data/",which_population,"/pup_production.dat",sep = ""),header = FALSE)
# Available fecundity data
fecundity <- read.table(paste("../Data/",which_population,"/fecundity.dat",sep = ""),header = FALSE)
# Birth ogives for different periods
Pdat <- read.table(paste("../Data/",which_population,"/wgharp.ogi",sep = ""),sep = "",header = TRUE)
# Which periods the various birth ogives applies to
Pper <- read.table(paste("../Data/",which_population,"/wgharp.ogp",sep = ""),header = TRUE)
#Pmatrix <- read.table("hessm.pma",header = FALSE)				#Birth ogives
#Priors used
priors <- read.table(paste("../Data/",which_population,"/priors.dat",sep = ""),header = FALSE)					#Priors used

FecAndP = buildPandF(Fdat = fecundity,cdata = catch_data,Fproj = Fproj,Pdat = Pdat,Pper = Pper)
Fdt = FecAndP$Fdt       #SJEKK VERDIEN PÅ DEN SISTE OG SAMMENLIKN MED WGHARP
Pmat = FecAndP$Pmatrix

# Prepare input to the model ------------
data <- list()
data$Amax = Amax													#Maximum age group
data$Cdata = as.matrix(catch_data)								#Add catch data
data$Nc = length(catch_data[,1])									#Number of years with catch data
data$pupProductionData = as.matrix(pup_production)					#Pup production estimates
data$Np = length(pup_production[,1])								#Number of years with pup production data
#data$Ftmp = as.matrix(fecundity)							#Observed fecundity rates
data$Ftmp = as.vector(Fdt)
#data$Nf = length(fecundity$V1)									#Number of years with fecundity data
data$Pmat = as.matrix(Pmat)									#Birt ogives
data$Npred = years_of_prediction									#Number of years to run projections
data$priors = as.matrix(priors)									#Priors for estimated parameters
data$Npriors = length(priors$V1)									#Number of priors
#data$CQuota = catch_quota											#Catch level in future projections


# Parameter to be estimated --------------------

#Initial values
Kinit = 2000000								#Initial population size
Minit = 0.09									#Natural adult mortality
M0init = 0.27								#Natural pup mortality

#Transform some parameters to ensure that estimates are > 0
parameters <- list()
parameters$logK= log(Kinit) #THESE SHOULD BE BOUNDED PARAMETERS, FIX THAT LATER
parameters$Mtilde= logit(Minit)
parameters$M0tilde= logit(M0init)


obj <- MakeADFun(data,parameters,DLL="harps_and_hoods_population_model")
#obj <- MakeADFun(data,parameters,random="u",DLL="hessm")

obj$fn()
obj$gr()

system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))

rep<-sdreport(obj, getJointPrecision=TRUE)
rep.matrix <- summary(rep)
rep.rnames = rownames(rep.matrix)
indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
indD1 = which(rep.rnames=="D1");
indD1New = which(rep.rnames=="D1New");   #NEED TO BE FIXED
indN0Current = which(rep.rnames=="N0CurrentYear");

yrs = 1946:2032

#Extract parameters
Kest = exp(opt$par[1])
Mest = ilogit(opt$par[2])
M0est = ilogit(opt$par[3])

D1 = rep.matrix[indD1,1]
D1New = rep.matrix[indD1New,1]
N0Current = rep.matrix[indN0Current,1]
D1.sd = rep.matrix[indD1,2]
D1New.sd = rep.matrix[indD1New,2]
N0Current.sd = rep.matrix[indN0Current,2]

#allsd<-sqrt(diag(solve(rep$jointPrecision)))
#plsd <- obj$env$parList(par=allsd)
quartz("",9,7)
plot(yrs,rep.matrix[indN1,1],type = "l",xlab = "Year",ylab = "Abundance",xlim = c(1946,2030),ylim = c(0,2000000),lwd = 3,col = "darkgreen")
lines(yrs,rep.matrix[indN0,1],col = "blue",lwd = 3)
lines(data$pupProductionData[,1],data$pupProductionData[,2],type = "p",col = "red",cex = 1.2,pch = 16)


