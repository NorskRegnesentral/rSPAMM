# Fitting state-space models to seal populations with scarce data
# Tor Arne Øigård and Hans J. Skaug
# ICES Journal of Marine Science, 2014
# Contact:  Tor Arne Øigård  (toro@imr.no) and Hans Julius Skaug (skaug@imr.no)

library(TMB)

#####################################
#Some functions needed below
ilogit1 <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
	return(log(x/(1-x)))
}

ilogit <- function(x){
	return(exp(x)/(1+exp(x)))
}
######################################

compile("hessm.cpp")
dyn.load(dynlib("hessm"))

###########
#Data part
###########
#Read in data
catch_data <- read.table("hessm.cat",header = FALSE)				#Catch data
pup_production <- read.table("hessm.est",header = FALSE)			#Pup production estimates
fecundity <- read.table("hessm.f",header = FALSE)				#Available fecundity data
Pmatrix <- read.table("hessm.pma",header = FALSE)				#Birth ogives
priors <- read.table("hessm.pri",header = FALSE)					#Priors used
years.of.prediction = 15											#Number of years to run projections

#Prepare input to the model
data <- list()
data$Amax = 20													#Maximum age group
data$SSstart = 1985												#State-space process starts from this year
data$Cdata = as.matrix(catch_data)								#Add catch data
data$Nc = length(catch_data[,1])									#Number of years with catch data
data$pup_production = as.matrix(pup_production)					#Pup production estimates
data$Np = length(pup_production[,1])								#Number of years with pup production data
data$Fecundity = as.matrix(fecundity)							#Observed fecundity rates
data$Nf = length(fecundity$V1)									#Number of years with fecundity data
data$Pmat = as.matrix(Pmatrix)									#Birt ogives
data$Npred = years.of.prediction									#Number of years to run projections
data$priors = as.matrix(priors)									#Priors for estimated parameters
data$Npriors = length(priors$V1)									#Number of priors
data$xmax = data$Cdata[data$Nc,1]+data$Npred+1-data$SSstart+1	#Length of state-space process
data$CQuota = c(0,0)											#Catch level in future projections

###############
#Parameter part
###############

#Initial values
Kinit = 2000000								#Initial population size
Minit = 0.09									#Natural adult mortality
M0init = 0.27								#Natural pup mortality
finit = 0.7									#Mean fecundity rate
ainit = 0.7									#AR(1) parameter
sigmainit = 0.85								#sigma

#Transform some parameters to ensure that estimates are > 0
parameters <- list()
parameters$logK= log(Kinit)  			
parameters$Mtilde= logit(Minit)				
parameters$M0tilde= logit(M0init)				
parameters$ftilde = logit(finit)
parameters$atilde = ainit
parameters$logSigma = log(sigmainit)
parameters$u=rep(0,length.out = (data$Nc+data$Npred-data$SSstart+data$Cdata[1,1]))


obj <- MakeADFun(data,parameters,random="u",DLL="hessm")

obj$fn()
obj$gr()

system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))

rep<-sdreport(obj, getJointPrecision=TRUE)
rep.matrix <- summary(rep)
rep.rnames = rownames(rep.matrix)
indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
indFt = which(rep.rnames=="Ft");indFt <- indFt[-1]
yrs = 1946:2030

#Extract parameters
Kest = exp(opt$par[1])          
Mest = ilogit(opt$par[2])       
M0est = ilogit(opt$par[3])      
fest = ilogit(opt$par[4])       
aest = opt$par[5]             
sigmaest = exp(opt$par[6])    

#allsd<-sqrt(diag(solve(rep$jointPrecision)))
#plsd <- obj$env$parList(par=allsd)
quartz("",9,7)
plot(1946:2030,rep.matrix[indN1,1],type = "l",xlab = "Year",ylab = "Abundance",xlim = c(1946,2030),ylim = c(0,2000000),lwd = 3,col = "darkgreen")
lines(1946:2030,rep.matrix[indN0,1],col = "blue",lwd = 3)
lines(data$pup_production[,1],data$pup_production[,2],type = "p",col = "red",cex = 1.2,pch = 16)

