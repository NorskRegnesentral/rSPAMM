#' Load and prepare data
#'
#' Load and prepare data to run the assessment model using TMB requirements.
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param Amax Maximum age group. Default is 20 years.
#' @param years_of_prediction Number of years in the future to project the model. Default is 15 years.
#' @param Fproj Which fecundity rate to use in future projections. Fproj = "mean" uses mean value of observed fecundity rates. Otherwise a fixed Fproj can be set.
#' @return data List of loaded data ready for TMB.
#' @keywords input, data
#' @export
#' @examples
#' load.data(population = "harpeast")

load.data <- function(population = "harpeast",Amax = 20,years_of_prediction = 15,Fproj = "mean")
{
  # Read in data ---------------

  # Catch data
  catch_data <- read.table(paste("../Data/",population,"/catch_data.dat",sep = ""),header = FALSE)
  # Pup production estimates
  pup_production <- read.table(paste("../Data/",population,"/pup_production.dat",sep = ""),header = FALSE)
  # Available fecundity data
  fecundity <- read.table(paste("../Data/",population,"/fecundity.dat",sep = ""),header = FALSE)
  # Birth ogives for different periods
  Pdat <- read.table(paste("../Data/",population,"/wgharp.ogi",sep = ""),sep = "",header = TRUE)
  # Which periods the various birth ogives applies to
  Pper <- read.table(paste("../Data/",population,"/wgharp.ogp",sep = ""),header = TRUE)
  #Priors used
  priors <- read.table(paste("../Data/",population,"/priors.dat",sep = ""),header = FALSE)					#Priors used

  #MAYBE ADD STEPWISE CHANGES IN FECUNDITY AND BIRTH OGIVES INSTEAD OF LINEAR TRANSITION
  FecAndP = build.PandF(Fdat = fecundity,cdata = catch_data,Fproj = Fproj,Pdat = Pdat,Pper = Pper,years = c(cdata[1,1],cdata[dim(cdata)[1],1]))
  Fdt = FecAndP$Fdt       #SJEKK VERDIEN PÅ DEN SISTE OG SAMMENLIKN MED WGHARP
  Pmat = FecAndP$Pmatrix

  # Prepare input to the model ------------
  data <- list()
  data$Amax = Amax													#Maximum age group
  data$Cdata = as.matrix(catch_data)								#Add catch data
  data$Nc = length(catch_data[,1])									#Number of years with catch data
  data$pupProductionData = as.matrix(pup_production)					#Pup production estimates
  data$Np = length(pup_production[,1])								#Number of years with pup production data
  data$Ftmp = as.vector(Fdt)                      # Fecundity rates
  data$Pmat = as.matrix(Pmat)									#Birt ogives
  data$Npred = years_of_prediction									#Number of years to run projections
  data$priors = as.matrix(priors)									#Priors for estimated parameters
  data$Npriors = length(priors$V1)									#Number of priors
  #data$CQuota = catch_quota											#Catch level in future projections

  return(data)
}


#' Construct time varying fecundity rates and birth ogives
#'
#' Prepare the fecundity vector and the birth ogive matrix used by the model. A linear
#' transition between years with missing data is used.
#' @param Fdat Observed fecundity rates.
#' @param Fproj Which fecundity rate used in projections. Default is mean value over all observed fecundity rates ("mean"). If not a specified fecundity rate can be used (just insert number).
#' @param Pdat Observed birth ogives.
#' @param Pper For which time periods the birth ogives are valid.
#' @param years Start and stop year for the model (without projections).
#' @return Pmatrix Time varying birth ogives.
#' @keywords Fecundity, birth ogive
#' @export
#' @examples
#' FvecandPmat <- build.PandF(fecundity,catch_data,Fproj,Pdat,Pper)

#MAYBE ADD STEPWISE CHANGES IN FECUNDITY AND BIRTH OGIVES INSTEAD OF LINEAR TRANSITION
build.PandF <- function(Fdat = fecundity,Fproj = Fproj,Pdat = Pdat,Pper = Pper,years = years)
{

  yr1 = years[1]
  yr2 = years[2]
  nyr = yr2 - yr1 + 1

  Fvec = array(0,nyr)
  Fvec[1:(Fdat$V1[1] - yr1 + 1)] = Fdat$V2[1]

  for (i in 2:length(Fdat$V1)){
    if(Fdat$V1[i]-Fdat$V1[i-1]==1){
      Fvec[Fdat$V1[i]-yr1+1] = Fdat$V2[i]
    } else {
      i1 = Fdat$V1[i-1] - yr1 + 1
      i2 = Fdat$V1[i] - yr1 + 1
      Fvec[i1:i2] = seq(Fdat$V2[i-1],Fdat$V2[i],length.out = (i2-i1+1))
    }
  }
  Fvec[i2:nyr] = Fdat$V2[i]

  if(class(Fproj) == "character") Fvec[(yr2-yr1+1):nyr] = mean(Fdat$V2)

  if(class(Fproj)=="numeric") Fvec[(yr2-yr1+1):nyr] = Fproj


  P = matrix(0,nyr,max(Pdat$Age))

  for(i in 1:(Pper$Pstop[1]-yr1+1)){
    P[i,] = Pdat[,2]
  }

  for(i in 2:length(Pper$Pstart)){
    i1 = Pper$Pstop[i-1]-yr1+1
    i2 = Pper$Pstart[i]-yr1+1
    i3 = Pper$Pstop[i]-yr1+1

    for(j in i2:i3){
      P[j,] = Pdat[,i+1]
    }

    for(age in 1:max(Pdat$Age)){
      P[i1:i2,age] = seq(P[i1,age],P[i2,age],length.out = (i2-i1+1))
    }

  }

  for(i in i3:nyr){
    P[i,] = Pdat[,length(Pper$Pstart)+1]
  }

  return(list(Fdt = Fvec,Pmatrix = P))
}

#' Load initial values
#'
#' Load initial values for which to start the optimization from. These are the values to be
#' estimated from the model.
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param fromFile Load initial values from file for specified population or set initial values manually.
#' @param Kinit Initial value for population size.
#' @param Minit Initial value for mortality of 1+ population, i.e., seals of age 1 or more.
#' @param M0init Initial value for pup mortality.
#' @return parameters List of parameters with initial values to be estimated fromt the model.
#' @keywords input, parameters, initial values
#' @export
#' @examples
#' load.initial.values(population = "harpeast")

load.initial.values <- function(population = "harpeast",fromFile = TRUE,Kinit = 2000000,Minit=0.09,M0init=0.27)
{

  if(fromFile == TRUE){
    initial_values <- read.table(paste("../Data/",population,"/initial_values.dat",sep = ""),header = FALSE)

    #Initial values
    Kinit = initial_values[1]								#Initial population size
    Minit = initial_values[2]									#Natural adult mortality
    M0init = initial_values[3]								#Natural pup mortality
  }


  #Transform some parameters to ensure that estimates are > 0
  parameters <- list()
  parameters$logK= log(Kinit) #THESE SHOULD BE BOUNDED PARAMETERS, FIX THAT LATER
  parameters$Mtilde= logit(Minit)
  parameters$M0tilde= logit(M0init)

  return(parameters)
  }
