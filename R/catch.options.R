#' Internal helper function for finding equilibrium quota
#'
#' Finding D-based equilibrium quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @param root.dir Root directory for population model files
#' @return minD Minimum found for difference between the depletion coefficient D and 1
#' @keywords population model
#' @export

eq.quota.helper.D <- function(Tot,population,quota, root.dir="C:/Users/a5406/Documents/Populasjonsmodellering/rSPAMM/rSPAMM-master")
{
  # Function called from find.eq.quota, that should be used!
  #setwd(root.dir)
  source('../R/inputFunctions.R')
  source('../R/run_assessment_models.R')
  source('../R/outputFunctions.R')
  source('../R/functions.R')
  data <- load.data(population=population, catch_quota=Tot*quota)
  parameters <- load.initial.values(population)
  obj <- load.model.object(data, parameters, template='harps_and_hoods_population_model2')
  opt <- nlminb(obj$par,obj$fn,obj$gr)
  rep<-sdreport(obj, getJointPrecision=TRUE)
  abs(1-rep$value[match('D1', names(rep$value))])
}


#' Internal helper function for finding equilibrium quota
#'
#' Finding N-based equilibrium quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @param root.dir Root directory for population model files
#' @return minN Minimum found for difference between current and projected 1+ population size
#' @keywords population model
#' @export

eq.quota.helper.N1 <- function(Tot,population,quota, root.dir="C:/Users/a5406/Documents/Populasjonsmodellering/rSPAMM/rSPAMM-master")
{
  # Function called from find.eq.quota, that should be used!
  #setwd(root.dir)
  source('../R/inputFunctions.R')
  source('../R/run_assessment_models.R')
  source('../R/outputFunctions.R')
  source('../R/functions.R')
  data <- load.data(population=population, catch_quota=Tot*quota)
  parameters <- load.initial.values(population)
  obj <- load.model.object(data, parameters, template='harps_and_hoods_population_model2')
  opt <- nlminb(obj$par,obj$fn,obj$gr)
  rep<-sdreport(obj, getJointPrecision=TRUE)
  N1Cur <- rep$value[match("N1CurrentYear", names(rep$value))]
  N1Proj <- rep$value[match("D1New", names(rep$value))]
  abs(N1Cur-N1Proj)
}


#' Main function for finding equilibrium quota
#'
#' Finding D-based or N-based equilibrium quota
#' @param MIN Minimum extent of search window for total catch 
#' @param MAX Maximum extent of search window for total catch 
#' @param quota Proportional catch of 0 and 1+ animals
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param method Set whether D-based (depletion coefficient D) or N-based (1+ population size) criterion should be used for optimisation 
#' @return qEq Optimum quota for achieving equilibrium projected population size
#' @keywords population model
#' @export

find.eq.quota <- function(MIN=1000,MAX=40000,quota=c(0,1),population="harpwest",method = "Dbased")
{
  # Function to find equilibrium quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
    tmp = optimize(eq.quota.helper.D,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  } else
  {
    tmp = optimize(eq.quota.helper.N1,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  }
  cat("Equilibrium quota for",population,"(pups,adults,total):",round(tmp$minimum*quota),sum(round(tmp$minimum*quota)),"\n")
  invisible(tmp$minimum*quota)
}


#' Internal helper function for finding N70 quota
#'
#' Finding N-based N70 quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @param root.dir Root directory for population model files
#' @return minN70 Minimum found for difference between N70 and the projected total population size
#' @keywords population model
#' @export

N70.helper.Nmax <- function(Tot,population,quota, root.dir="C:/Users/a5406/Documents/Populasjonsmodellering/rSPAMM/rSPAMM-master")
{
  setwd(root.dir)
  source('../R/inputFunctions.R')
  source('../R/run_assessment_models.R')
  source('../R/outputFunctions.R')
  source('../R/functions.R')
  data <- load.data(population=population, catch_quota=Tot*quota)
  parameters <- load.initial.values(population)
  obj <- load.model.object(data, parameters, template='harps_and_hoods_population_model2')
  opt <- nlminb(obj$par,obj$fn,obj$gr)
  rep <- sdreport(obj, getJointPrecision=TRUE)

  indNTot <- which(names(rep$value)=="NTot")
  indNTot <- indNTot[-1]
  indCur <- diff(range(data$Cdata[,1]))+1
  NTot <- rep$value[indNTot]
  NTotSD <- rep$sd[indNTot]
  N70 <- 0.7*max(NTot[c(1:indCur)])
  
  ##Npred <- NTot[indCur+10]-qnorm(1-0.1)*NTotSD[indCur+10]
  ## Above line from original code. Changed to qnorm(1-0.2) 
  ## to reflect criterion of 80% rather than 90% probability
  
  Npred <- NTot[indCur+15]-qnorm(1-0.1)*NTotSD[indCur+15]
  
  if(Npred>0) {
    abs(N70-Npred)
  } else {
    99999
  }  
}


#' Internal helper function for finding N70 quota
#'
#' Finding D-based N70 quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @param root.dir Root directory for population model files
#' @return minD07 Minimum found for difference between the depletion coefficient D and 0.7
#' @keywords population model
#' @export

N70.helper.D <- function(Tot,population,quota, root.dir="C:/Users/a5406/Documents/Populasjonsmodellering/rSPAMM/rSPAMM-master")
{
  # Function called from find.N70 that should be used!
  #setwd(root.dir)
  source('../R/inputFunctions.R')
  source('../R/run_assessment_models.R')
  source('../R/outputFunctions.R')
  source('../R/functions.R')
  data <- load.data(population=population, catch_quota=Tot*quota)
  parameters <- load.initial.values(population)
  obj <- load.model.object(data, parameters, template='harps_and_hoods_population_model2')
  opt <- nlminb(obj$par,obj$fn,obj$gr)
  rep <- sdreport(obj, getJointPrecision=TRUE)
  
  indNTot <- which(names(rep$value)=="NTot")
  indNTot <- indNTot[-1]
  indCur <- diff(range(data$Cdata[,1]))+1
  NTot <- rep$value[indNTot]
  NTotSD <- rep$sd[indNTot]
  N70 <- 0.7*max(NTot[c(1:indCur)])
  
  D1 <- rep$value[match('D1', names(rep$value))]
  D1SD <- rep$sd[match('D1', names(rep$value))]
  ## Dpred = D1New-qnorm(1-0.1)*std$Dnew.std
  ## Above line from original code. Changed to qnorm(1-0.2) 
  ## to reflect criterion of 80% rather than 90% probability
  Dpred = D1-qnorm(1-0.1)*D1SD
  
  if(Dpred>0) {
    abs(0.7-Dpred)
  } else {
    99999
  }  
}


#' Main function for finding N70 quota
#'
#' Finding D-based or N-based N70 quota
#' @param MIN Minimum extent of search window for total catch 
#' @param MAX Maximum extent of search window for total catch 
#' @param quota Proportional catch of 0 and 1+ animals
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param method Set whether D-based (depletion coefficient D) or N-based (total population size) criterion should be used for optimisation 
#' @return q70 Optimum quota for achieving projected size of 70% of maximum population size
#' @keywords population model
#' @export

find.N70.quota <- function(MIN=5000,MAX=50000,quota=c(0,1),population="harpwest",method = "Dbased")
{
  # Function to find 70% quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
    tmp = optimize(N70.helper.D,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  }
  if (method == "Nbased") {
    tmp = optimize(N70.helper.Nmax,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  }
  cat("N70 quota for",population,"(pups,adults,total):",round(tmp$minimum*quota),sum(round(tmp$minimum*quota)),"\n")
  
  tmp$minimum*quota
}


#' Function for finding PBR quota
#'
#' Finding Potential Biological Removal (PBR) quotas
#' @param n0 Estimated current population size of 0 animals
#' @param n1 Estimated current population size of 1+ animals
#' @param se0 Standard error of estimate for n0
#' @param se1 Standard error of estimate for n1
#' @param rMax Maximum rate of population increase (by default 0.12, commonly used for pinnipeds)
#' @param Fr Assumed recovery factor
#' @param quota Proportional catch of 0 and 1+ animals
#' @return Returns Minimum projected population size (Nmin) and Potential Biological Removal (PBR),
#' along with the PBR divided by 0 and 1+ animals given specified quota 
#' @keywords population model
#' @export

PBR <- function(n0=partab[4,3], n1=partab[5,3], se0=partab[4,4], se1=partab[5,4],
                rMax=0.12, Fr=0.5, quota=c(0,1), cv=NA) {
  
  if(is.na(cv)) cv <- sqrt((se0^2) + (se1^2)+(2*se0*se1))/(n0+n1)
  
  Nmin <- round((n0+n1) / exp(0.842*sqrt(log(1+cv^2))))
  pbr <- round(0.5 * rMax * Fr * Nmin)
  
  quota <- as.vector(quota)
  
  list(Nmin=Nmin, CV=cv, PBR=pbr, p0=round(pbr*quota[1]), p1=round(pbr*quota[2]))
}
