#' Internal helper function for finding equilibrium quota
#'
#' Finding D-based equilibrium quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @return minD Minimum found for difference between the depletion coefficient D and 1
#' @keywords population model
#' @export

eq.quota.helper.D <- function(Tot,population,quota)
{
  # Function called from find.eq.quota, that should be used!
  dataD <- load.data(population=population, catch_quota=Tot*quota)
  parametersD <- load.initial.values(population)
  objD <- load.model.object(dataD, parametersD)
  optD <- nlminb(objD$par,objD$fn,objD$gr)
  repD<-sdreport(objD, getJointPrecision=TRUE)
  abs(1-repD$value[match('D1', names(repD$value))])
}


#' Internal helper function for finding equilibrium quota
#'
#' Finding N-based equilibrium quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @return minN Minimum found for difference between current and projected 1+ population size
#' @keywords population model
#' @export

eq.quota.helper.N1 <- function(Tot,population,quota)
{
  # Function called from find.eq.quota, that should be used!
  dataN1 <- load.data(population=population, catch_quota=Tot*quota)
  parametersN1 <- load.initial.values(population)
  objN1 <- load.model.object(dataN1, parametersN1)
  optN1 <- nlminb(objN1$par,objN1$fn,objN1$gr)
  repN1<-sdreport(objN1, getJointPrecision=TRUE)
  N1Cur <- repN1$value[match("N1CurrentYear", names(repN1$value))]
  N1Proj <- repN1$value[match("D1New", names(repN1$value))]
  abs(N1Cur-N1Proj)
}


#' Main function for finding equilibrium quota
#'
#' Finding D-based or N-based equilibrium quota
#' @param MIN Minimum extent of search window for total catch 
#' @param MAX Maximum extent of search window for total catch 
#' @param quota Proportional catch of 0 and 1+ animals
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param method Set whether D-based (depletion coefficient D) or N-based (1+ population size) criterion should be used for optimisation  (Dbased/Nbased)
#' @return qEq Optimum quota for achieving equilibrium projected population size
#' @keywords population model
#' @export

find.eq.quota <- function(MIN=1000,MAX=40000,quota=c(0,1),population="harpeast",method = "Dbased")
{
  # Function to find equilibrium quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
    tmp = optimize(eq.quota.helper.D,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  } else
  {
    tmp = optimize(eq.quota.helper.N1,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  }
  cat("-----------------------------------------------------\n\n")
  cat("Equilibrium quota for",population,":\n")
  cat("Pups:  ",round(tmp$minimum*quota)[1],"\n")
  cat("Adults:",round(tmp$minimum*quota)[2],"\n")
  cat("Total :",sum(round(tmp$minimum*quota)),"\n\n")
  cat("-----------------------------------------------------")
  invisible(tmp$minimum*quota)
}


#' Internal helper function for finding N70 quota
#'
#' Finding N-based N70 quota
#' @param Tot Total annual catch for population projections 
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param quota Proportional catch of 0 and 1+ animals
#' @return minN70 Minimum found for difference between N70 and the projected total population size
#' @keywords population model
#' @export

N70.helper.Nmax <- function(Tot,population,quota)
{
  dataNmax <- load.data(population=population, catch_quota=Tot*quota)
  parametersNmax <- load.initial.values(population)
  
  objNmax <- load.model.object(dataNmax, parametersNmax)
  optNmax <- nlminb(objNmax$par,objNmax$fn,objNmax$gr)
  repNmax <- sdreport(objNmax, getJointPrecision=TRUE)

  indNTot <- which(names(repNmax$value)=="NTot")
  indNTot <- indNTot[-1]
  indCur <- diff(range(dataNmax$Cdata[,1]))+1
  NTot <- repNmax$value[indNTot]
  NTotSD <- repNmax$sd[indNTot]
  N70 <- 0.7*max(NTot[c(1:indCur)])
  
  Npred <- NTot[indCur+15]-qnorm(1-0.1)*NTotSD[indCur+15]
  
  if(Npred>0) {
    return(abs(N70-Npred))
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
#' @return minD07 Minimum found for difference between the depletion coefficient D and 0.7
#' @keywords population model
#' @export

N70.helper.D <- function(Tot,population,quota)
{
  # Function called from find.N70 that should be used!
  dataD <- load.data(population=population, catch_quota=Tot*quota)
  parametersD <- load.initial.values(population)
  objD <- load.model.object(dataD, parametersD)
  optD <- nlminb(objD$par,objD$fn,objD$gr)
  repD <- sdreport(objD, getJointPrecision=TRUE)
  
  
  DNmax <- repD$value[match('DNmax', names(repD$value))]
  DNmaxSD <- repD$sd[match('DNmax', names(repD$value))]
  Dpred = DNmax-qnorm(1-0.1)*DNmaxSD
  
  if(Dpred>0) {
    return(abs(0.7-Dpred))
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
#' @param method Set whether D-based (depletion coefficient D) or N-based (total population size) criterion should be used for optimisation (Dbased,Nbased)
#' @return q70 Optimum quota for achieving projected size of 70% of maximum population size
#' @keywords population model
#' @export

find.N70.quota <- function(MIN=5000,MAX=50000,quota=c(0,1),population="harpeast",method = "Nbased")
{
  # Function to find 70% quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
    tmp = optimize(N70.helper.D,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  }
  if (method == "Nbased") {
    tmp = optimize(N70.helper.Nmax,lower=MIN,upper=MAX,population=population,quota=quota,tol=50)
  }
  #cat("N70 quota for",population,"(pups,adults,total):",round(tmp$minimum*quota),sum(round(tmp$minimum*quota)),"\n")
  cat("-----------------------------------------------------\n\n")
  cat("N70 quota for",population,":\n")
  cat("Pups:  ",round(tmp$minimum*quota)[1],"\n")
  cat("Adults:",round(tmp$minimum*quota)[2],"\n")
  cat("Total :",sum(round(tmp$minimum*quota)),"\n\n")
  cat("-----------------------------------------------------")
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

PBR <- function(n0=n0, n1=n1, se0=se0, se1=se1,
                rMax=0.12, Fr=0.5, quota=c(0.15,1-0.15), cv=NA) {
  
  if(is.na(cv)) cv <- sqrt((se0^2) + (se1^2)+(2*se0*se1))/(n0+n1)
  
  Nmin <- round((n0+n1) / exp(0.842*sqrt(log(1+cv^2))))
  pbr <- round(0.5 * rMax * Fr * Nmin)
  
  quota <- as.vector(quota)
  
  list(Nmin=Nmin, CV=cv, PBR=pbr, p0=round(pbr*quota[1]), p1=round(pbr*quota[2]))
}
