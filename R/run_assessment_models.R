#' Load the population model for harp seals and hooded seals
#'
#' Load the model object to be optimized.
#' @param population Choose which population to run the model on (harpwest,harpeast,hooded).
#' @param Amax Maximum age group. Default is 20 years.
#' @param years_of_prediction Number of years in the future to project the model. Default is 15 years.
#' @param Fproj Which fecundity rate to use in future projections. Fproj = "mean" uses mean value of observed fecundity rates. Otherwise a fixed Fproj can be set.
#' @return obj Model object to be optimized.
#' @keywords population model
#' @export
#' @examples
#' load.model.object(population = "harpeast")

load.model.object <- function(dat = data,par = parameters,template='harps_and_hoods_population_model')
{
  require(TMB)
  cat('Compiling model dll, be patient....\n\n')
  flush.console()
  
  compile(paste0("./R/", template, ".cpp"),"-O1 -g",DLLFLAGS="")
 
  cat('\n\nDone!\n')
  flush.console()
  
  cat('\nLoading model dll....\n')
  flush.console()

  dyn.load(dynlib(paste0("./R/", template)))

  cat('\n\nDone!\n')
  flush.console()
  
  cat('\nBuilding model object....\n')
  flush.console()
  
  obj <- MakeADFun(data=dat,parameters=par,DLL=template)
  
  cat('Done!\n')
  flush.console()
  
  obj
}

#' Run the loaded population model for harp seals and hooded seals
#'
#' Optimize the model object.
#' @param obj Model object to be optimized.
#' @return Output from optimized model.
#' @keywords population model
#' @export
#' @examples
#' run.model()


run.model <- function(object=obj)
{
  nlminb(object$par,object$fn,object$gr)
}


#' Get results from optimized population model for harp seals and hooded seals
#'
#' Get model results.
#' @param dat Data object that was used to fit the model
#' @param object The TMB model object
#' @param optimized Output from the optimization
#' @param to.file Logical, whether to return results to the workspace (default) or save to file. NOT YET IMPLEMENTED!
#' @return results Results returned to the workspace or saved to file.
#' @keywords population model
#' @export
#' @examples
#' load.model.object(population = "harpeast")

model.results <- function(dat=data, object=obj, optimized=opt) 
{
  rep<-sdreport(object, getJointPrecision=TRUE)
  rep.matrix <- summary(rep)
  rep.rnames = rownames(rep.matrix)
  indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
  indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
  indNTot = which(rep.rnames=="NTot");indNTot <- indNTot[-1]
  indD1 = which(rep.rnames=="D1");
  indD1New = which(rep.rnames=="D1New");   #NEED TO BE FIXED

##  indN0Current = which(rep.rnames=="N0CurrentYear");
##  indN1Current = which(rep.rnames=="N1CurrentYear");
##  indNTotCurrent = which(rep.rnames=="NTotCurrentYear");
## Above lines replaced, seems TMB has wrong indexes for 
## N0Current and NTotCurrent
  
  yrs = c(min(data$Cdata[,1]):(max(dat$Cdata[,1])+data$Npred+1))
  
  curYr <- match(max(dat$Cdata[,1]), yrs)
  indN0Current <- indN0[curYr]
  indN1Current <- indN1[curYr]
  indNTotCurrent <- indNTot[curYr]
  
  #Extract parameters
  Kest = exp(optimized$par[1])
  Mest = ilogit(optimized$par[2])
  M0est = ilogit(optimized$par[3])
  
  D1 = rep.matrix[indD1,1]
  D1New = rep.matrix[indD1New,1]
  N0Current = rep.matrix[indN0Current,1]
  N1Current = rep.matrix[indN1Current,1]
  NTotCurrent = rep.matrix[indNTotCurrent,1]
  D1.sd = rep.matrix[indD1,2]
  D1New.sd = rep.matrix[indD1New,2]
  N0Current.sd = rep.matrix[indN0Current,2]
  N1Current.sd = rep.matrix[indN1Current,2]
  NTotCurrent.sd = rep.matrix[indNTotCurrent,2]
  
  ## Some suggestions:
  cur.yr <- dim(dat$Cdata)[1]
  N1Current <- rep.matrix[indN1[cur.yr],1]
  N1Current.sd <- rep.matrix[indN1[cur.yr],2]
    
  list(rep=rep, rep.matrix=rep.matrix, rep.rnames=rep.rnames, indN0=indN0,
       indN1=indN1, indNTot=indNTot, indD1=indD1, indD1New=indD1New, indN0Current=indN0Current,
       indN1Current=indN1Current, indNTotCurrent=indNTotCurrent, 
       years=yrs, Kest=Kest, Mest=Mest, M0est=M0est, 
       D1=D1, D1New=D1New, N0Current=N0Current, N1Current=N1Current, 
       NTotCurrent=NTotCurrent, D1.sd=D1.sd, D1New.sd=D1New.sd,
       N0Current.sd=N0Current.sd, N1Current.sd=N1Current.sd, NTotCurrent.sd=NTotCurrent.sd)
}


