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

load.model.object <- function(data = data,parameters = parameters)
{
  compile("harps_and_hoods_population_model.cpp","-O1 -g",DLLFLAGS="")
  dyn.load(dynlib("harps_and_hoods_population_model"))

  obj <- MakeADFun(data,parameters,DLL="harps_and_hoods_population_model")


}
