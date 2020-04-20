library(rSPAMM)

data("harpeastDemo")

data = harpeast$data
parameters = harpeast$parameters

optobj <- run.model(data = data, par = parameters)

res <- model.results(data = data,optobject = optobj)


component = c("N0","Ntot")
