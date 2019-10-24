library(rSPAMM)

#Set population
population = 'harpeast'
#Load the data
data <- load.data(population = population)

#Prepare parameters to be estimated
parameters <- load.initial.values(population = population)

#Load the model object
obj = load.model.object(dat = data,par = parameters,template='harps_and_hoods_population_model2')

#Run the assessment model
opt <- run.model(object = obj)

#Get the model results
res <- model.results(dat = data,object = obj,optimized = opt)

#Catch quotas
## Find equilibrium quota: 
## 0% pups:
qEq <- find.eq.quota(population='harpeast', method='Dbased')
datEq <- load.data('harpeast', catch_quota=qEq)
objEq <- load.model.object(dat=datEq, template='harps_and_hoods_population_model2')
optEq <- run.model(objEq)
resEq <- model.results(datEq, objEq, optEq)
partEq <- par.table(resEq, datEq) 

## Find N70 quota: 

## 0% pups:
q70 <- find.N70.quota(population='harpeast', method='Dbased')
dat70 <- load.data('harpeast', catch_quota=q70)
obj70 <- load.model.object(dat=dat70, template='harps_and_hoods_population_model2')
opt70 <- run.model(obj70)
res70 <- model.results(dat70, obj70, opt70)
part70 <- par.table(res70, dat70) 

## PBR quotas:

pbr.05 <- PBR(quota=part70[c(4,5),3]/part70[6,3])


#######
# Plot catch data
# Plot fecundity data
# Update manual info for all functions and adding
# examples to each function.