#Sample script for a complete analysis of the 
#harps seal population in the East Ice (in the White Sea)
#Using the rSPAMM package

#Load the rSPAMM package
library(rSPAMM)

###################
# Data
###################
# You can choose to use either the demo data which is 
# included in the rSPAMM package or loading the full data

#Use this if analysing the demo data.
#Loading demo data and parameters
data("harpeastDemo")
data = harpeast$data
parameters = harpeast$parameters


#Use this if analysing latest version of the
#complete data set
#Download full data. Default the function will
#ask to create a new folder. If you already have 
#a folder you can add the option 
#"chooseFolder = FALSE".
#Note the Working Directory has to be set to the
#root folder of the downloaded data, i.e., not 
#the Data or the Scripts folder.
downloadData()
#Loading the data
data = load.data(population = "harpeast")
#Loading the parameters
parameters = load.initial.values(population = "harpeast")


##################
# Model fitting
##################
#Run the model
optobj = run.model(data = data, par = parameters)

#Obtain the results
res = model.results(data = data,optobject = optobj)

#Create a nice table with the results
partab = par.table(results=res, dat=data) 

#Plot the results (both the pup model fit and the N1+ population)
plotRes(res,data,grDev = TRUE)


##################
# Catch options
##################

#Find the equilibrium catch level
#In this example we are assuming 15% pups 
#and 85% 1+ animals in the catch
EquilibriumCatch = find.eq.quota(data = data,
                                  parameters = parameters, 
                                  quota = c(0.15,0.85))

#Rerun the model using the estimated equilibrium catch level
data$CQuota = EquilibriumCatch
optEq = run.model(data = data, par = parameters)
resEq = model.results(data = data,optobject = optEq)

#Plot the estimated future trajectory using equilibrium catch level
plotRes(resEq,data,grDev = TRUE)

#--------------------------

#Find the N70 catch level
catchN70 = find.N70.quota(data = data,
                          parameters = parameters, 
                          quota = c(0,1))

#For this population it turns out that the current population
#size is below N70. Because of this there is no point estimating
#the N70 catch level.

#If the current population size was above N70 you could
#rerun the model using the estimated N70 catch level
data$CQuota = catchN70
optN70 = run.model(data = data, par = parameters)
resN70 = model.results(data = data,optobject = optN70)

#Plot the estimated future trajectory using the N70 catch level
plotRes(resN70,data,grDev = TRUE)

#--------------------------

#Find the PBR catch level
#In this example we are assuming 14% pups 
#and 86% 1+ animals in the catch
pbrCatch = PBR(n0=partab[4,3], 
               n1=partab[5,3], 
               se0=partab[4,4], 
               se1=partab[5,4])

#Re-run the model using the PBR catch level
data$CQuota = c(pbrCatch$n0catch,pbrCatch$n1catch)
optPBR = run.model(data = data, par = parameters)
resPBR = model.results(data = data,optobject = optPBR)

#Plot the estimated future trajectory using PBR catch level
plotRes(resPBR,data,grDev = TRUE)
