source('R/inputFunctions.R')
source('R/run_assessment_models.R')
source('R/outputFunctions.R')
source('R/functions.R')
source('R/catch.options.R')


#Set population
population = 'harpeast'
#Load the data
data <- load.data(population = population)

#Prepare parameters to be estimated
parameters <- load.initial.values(population = population)

#Run the assessment model
opt <- run.model()

#Get the model results
res <- model.results()

#######
# Plot catch data
# Plot fecundity data