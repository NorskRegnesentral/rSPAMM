data <- load.data('harpwest')
parameters <- load.initial.values('harpwest')
obj <- load.model.object(template='harps_and_hoods_population_model2')
opt <- run.model()
res <- model.results()
partab <- par.table() 
LCI <- partab$Mean - (1.96*partab$SD)
UCI <- partab$Mean + (1.96*partab$SD)

## Zero catch:
datZ <- load.data('harpwest', catch_quota=c(0,0))
objZ <- load.model.object(dat=datZ, template='harps_and_hoods_population_model2')
optZ <- run.model(objZ)
resZ <- model.results(datZ, objZ, optZ)
partZ <- par.table(resZ, datZ) 


curCatch <- apply(tail(data$Cdata, 5)[,-1], 2, mean)
## Find equilibrium quota: 
## 0% pups:
qEq <- find.eq.quota(method='Dbased')
datEq <- load.data('harpwest', catch_quota=qEq)
objEq <- load.model.object(dat=datEq, template='harps_and_hoods_population_model2')
optEq <- run.model(objEq)
resEq <- model.results(datEq, objEq, optEq)
partEq <- par.table(resEq, datEq) 

## Current pup quota:
qEqM <- find.eq.quota(method='Dbased', quota=curCatch/sum(curCatch))
datEqM <- load.data('harpwest', catch_quota=qEqM)
objEqM <- load.model.object(dat=datEqM, template='harps_and_hoods_population_model2')
optEqM <- run.model(objEqM)
resEqM <- model.results(datEqM, objEqM, optEqM)
partEqM <- par.table(resEqM, datEqM) 


## Find N70 quota: 

## 0% pups:
q70 <- find.N70.quota(method='Dbased')
dat70 <- load.data('harpwest', catch_quota=q70)
obj70 <- load.model.object(dat=dat70, template='harps_and_hoods_population_model2')
opt70 <- run.model(obj70)
res70 <- model.results(dat70, obj70, opt70)
part70 <- par.table(res70, dat70) 

## Current pup quota:
q70M <- find.N70.quota(method='Dbased', quota=curCatch/sum(curCatch))
dat70M <- load.data('harpwest', catch_quota=q70M)
obj70M <- load.model.object(dat=dat70M, template='harps_and_hoods_population_model2')
opt70M <- run.model(obj70M)
res70M <- model.results(dat70M, obj70M, opt70M)
part70M <- par.table(res70M, dat70M) 


## Same with data available in 2016:

dat16 <- load.data('harpwest2016')
parameters16 <- load.initial.values('harpwest2016')
obj16 <- load.model.object(dat=dat16, template='harps_and_hoods_population_model2')
opt16 <- run.model(obj16)
res16 <- model.results(dat16, obj16, opt16)
part16 <- par.table(res16, dat16) 

## Zero catch:
datZ16 <- load.data('harpwest2016', catch_quota=c(0,0))
objZ16 <- load.model.object(dat=datZ16, template='harps_and_hoods_population_model2')
optZ16 <- run.model(objZ16)
resZ16 <- model.results(datZ16, objZ16, optZ16)
partZ16 <- par.table(resZ16, datZ16) 


curCatch <- apply(tail(data$Cdata, 5)[,-1], 2, mean)
## Find equilibrium quota: 
## 0% pups:
qEq16 <- find.eq.quota(population='harpwest2016', method='Dbased')
datEq16 <- load.data('harpwest2016', catch_quota=qEq16)
objEq16 <- load.model.object(dat=datEq16, template='harps_and_hoods_population_model2')
optEq16 <- run.model(objEq16)
resEq16 <- model.results(datEq16, objEq16, optEq16)
partEq16 <- par.table(resEq16, datEq16) 

## Current pup quota:
qEqM16 <- find.eq.quota(population='harpwest2016', method='Dbased', quota=curCatch/sum(curCatch))
datEqM16 <- load.data('harpwest2016', catch_quota=qEqM16)
objEqM16 <- load.model.object(dat=datEqM16, template='harps_and_hoods_population_model2')
optEqM16 <- run.model(objEqM16)
resEqM16 <- model.results(datEqM16, objEqM16, optEqM16)
partEqM16 <- par.table(resEqM16, datEqM16) 


## Find N70 quota: 

## 0% pups:
q7016 <- find.N70.quota(population='harpwest2016', method='Dbased')
dat7016 <- load.data('harpwest2016', catch_quota=q7016)
obj7016 <- load.model.object(dat=dat7016, template='harps_and_hoods_population_model2')
opt7016 <- run.model(obj7016)
res7016 <- model.results(dat7016, obj7016, opt7016)
part7016 <- par.table(res7016, dat7016) 

## Current pup quota:
q70M16 <- find.N70.quota(population='harpwest2016', method='Dbased', quota=curCatch/sum(curCatch))
dat70M16 <- load.data('harpwest2016', catch_quota=q70M16)
obj70M16 <- load.model.object(dat=dat70M16, template='harps_and_hoods_population_model2')
opt70M16 <- run.model(obj70M16)
res70M16 <- model.results(dat70M16, obj70M16, opt70M16)
part70M16 <- par.table(res70M16, dat70M16) 

