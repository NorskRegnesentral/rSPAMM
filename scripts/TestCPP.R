#compile("harp_and_hoods.cpp")
compile("harps_and_hoods_population_model2.cpp","-O1 -g",DLLFLAGS="")
dyn.load(dynlib("harps_and_hoods_population_model2"))
obj <- MakeADFun(data,parameters,DLL="harps_and_hoods_population_model2")


system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))

rep<-sdreport(obj, getJointPrecision=TRUE)
rep.matrix <- summary(rep)
rep.rnames = rownames(rep.matrix)
indN0 = which(rep.rnames=="N0");indN0 <- indN0[-1]
indN1 = which(rep.rnames=="N1");indN1 <- indN1[-1]
indD1 = which(rep.rnames=="D1");
indD1New = which(rep.rnames=="D1New");   #NEED TO BE FIXED
indN0Current = which(rep.rnames=="N0CurrentYear");
indNTot = which(rep.rnames=="NTot")
indNTotmax = which(rep.rnames=="NTotmax")
indDNmax = which(rep.rnames=="DNmax")
indNTotCurrent = which(rep.rnames=="NTotCurrentYear")

D1 = rep.matrix[indD1,1]
D1New = rep.matrix[indD1New,1]
N0Current = rep.matrix[indN0Current,1]
D1.sd = rep.matrix[indD1,2]
D1New.sd = rep.matrix[indD1New,2]
N0Current.sd = rep.matrix[indN0Current,2]

NTot = rep.matrix[indNTot,1]
NTotmax = rep.matrix[indNTotmax,1]
NTotcurrent = rep.matrix[indNTotCurrent,1]
DNmax = rep.matrix[indDNmax,1]

NTot
NTotmax
