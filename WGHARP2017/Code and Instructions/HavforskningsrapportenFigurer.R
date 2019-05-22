#xha2 = run.model("harpeast",c(0,0),T,T,years.of.prediction=12,Pval = 2, Fval = 2)
setwd("~/Work/HarpsHoods/WGHARPModel/WGHARP2013A/")
source("functions.R")
yrs = 1946:2014
i1 = 1
i2 = length(yrs)		#last year for catch




# sigmaha = 21062
# sigmaho = 1665

# xha$Ntot = xha$N0+xha$N1
# xha2$Ntot = xha2$N0+xha2$N1
# xho$Ntot = xho$N0+xho$N1

# #xha$lowerN0 = xha$N0 - 1.96*sigmaha
# #xha$upperN0 = xha$N0 + 1.96*sigmaha

# xha$lowerNtot = xha$N0+xha$lowerN1
# xha$upperNtot = xha$N0+xha$upperN1

# xho$lowerNtot = xho$N0+xho$lowerN1
# xho$upperNtot = xho$N0+xho$upperN1



# Harp West - Norwegian
#xha = run.model("harpwest",c(0,0),T,T,years.of.prediction=15,Pval = 2, Fval = 2)
xha = run.model("harpwest",plotfigs = FALSE)

data <- read.table("~/Work/HarpsHoods/WGHARPModel/WGHARP2013A/harpwest/harpwestcatch.txt",sep = "")
  estha = read.survey.pup.production()
  Nest = dim(estha)
  estha$lower = estha$pup*(1-1.96*estha$cv)
  estha$upper = estha$pup*(1+1.96*estha$cv)
    

quartz("",12,8)
par(mar=c(6,5,4,5),bg = "white")
par(xpd=T)

scalef = 100000
plot(yrs[i1:i2],xha$total[i1:i2]/scalef,type = "n",ylim = c(-6,10))
plot(yrs[i1:i2],xha$total[i1:i2]/scalef,type = "l",col = "orange",lwd = 4,ylim = c(-6,10),yaxt = "n",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l",xlab = "År",ylab="",main = "Grønlandssel i Vesterisen")
axis(2,at=c(0,2,4,6,8,10),labels=c(0,2,4,6,8,10),cex.axis=1.5)

newcatch = array(0,length(data$V1))
for (i in 1:length(newcatch)){
	newcatch[i] = data$V2[i]+data$V3[i]
	newcatch[i] = newcatch[i] / 10000
	lines(data$V1[i]*array(1,10),seq(-6,(-6+newcatch[i]),length.out=10),type = "l",lwd = 5,col = "lightseagreen")
	}

axis(4,at=c(-6,-4,-2,0,2),labels=c(0,2,4,6,8),cex.axis=1.5)
text(1939,5,"Bestandsstørrelse (i 100 000)",srt = 90,cex=1.5)
lines(estha$year,estha$pup/scalef,type = "p",col = "red",ps = 5,lwd = 3,pch=19)

for (i in 1:Nest[1]){
  	 lines(seq(estha$year[i]-0.3,estha$year[i]+0.3,length.out = 10),rep(estha$upper[i]/scalef,10),type = "l",col = "red")
	  
	  lines(seq(estha$year[i]-0.3,estha$year[i]+0.3,length.out = 10),rep(estha$lower[i]/scalef,10),type = "l",col = "red")
	  lines(rep(estha$year[i],10),seq(estha$lower[i]/scalef,estha$upper[i]/scalef,length.out = 10),type = "l",col = "red")
  	}  

lines(yrs[i1:i2],xha$N0[i1:i2]/scalef,lwd = 4,lty = 1,col = "royalblue")
legend(1945,10,legend=c("Totalbestand","Årsunger","Fangstnivå","Merke- og telleestimater av årsunger"),col = c("orange","royalblue","lightseagreen","red"),lwd = c(4,4,4,NA),lty = c(1,1,1,NA),pch = c(NA,NA,NA,19),bty = "n",cex = 1.2)
text(2021,-2,"Fangst (i 10 000)",srt = 90,cex=1.5)




# Hooded - Norwegian
xho = run.model("hooded",plotfigs = FALSE)

data <- read.table("~/Work/HarpsHoods/WGHARPModel/WGHARP2013A/hooded/hoodedcatch.txt",sep = "")
 estho = read.survey.pup.production()
  Nest = dim(estha)
  estho$lower = estho$pup*(1-1.96*estho$cv)
  estho$upper = estho$pup*(1+1.96*estho$cv)
    

quartz("",12,8)
par(mar=c(6,5,4,5),bg = "white")
par(xpd=T)

scalef = 100000
plot(yrs[i1:i2],xho$total[i1:i2]/scalef,type = "l",col = "orange",lwd = 4,ylim = c(-8,14),yaxt = "n",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l",xlab = "År",ylab="",main = "Klappmyss i Vesterisen")
axis(2,at=c(0,2,4,6,8,10,12,14),labels=c(0,2,4,6,8,10,12,14),cex.axis=1.5)

newcatch = array(0,length(data$V1))
for (i in 1:length(newcatch)){
	newcatch[i] = data$V2[i]+data$V3[i]
	newcatch[i] = newcatch[i] / 10000
	lines(data$V1[i]*array(1,10),seq(-8,(-8+newcatch[i]),length.out=10),type = "l",lwd = 5,col = "lightseagreen")
	}

text(1939,7.5,"Bestandsstørrelse (i 100 000)",srt = 90,cex=1.5)
lines(yrs[i1:i2],xho$N0[i1:i2]/scalef,lwd = 4,lty = 1,col = "royalblue")

lines(estho$year,estho$pup/scalef,type = "p",col = "red",ps = 5,lwd = 3,pch=19)

for (i in 1:Nest[1]){
  	 lines(seq(estho$year[i]-0.3,estho$year[i]+0.3,length.out = 10),rep(estho$upper[i]/scalef,10),type = "l",col = "red")
	  
	  lines(seq(estho$year[i]-0.3,estho$year[i]+0.3,length.out = 10),rep(estho$lower[i]/scalef,10),type = "l",col = "red")
	  lines(rep(estho$year[i],10),seq(estho$lower[i]/scalef,estho$upper[i]/scalef,length.out = 10),type = "l",col = "red")
  	}  

axis(4,at=c(-8,-6,-4,-2,0),labels=c(0,2,4,6,8),cex.axis=1.5)
legend(1990,14,legend=c("Totalbestand","Årsunger","Fangstnivå","Telleestimater av årsunger"),col = c("orange","royalblue","lightseagreen","red"),lwd = c(4,4,4,NA),lty = c(1,1,1,NA),pch = c(NA,NA,NA,19),bty = "n",cex = 1.2)
text(2021,-4,"Fangst (i 10 000)",srt = 90,cex=1.5)


# Harp East - Norwegian
xha2 = run.model("harpeast",plotfigs = FALSE)

pupest <- read.table("~/Work/HarpsHoods/WGHARPModel/WGHARP2013A/harpeast/pupestimates.txt",sep = "")
data <- read.table("~/Work/HarpsHoods/WGHARPModel/WGHARP2013A/harpeast/harpeastcatch2.txt",sep = "")
pupest$std <- pupest$V2*pupest$V3

quartz("",12,8)
par(mar=c(6,5,4,5),bg = "white")
scalef = 100000
#plot(yrs[i1:i2],xha2$Ntot[i1:i2]/scalef,type = "n",ylim = c(-8,15))
plot(yrs[i1:i2],xha2$total[i1:i2]/scalef,type = "l",col = "orange",lwd = 4,ylim = c(-6,20),yaxt = "n",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l",xlab = "År",ylab="",main = "Grønlandssel i Østisen")
axis(2,at=c(0,4,8,12,16,20),labels=c(0,4,8,12,16,20),cex.axis=1.5)

newcatch = array(0,length(data$V1))
for (i in 1:length(newcatch)){
	newcatch[i] = data$V2[i]+data$V3[i]
	newcatch[i] = data$V4[i]
	newcatch[i] = newcatch[i] / 20000
	lines(data$V1[i]*array(1,10),seq(-8,(-8+newcatch[i]),length.out=10),type = "l",lwd = 5,col = "lightseagreen")
	}

axis(4,at=c(-8,-6,-4,-2,0,2),labels=c(0,4,8,12,16,20),cex.axis=1.5)

lines(yrs[i1:i2],xha2$N0[i1:i2]/scalef,lwd = 4,lty = 1,col = "royalblue")
lines(pupest$V1,pupest$V2/scalef,type = "p",col = "red",ps = 5,lwd = 3,pch=19)
#legend(1945,10,legend=c("Totalbestand","Årsunger","Fangstnivå","Telleestimater av årsunger"),col = c("orange","royalblue","lightseagreen","red"),lwd = c(4,4,4,NA),lty = c(1,1,1,NA),pch = c(NA,NA,NA,19),bty = "n",cex = 1.2)
legend("topright",legend=c("Totalbestand","Årsunger","Fangstnivå","Telleestimater av årsunger"),col = c("orange","royalblue","lightseagreen","red"),lwd = c(4,4,4,NA),lty = c(1,1,1,NA),pch = c(NA,NA,NA,19),bty = "n",cex = 1.2)

par(xpd=T)
text(2021,-2,"Fangst (i 10 000)",srt = 90,cex=1.5)
text(1939,10,"Bestandsstørrelse (i 100 000)",srt = 90,cex=1.5)



