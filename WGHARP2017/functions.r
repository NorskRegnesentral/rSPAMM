################################################################
# Abundance estimation of harp and hodded seals in White Sea and Greenland Sea
# ICES/WGHARP
#
# Functions to manage input and output files used by wgharp.exe and run wgharp.exe.
# wgharp.exe is the implementation of the model described in ICES CM 2004/ACFM::06
# wgharp.tpl written by Hans Julius Skaug (skaug@imr.no)
#
# Date: 26.08.03
# Modified 31.05.17
# Author: Gjermund Beothun (gjermund@imr.no)
# Edited: Tor Arne Øigård (oigard@nr.no)
################################################################

################################################################
#Copyright (C) 2004  Gjermund Boethun
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#Or visit: http://www.gnu.org/copyleft/gpl.html
################################################################


run.model <- function(directory="harpeast",QUOTA=c(0,0),plot.pred=TRUE,conf.bound=TRUE,years.of.prediction=15,legplace = "topleft",plotfigs = TRUE,mcmc.check = FALSE,nsims=1e5,nthin=10,shiftymax = 0,colour = TRUE,cilevel = 95,Fproj = "mean")
{
  # Main function to run the model.
  # quota:     catch quota to be used in prediction of population size.
  #            Vector c(yearly catch of pups, yearly chatch of adults), default c(0,0)
  # directory: Path to directory containing input files for the model.
  # conf.hist: logical variable for plotting confidence bound for historical data or not.
  # conf.pred: logical variable for plotting confidence bound for predicted data or not.
  # years.of.prediction: Number of years to make prognosis for.
  # Fproj - Set Fecundity rate for projections, NA (default) use last observed Fecundity, "mean" use the average value of all historical Fecundity rates, or set it manually.

  #Catch year t provides estimate for year t+1 (not implemented in tpl file)
  years.of.prediction = years.of.prediction + 1
  copy.input.files(directory)
  buildPandF(Fproj)
  write.quota(QUOTA,years.of.prediction)
  run.wgharp(runmcmc)
  REP=read.abundance("./")
  std=read.standard.deviation("./","wgharp.std",years.of.prediction)
  #Making text and graphical output to screen.
  x=Plot.abundance(REP,std,years.of.prediction,plot.pred,conf.bound,directory,legplace,plotfigs,shiftymax,colour,cilevel)
  N70 = 0.70*max(x$N1[1:(length(x$N0)-years.of.prediction+1)]+x$N0[1:(length(x$N0)-years.of.prediction+1)])
  if(mcmc.check == TRUE){
  	mcmc.eval(nsims = nsims,nthin = nthin,N70 = N70,pos.year = scan("wgharp.cat",n=1,com="#",quiet=T)+11)
  }
  invisible(x)

}


################### Functions for plotting #################
Plot.abundance<-function(rep,std,years.of.prediction,plot.pred=T,conf.bound=F,directory,legplace = "topleft",plotfigs,shiftymax,colour = colour,cilevel = 80)
{
  # rep: abundance
  # std: standard deviation
  # years.of.prediction : number of yers to predict
  # conf.bound.hist: plot 95% confidence bound for historical data
  # conf.bond.pred: plot 95% confidence bound for predicted data
  ciquantile = qnorm(1-(100-cilevel)/(200))
  x=NULL
  x$year = rep$V1
  x$lowerN1 = rep$V3 - ciquantile*std$N1.std    #1.64 eqv 90% conf.bound
  x$upperN1 = rep$V3 + ciquantile*std$N1.std    #1.96 eqv 95% conf.bound qnorm(1-0.025)=1.96
  x$lowerN0 = rep$V2 - ciquantile*std$N0.std    #1.64 eqv 90% conf.bound
  x$upperN0 = rep$V2 + ciquantile*std$N0.std    #1.96 eqv 95% conf.bound qnorm(1-0.025)=1.96
  x$N1 = rep$V3
  x$N0 = rep$V2
  x$total <- x$N1 + x$N0
  x$totalSD <- sqrt(std$N0.std^2+std$N1.std^2)
  x$lowertotal = x$total - ciquantile*x$totalSD    #1.64 eqv 90% conf.bound
  x$uppertotal = x$total + ciquantile*x$totalSD    #1.96 eqv 95% conf.bound qnorm(1-0.025)=1.96
  #x$f = rep$V4

  N70 = 0.70*max(x$N1[1:(length(x$N0)-years.of.prediction+1)]+x$N0[1:(length(x$N0)-years.of.prediction+1)])


  # Text output
  pos.today=scan("wgharp.cat",n=1,com="#",quiet=T)  #in 2003: x$year[pos.today]=2003 (pos.today=58)
  pos.today = pos.today + 1
  #pos.today = 62 #2007 estimate
  cat("\n ###############",directory,"##################")
  cat("\nD1+=",std$D," std(D1+)=",std$D.std,", 95% CI: (",std$D-1.96*std$D.std," - ",std$D+1.96*std$D.std,")\n")
  cat("Dnmax=",std$Dnew," std(Dnmax)=",std$Dnew.std,", 80% CI: (",std$Dnew-0.84*std$Dnew.std," - ",std$Dnew+0.84*std$Dnew.std,")\n")

  cat("95% CI N0 in",x$year[pos.today],": (",std$N0.2003-1.96*std$N0.2003.std,"-",std$N0.2003+1.96*std$N0.2003.std,")\n")
  cat("N70 = ",N70,"\n")
  print(as.data.frame(x)[(x$year==x$year[(pos.today-1)]|x$year==x$year[(pos.today)]),c("year","lowerN1","N1","upperN1","lowerN0","N0","upperN0","lowertotal","total","uppertotal","totalSD")])
  cat("\n\n")
  # D=N_(today+10)/N_today , would be better: D=N_(today+15)/N_(today+5)


  # Graphical output
  # Historical data
  if (plotfigs){
  X11("",width = 9,height = 7)
  par(mfrow=c(3,1))

  # (1,1) Text
  plot(1:10,1:10,type="n",axes=F,xlab="",ylab="",main=directory)
  mtext("Estimates (Prior)",3)
  tmp = scan("wgharp.par",com="#",quiet=T)
  pri = scan("wgharp.pri",com="#",quiet=T)

  text(3,9,paste("K =",round(tmp[1]),"(",pri[1],",",pri[2],")"),pos=4)
  text(3,8,paste("M =",round(tmp[2],3),"(",pri[3],",",pri[4],")"),pos=4)
  text(3,7,paste("M0=",round(tmp[3],3),"(",pri[5],",",pri[6],")"),pos=4)
  #text(3,6,paste("f =",round(tmp[4],4),"(",pri[7],",",pri[8],")"),pos=4)

  tmpN1 = round(as.data.frame(x)[pos.today,c("lowerN1","upperN1","N1")],-1)
  tmpN0 = round(as.data.frame(x)[pos.today,c("lowerN0","upperN0","N0")],-1)
  tmptotal = round(as.data.frame(x)[pos.today,c("lowertotal","uppertotal","total")],-1)

  stringN1 = paste("N1 (",tmpN1[1,1],tmpN1[1,2],"),",tmpN1[1,3],", std = ",std$N1.std[pos.today])
  stringN0 = paste("N0 (",tmpN0[1],tmpN0[2],"),",tmpN0[3],", std = ",std$N0.2003.std)
  stringNtot = paste("Total (",tmptotal[1],tmptotal[2],"),",tmptotal[3],", std = ",round(x$totalSD[pos.today]))

  text(3,4,paste(x$year[pos.today],"95% ci, point est."),pos=4,cex=0.8)
  text(3,3,stringN1,pos=4,cex=0.8)
  text(3,2,stringN0,pos=4,cex=0.8)
  text(3,1,stringNtot,pos=4,cex=0.8)

  # (1,2) Adult
  scala = 1000000 #100 000
  x.scale = 1:pos.today #which(x$year <= 2003)
  #Setting  axes
  if (plot.pred){ Xlim=range(x$year); Ymax=max(x$N1); Ymin=min(x$N1) } else          { Xlim=range(x$year[x.scale]); Ymax=max(x$N1[x.scale]); Ymin=min(x$N1[x.scale]) }

  if (conf.bound) Ylim=c(0,ifelse(plot.pred,x$upperN1[pos.today+years.of.prediction],x$upperN1[pos.today])) else Ylim=c(Ymin,Ymax)

 if (conf.bound) Ylim=c(0,max(x$upperN1)) else Ylim = c(Ymin,Ymax)

  #Trajector
  plot(x$year[x.scale],x$N1[x.scale]/scala,xlim=Xlim,ylim=Ylim/scala,type="l",xlab="Year",ylab="Abundance [10^6]",main="Adults")
  lines(x$year[-x.scale],x$N1[-x.scale]/scala,lty=2) #prediction N1

  # 95% conf lines
  if (conf.bound){
      #tmp=rep(0,length(x$lowerN1[x.scale]))
      #tmp2=rep(0,length(x$lowerN1[-x.scale]))
      lines(x$year[x.scale],x$upperN1[x.scale]/scala,col=3)         #Upper 95%
      lines(x$year[x.scale],x$lowerN1[x.scale]/scala,col=2)         #Lower 95%
      lines(x$year[-x.scale],x$upperN1[-x.scale]/scala,col=3,lty=2) #Upper, perdiction
      lines(x$year[-x.scale],x$lowerN1[-x.scale]/scala,col=2,lty=2) #Lower, perdiction
  }
  if (exists("Yeld")) abline(h=Yeld/scala,lty=3)
  if (exists("Year")) abline(v=Year,lty=3)



  # (2,1) Model fit
  # (2,2) f-plot
  est = read.survey.pup.production()
  tmp=list()
  tmp$x=x
  Plot.fit(tmp,est,F,T,T)


########################
  ## Figure 2

 if(colour == TRUE){
  	coltot = rgb(32,113,178,maxColorValue = 255)
  	coltotci = "aliceblue"
  	colpup = "lightseagreen"
  	colpupci = rgb(217,248,246,maxColorValue = 255)
  	colpupest = "red"
  } else {
  	coltot = "grey20"
  	coltotci = "grey90"
  	colpup = "grey60"
  	colpupci = "grey90"
  	colpupest = "black"
  }


  X11("",10,7)
  par(mar=c(6,5,4,5),bg = "white")
  plot(x$year[x.scale],x$N0[x.scale],xlim=Xlim,ylim=c(0,(Ylim[2]/4)/scala),type="n",col = "darkmagenta",lwd = 4,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l",xlab="Year",ylab="Population size (in 1 000 000)",main="Modelled pup abundance")

if (conf.bound){
  	if(plot.pred){
  	  polygon(c(x$year,rev(x$year)),c(x$lowerN0/scala,rev(x$upperN0/scala)),col = colpupci,border = NA)
      } else {
     	  polygon(c(x$year[x.scale],rev(x$year[x.scale])),c(x$lowerN0[x.scale]/scala,rev(x$upperN0[x.scale]/scala)),col = colpupci,border = NA)
      }
}


  lines(x$year[x.scale],x$N0[x.scale]/scala,lty=1,col = colpup,lwd = 4)


   if(plot.pred){
     lines(x$year[-x.scale],x$N0[-x.scale]/scala,lty=2,col = colpup,lwd = 4)
   }

  est = read.survey.pup.production()
  lines(est$year,est$pup/scala,type = "p",col = colpupest,pch=16)

  Nest = dim(est)
  est$lower = est$pup*(1-1.96*est$cv)
  est$upper = est$pup*(1+1.96*est$cv)

for (i in 1:Nest[1]){
  	 lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$upper[i]/scala,10),type = "l",col = colpupest)

	  lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$lower[i]/scala,10),type = "l",col = colpupest)
	  lines(rep(est$year[i],10),seq(est$lower[i]/scala,est$upper[i]/scala,length.out = 10),type = "l",col = colpupest)
  	}


########################
## Figure 3

  X11("",10,7)
  par(mar=c(6,5,4,5),bg = "white")
  plot(x$year[x.scale],x$N1[x.scale]/scala,xlim=Xlim,ylim=c(0,Ylim[2]/scala),type="n",col = "darkmagenta",lwd = 4,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l",xlab="Year",ylab="Population size (in 1 000 000)",main="Modelled 1+ population abundance")

  if (conf.bound){
  	if(plot.pred){
  	  polygon(c(x$year,rev(x$year)),c(x$lowerN1/scala,rev(x$upperN1/scala)),col = coltotci,border = NA)   #rgb(200,225,245,maxColorValue = 255)

  	  polygon(c(x$year,rev(x$year)),c(x$lowerN0/scala,rev(x$upperN0/scala)),col = colpupci,border = NA)
      } else {
      polygon(c(x$year[x.scale],rev(x$year[x.scale])),c(x$lowerN1[x.scale]/scala,rev(x$upperN1[x.scale]/scala)),col = coltotci,border = NA)
	  polygon(c(x$year[x.scale],rev(x$year[x.scale])),c(x$lowerN0[x.scale]/scala,rev(x$upperN0[x.scale]/scala)),col = colpupci,border = NA)
      }


  }

  lines(x$year[x.scale],x$N1[x.scale]/scala,xlim=Xlim,ylim=c(0,Ylim[2]/scala),col = coltot,lwd = 4)
  lines(x$year[x.scale],x$N0[x.scale]/scala,lty=1,col = colpup,lwd = 4)


   if(plot.pred){
   lines(x$year[-x.scale],x$N1[-x.scale]/scala,lty=2,col = coltot,lwd = 4) #prediction N1

   lines(x$year[-x.scale],x$N0[-x.scale]/scala,lty=2,col = colpup,lwd = 4)
   }



  est = read.survey.pup.production()
  lines(est$year,est$pup/scala,type = "p",col = colpupest,pch=16)
  Nest = dim(est)
  est$lower = est$pup*(1-1.96*est$cv)
  est$upper = est$pup*(1+1.96*est$cv)

  for (i in 1:Nest[1]){
  	 lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$upper[i]/scala,10),type = "l",col = colpupest)

	  lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$lower[i]/scala,10),type = "l",col = colpupest)
	  lines(rep(est$year[i],10),seq(est$lower[i]/scala,est$upper[i]/scala,length.out = 10),type = "l",col = colpupest)
  	}
  legend(legplace,legend=c("N1+ population","N0 population","Pup estimates"),col = c(coltot,colpup,colpupest),lwd = c(4,4,NA),lty = c(1,1,NA),pch=c(NA,NA,16),bty = "n",cex = 1.5)



  #Figure 4
  scala = 1000000 #100 000
  x.scale = 1:pos.today #which(x$year <= 2003)
  #Setting  axes
  if (plot.pred){ Xlim=range(x$year); Ymax=max(x$N1+x$N0); Ymin=min(x$N1) } else          { Xlim=range(x$year[x.scale]); Ymax=max(x$N1[x.scale]+x$N0[x.scale]); Ymin=min(x$N1[x.scale]) }

  if (conf.bound) Ylim=c(0,ifelse(plot.pred,x$upperN1[pos.today+years.of.prediction],x$upperN1[pos.today])) else Ylim=c(Ymin,Ymax)

 #if (conf.bound) Ylim=c(0,max(x$upperN1)) else Ylim = c(Ymin,Ymax)

  if (plot.pred){
  	if (conf.bound) Ylim = c(0,max(x$uppertotal)) else Ylim = c(0,max(x$total))
  } else
  {
  	if (conf.bound) Ylim = c(0,max(x$uppertotal[x.scale])) else Ylim = c(0,max(x$total[x.scale]))
  }



  X11("",10,7)
  par(mar=c(6,5,4,5),bg = "white")
  par(xpd=T)

  plot(x$year[x.scale],(x$total[x.scale])/scala,xlim=Xlim,ylim=c(0,(Ylim[2]+shiftymax)/scala),type="n",col = coltot,lwd = 4,cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l",xlab="Year",ylab="Population size (in 1 000 000)",main="Modelled total population abundance")

  if (conf.bound){
  	if(plot.pred){
  	  polygon(c(x$year,rev(x$year)),c(x$lowertotal/scala,rev(x$uppertotal/scala)),col = coltotci,border = NA)

  	  polygon(c(x$year,rev(x$year)),c(x$lowerN0/scala,rev(x$upperN0/scala)),col = colpupci,border = NA)
      } else {
      polygon(c(x$year[x.scale],rev(x$year[x.scale])),c(x$lowertotal[x.scale]/scala,rev(x$uppertotal[x.scale]/scala)),col = coltotci,border = NA)
	  polygon(c(x$year[x.scale],rev(x$year[x.scale])),c(x$lowerN0[x.scale]/scala,rev(x$upperN0[x.scale]/scala)),col = colpupci,border = NA)
      }


  }
   lines(x$year[x.scale],(x$total[x.scale])/scala,col = coltot,lwd = 4)
   lines(x$year[x.scale],x$N0[x.scale]/scala,lty=1,col = colpup,lwd = 4)

   #prediction N1
   if(plot.pred){
   lines(x$year[-x.scale],x$total[-x.scale]/scala,lty=2,col = coltot,lwd = 4)
   lines(x$year[-x.scale],x$N0[-x.scale]/scala,lty=2,col = colpup,lwd = 4)

   N70 = 0.70*max(x$N1[x.scale]+x$N0[x.scale])
   N50 = 0.50*max(x$N1[x.scale]+x$N0[x.scale])
   N30 = 0.30*max(x$N1[x.scale]+x$N0[x.scale])

   lines(x$year,array(N70/scala,length(x$year)),lty = 1,lwd = 1,col = "grey80")
   lines(x$year,array(N50/scala,length(x$year)),lty = 1,lwd = 1,col = "grey80")
   lines(x$year,array(N30/scala,length(x$year)),lty = 1,lwd = 1,col = "grey80")
   text(x$year[length(x$year)]+5,N70/scala,expression(N[70]))
   text(x$year[length(x$year)]+5,N50/scala,expression(N[50]))
   text(x$year[length(x$year)]+5,N30/scala,expression(N[lim]))
   } else {
   N70 = 0.70*max(x$N1[x.scale]+x$N0[x.scale])
   N50 = 0.50*max(x$N1[x.scale]+x$N0[x.scale])
   N30 = 0.30*max(x$N1[x.scale]+x$N0[x.scale])

   lines(x$year[x.scale],array(N70/scala,length(x.scale)),lty = 1,lwd = 1,col = "grey80")
   lines(x$year[x.scale],array(N50/scala,length(x.scale)),lty = 1,lwd = 1,col = "grey80")
   lines(x$year[x.scale],array(N30/scala,length(x.scale)),lty = 1,lwd = 1,col = "grey80")
   text(x$year[length(x.scale)]+5,N70/scala,expression(N[70]))
   text(x$year[length(x.scale)]+5,N50/scala,expression(N[50]))
   text(x$year[length(x.scale)]+5,N30/scala,expression(N[lim]))

   }

# if (conf.bound){
	# if(plot.pred){
      # lines(x$year,x$uppertotal/scala,col=coltot,lty = 3,lwd = 2) #Upper, perdiction
      # lines(x$year,x$lowertotal/scala,col=coltot,lty = 3,lwd = 2) #Lower, perdiction
      # lines(x$year,x$upperN0/scala,col= colpup,lty = 3,lwd = 2)         #Upper 95%
      # lines(x$year,x$lowerN0/scala,col= colpup,lty = 3,lwd = 2)         #Lower 95%
      # } else {
      # lines(x$year[x.scale],x$uppertotal[x.scale]/scala,col= coltot,lty = 3,lwd = 2)         #Upper 95%
      # lines(x$year[x.scale],x$lowertotal[x.scale]/scala,col= coltot,lty = 3,lwd = 2)         #Lower 95%
      # lines(x$year[x.scale],x$upperN0[x.scale]/scala,col= colpup,lty = 3,lwd = 2)         #Upper 95%
      # lines(x$year[x.scale],x$lowerN0[x.scale]/scala,col= colpup,lty = 3,lwd = 2)         #Lower 95%
      # }
  # }

  est = read.survey.pup.production()
  lines(est$year,est$pup/scala,type = "p",col = colpupest,pch=16)
  Nest = dim(est)
  est$lower = est$pup*(1-1.96*est$cv)
  est$upper = est$pup*(1+1.96*est$cv)

  for (i in 1:Nest[1]){
  	 lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$upper[i]/scala,10),type = "l",col = colpupest)

	  lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$lower[i]/scala,10),type = "l",col = colpupest)
	  lines(rep(est$year[i],10),seq(est$lower[i]/scala,est$upper[i]/scala,length.out = 10),type = "l",col = colpupest)
  	}
 legend(legplace,legend=c("Total population",expression(paste(N[0]," population")),"Pup estimates"),col = c(coltot,colpup,colpupest),lwd = c(4,4,NA),lty = c(1,1,NA),pch=c(NA,NA,16),bty = "n",cex = 1.5)


 ################################


 }
  invisible(as.data.frame(x))
}


Plot.fit <- function(x,est,f.plot=F,fit.plot=T,color=F)
{
# Plots survey estimates of pup production and pup production from model
# Can overlay plots for diferent parameter settings in one plot
# x objects from run.model
# x=NULL
# x$m0.10=run.model(...Different...)
# x$m0.11=run.model(...parameter...)
# x$m0.12=run.model(...settings....)
#
# est = read.survey.pup.production()

  #if (f.plot & fit.plot) par(mfrow=c(2,1)) else par(mfrow=c(1,1))
  ColName = sort(names(x))
  lty = c(1,3,6,2,4,5) #Max 6 senarios in one plot!
  #lty = c(1,18,8424)
  if (color) col = (1:length(ColName)) else col = rep(1,length(ColName))

if (f.plot) {
  # mean birth rate among all females
  # Finding min/max value among all datasets
  tmp = 1:scan("wgharp.cat",n=1,com="#",quiet=T)
  all.y = NULL
  for (i in 1:length(ColName)) {
    all.y = c(all.y,x[[ColName[i]]]$f[tmp])
  }
  #browser()
  plot(x[[ColName[1]]]$year[tmp],x[[ColName[1]]]$f[tmp],
       ylim=range(all.y),
       type="l",xlab="Year",ylab="f")
  # Lines for additional models
  if (length(ColName) > 1){
    for ( i in 2:length(ColName)) {
      lines(x[[ColName[i]]]$year[tmp],x[[ColName[i]]]$f[tmp],
         lty=lty[i],col=col[i])
    } #End for
    cat("Click in plot to place legend\n"); flush.console()
    legend(locator(1),legend=ColName,lty=(lty[1:length(ColName)]),col=col)
  } #End if
  title("Mean birth rate among all females")
}

  # Pup production model vs survey (model fit)
  if (fit.plot) {
    # Calculating 95% confidence bound
    est$lower = est$pup*(1-1.96*est$cv)
    est$upper = est$pup*(1+1.96*est$cv)

    # setting up x scale for plotting
    Xmax=0;Xmin=Inf
    for ( i in ColName){
      Xmax=max( c(Xmax,x[[i]]$N0[x[[i]]$year %in% est$year]) )
      Xmin=min( c(Xmin,x[[i]]$N0[x[[i]]$year %in% est$year]))
    }
    x.scale = (1:length(x[[ColName[1]]]$year))[x[[ColName[1]]]$year %in% est$year]
    x.scale = c(x.scale[1]-1,x.scale,x.scale[length(x.scale)]+2)

    plot(x[[ColName[1]]]$year[x.scale],x[[ColName[1]]]$N0[x.scale],
       	 ylim=c(min(c(Xmin,est$lower)),max(c(Xmax,est$upper))),
 	       type="n",xlab="Year",ylab="Abundance",axes=F)
    axis(1,at=( (min(est$year)-1) : (max(est$year)+2) ) )
    axis(2)
    box()
    dummy=1
    for (i in ColName) {
      lines(x[[i]]$year[x.scale[1]:x.scale[length(x.scale)]],
            x[[i]]$N0[x.scale[1]:x.scale[length(x.scale)]],
            lty=lty[dummy],col=col[dummy])
      dummy=dummy+1
    }
    #Vertical lines (survey estimates with 95% C.I.
    #test1 = matrix(0,1,length(est$year))
    #test2 = test1

    for ( i in 1:length(est$year)) {
      lines(rep(est$year[i],2),c(est$lower[i],est$upper[i]),lwd=2)
      points(est$year[i],est$pup[i],pch=1)
      #test1[i] = est$year[i]
      #test2[i] = est$pup[i]*est$cv[i]
    }
    #library(gplots)
    #plotCI(x=test1,uiw=test2,col="black", barcol="blue", lwd=1,main = "",add=T, ylab = "",xlab="")

  }
  #cat("Click in plot to place legend\n"); flush.console()
  #legend(locator(1),ColName,lty=lty[1:length(ColName)],col=col)
  title("Pups")
  mtext("Model vs Survey",3)
}

################### Misc functions #########################
run.wgharp <-function(runmcmc)
{
  # start the binary code compiled with AD-model builder, see wgharp.tpl for full code.
  #Run the model
		switch(Sys.info()[['sysname']],
   			 Windows= {system("./wgharp.exe")},
    		 Darwin = {system("./wgharp")},Linux ={system("./wgharp")})
}


copy.input.files<-function(directory="harpwest")
{
  to = list.files(directory)
  from = paste(directory,"//",to,sep="")
  Check = file.copy(from,to,overwrite=T)
  if (any(Check==FALSE)) warnings("Could not copy all files!")
  #system(paste("copy_wgharp.bat",from,to))
}

eq.quota.helper.original <- function(Tot,directory,quota)
{
  # Function called from find.eq.quota, that should be used!
  x=run.model(directory,Tot*quota,plot.pred=TRUE,conf.bound=FALSE,years.of.prediction=50,plotfigs = FALSE)
  out=x$N1[x$year==2043]-x$N1[x$year==2015]
  #ev gjoer regresjon pÂ 2010:2040 og gi ut stigningstall?
  #cat("#",Tot*quota,"\n")
  return(abs(out))
}

eq.quota.helper.N1 <- function(Tot,directory,quota)
{
  # Function called from find.eq.quota, that should be used!
    years.of.prediction = 50
  x=run.model(directory,Tot*quota,plot.pred=TRUE,conf.bound=FALSE,years.of.prediction=years.of.prediction,plotfigs = FALSE)
  ind_today = which(x$year == max(x$year-years.of.prediction+1))

  out=x$N1[ind_today]-x$N1[ind_today+30]
  #ev gjoer regresjon pÂ 2010:2040 og gi ut stigningstall?
  #cat("#",Tot*quota,"\n")
  return(abs(out))
}

eq.quota.helper.D <- function(Tot,directory,quota)
{
  # Function called from find.eq.quota, that should be used!
  years.of.prediction = 50
  x=run.model(directory,Tot*quota,plot.pred=TRUE,conf.bound=FALSE,years.of.prediction=years.of.prediction,plotfigs = FALSE)
  std=read.standard.deviation("./","wgharp.std",years.of.prediction)
  out =(1-std$D)
  return(abs(out))
}

find.eq.quota <- function(MIN=180000,MAX=230000,quota=c(0.93,0.07),directory="harpeast",method = "Dbased")
{
  # Function to find equabrilum quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
  tmp = optimize(eq.quota.helper.D,lower=MIN,upper=MAX,directory=directory,quota=quota,tol=50)
  } else
  {
  	tmp = optimize(eq.quota.helper.N1,lower=MIN,upper=MAX,directory=directory,quota=quota,tol=50)
  }
  cat("Equabriluim quota for",directory,"(pups,adults,total):",round(tmp$minimum*quota),sum(round(tmp$minimum*quota)),"\n")
  invisible(tmp$minimum*quota)
}

#This one will can cause equlibrium N70 quota to be smaller than eq. quota because of wide confidence intervals in model predictions
N70.helper.Nmax <- function(Tot,directory,quota)
{
  # Function called from find.N70 that should be used!
  years.of.prediction = 50
  x=run.model(directory,Tot*quota,plot.pred=TRUE,conf.bound=FALSE,years.of.prediction=years.of.prediction,plotfigs = FALSE)
  ind = which(x$total == max(x$total[1:(length(x$total)-years.of.prediction+1)]))
  N70 = 0.7*x$total[ind]
  std=read.standard.deviation("./","wgharp.std",years.of.prediction)
  ind_today = which(x$year == max(x$year-years.of.prediction+1))
  Npred = x$total[ind_today+10]-qnorm(1-0.1)*x$totalSD[ind_today+10]
  out = (N70-Npred)
  return(abs(out))
}

N70.helper.D <- function(Tot,directory,quota)
{
  # Function called from find.N70 that should be used!
  years.of.prediction = 50
  x=run.model(directory,Tot*quota,plot.pred=TRUE,conf.bound=FALSE,years.of.prediction=years.of.prediction,plotfigs = FALSE)
  std=read.standard.deviation("./","wgharp.std",years.of.prediction)
  Dpred = std$Dnew-qnorm(1-0.1)*std$Dnew.std
  out =(0.7-Dpred)
  return(abs(out))
}

N70.helper.mcmc <- function(Tot,directory,quota,nsims,nthin,pos.pred.year)
{
  # Function called from find.N70 that should be used!
  years.of.prediction = 15
  x=run.model(directory,Tot*quota,plot.pred=TRUE,conf.bound=FALSE,years.of.prediction=years.of.prediction,plotfigs = FALSE,mcmc.check=FALSE,nsims,nthin)

  ind = which(x$total == max(x$total[1:(length(x$total)-years.of.prediction+1)]))
  N70 = 0.7*x$total[ind]

  system(paste("./wgharp -mcmc",nsims,"-mcsave",nthin))
  #system("./wgharp -mcmc 1e6 -mcsave 1000")

  system("./wgharp -mceval")

  MCMCSimulations = read.table("MCMC_chains.txt",header = FALSE)
  MCMCSimulations <- as.matrix(MCMCSimulations)
  MCdim = dim(MCMCSimulations)
  Nvar = 2
  N0sims = matrix(NA,MCdim[1],MCdim[2]/Nvar)
  N1sims = matrix(NA,MCdim[1],MCdim[2]/Nvar)

  for (i in 1:MCdim[1]){
		N0sims[i,] = MCMCSimulations[i,1:(MCdim[2]/Nvar)]
		N1sims[i,] = MCMCSimulations[i,(MCdim[2]/Nvar+1):MCdim[2]]
  }

  N0sims = N0sims[-1,]
  N1sims = N1sims[-1,]
  Ntotsims = N0sims + N1sims
  Quantile10 = quantile(Ntotsims[,pos.pred.year],0.10)

  out = (Quantile10 - N70)
  return(abs(out))
}



find.N70.quota <- function(MIN=180000,MAX=230000,quota=c(0.93,0.07),directory="harpeast",method = "Dbased",nsims = 1e6,nthin = 100)
{
  # Function to find equabrilum quota
  quota = quota/sum(quota)
  if (method == "Dbased"){
  tmp = optimize(N70.helper.D,lower=MIN,upper=MAX,directory=directory,quota=quota,tol=50)
  }
  if (method == "Nbased") {
  	tmp = optimize(N70.helper.Nmax,lower=MIN,upper=MAX,directory=directory,quota=quota,tol=50)
  }
  if (method == "mcmc") {
  	pos.year = scan("wgharp.cat",n=1,com="#",quiet=T)+11

  	tmp = optimize(N70.helper.mcmc,lower=MIN,upper=MAX,directory=directory,quota=quota,tol=50,nsims = nsims,nthin = nthin,pos.pred.year = pos.year)
  }
  cat("Equabriluim quota for",directory,"(pups,adults,total):",round(tmp$minimum*quota),sum(round(tmp$minimum*quota)),"\n")

  invisible(tmp$minimum*quota)
}


################### WRITE functions ########################
write.dat <- function(sti="./",
                      A=20,
                      p=c(0.0,0.0,0.0,0.1,0.3,0.5,0.8,0.9,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0) ,
                      K= c(1, 100001, 10000001),
                      M= c(-1, 0.01, 0.12),
                      M0=c(-1, 0.2, 0.8),
                      f= c(-1, 0.7, 1.0))
{

  cat(
       #file=paste(sti,"wgharp.dat",sep=""),
       "#A\n",A,
       "\n\n#P\n",p,"\n\n\n",
       K,"\t#K\n",
       M,"\t#M\n",
       M0,"\t#M0\n",
       f,"\t#f\n",
sep=" ",append=F)
  invisible()
}

write.pin <-function(path="./",
                     K=900001,
                     M=0.10,
                     M0=0.30,
                     f=0.84)
{
  cat(file=paste(path,"wgharp.pin",sep=""),
      K[1],M[1],M0[1],f[1],sep="\n")
  invisible()
}

write.quota <- function(quota=c(0,0),years.of.prediction=50,path="./")
{
  # write input file with catch level for the prediction and years of prediction

  cat(file=paste(path,"wgharp.quo",sep=""),"# Quota file\n",years.of.prediction,"\n",
      formatC(quota[1],format="fg"),formatC(quota[2],format="fg"),append=F)
  invisible(years.of.prediction)
}

write.prior <-function(path="./",
                       K=c(2000000,30000000000),
                       M=c(0.1,0.015),
                       M0=c(0.1,0.015),
                       f=c(0.84,0.2))
{
  cat(file=paste(path,"wgharp.pri",sep=""),
                 "#Prior distribution, Mean, std",
                 "\n#K\n",
                 #formatC(K[1],format="fg")," ",formatC(K[2],format="fg"),
                 K[1]," ",K[2],
                 "\n#M\n",
                 M[1]," ",M[2],
                 "\n#M0\n",
                 M0[1]," ",M0[2],
                 "\n#f\n",
                 f[1]," ",f[2],"\n",
                 sep="")
}



################### READ functions #########################

read.survey.pup.production <- function(path=".\\")
{
  tmp=read.table(file=paste("wgharp.est",sep=""),skip=3)
  colnames(tmp) <- c("year","pup","cv")
  invisible(tmp)
}

read.abundance<-function(path="c:/r/sel/ark03/analyse/",File="wgharp.rep")
{
  # Read abundance trajector for N1, N0 (Adults & pups)
  # colums: year, N1, N0, f

  Data=read.table(paste(path,File,sep=""))
  return(Data)
}

read.standard.deviation<-function(path="c:/r/sel/ark03/analyse/",File,years.of.prediction)
{
  # Colums: index, variable name, value, standard deviation
  # names: V1 ... V4
  # Row:
  # First lines K,M,M0,f if estimated!
  # second last line N0
  # Last line D1+ statistics

  pos.today=scan("wgharp.cat",n=1,com="#",quiet=T)  #in 2003: x$year[pos.today]=2003 (pos.today=58)
  pos.today = pos.today + 1

  tmp=read.table(paste(path,File,sep=""),skip=1)
  ut = list()
  ut$N0 = tmp$V3[5:(3+pos.today+years.of.prediction)]
  ut$N0.std = tmp$V4[5:(3+pos.today+years.of.prediction)]
  ut$N1 = tmp$V3[(5+pos.today+years.of.prediction):(3+2*pos.today+2*years.of.prediction)]
  ut$N1.std = tmp$V4[(5+pos.today+years.of.prediction):(3+2*pos.today+2*years.of.prediction)]
  # N0, N0.2003, D, D,std are temporarly hard coded in wgharp.exe
  ut$N0.2003 = tmp$V3[(length(tmp$V1)-2)]
  ut$N0.2003.std = tmp$V4[(length(tmp$V1)-2)]
  ut$D=tmp$V3[length(tmp$V1)-1]
  ut$D.std=tmp$V4[length(tmp$V1)-1]
  ut$Dnew=tmp$V3[length(tmp$V1)]
  ut$Dnew.std=tmp$V4[length(tmp$V1)]
  #browser()
  return(ut)
}


mcmc.eval <- function(newchain = TRUE,filename = "MCMC_chains.txt",nsims = 1e5,nthin = 10,N70=500,pos.year = 78){
	if(newchain){
		system(paste("./wgharp.exe -mcscale -mcmc",nsims,"-mcsave",nthin))
		system("./wgharp.exe -mceval")
	}
	MCMCSimulations = read.table(filename,header = FALSE)
    MCMCSimulations <- as.matrix(MCMCSimulations)
    MCdim = dim(MCMCSimulations)


	Nvar = 2


	N0sims = matrix(NA,MCdim[1],MCdim[2]/Nvar)
	N1sims = matrix(NA,MCdim[1],MCdim[2]/Nvar)

	for (i in 1:MCdim[1]){
			N0sims[i,] = MCMCSimulations[i,1:(MCdim[2]/Nvar)]
			N1sims[i,] = MCMCSimulations[i,(MCdim[2]/Nvar+1):MCdim[2]]
	}

	N0sims = N0sims[-1,]
	N1sims = N1sims[-1,]

	Ntotsims = N0sims + N1sims



    Quantile10 = quantile(Ntotsims[,pos.year],0.10)
    Quantile90 = quantile(Ntotsims[,pos.year],0.90)

    cat("80% credible interval for year 2023: (",round(Quantile10)," - ",round(Quantile90),")\n")
    if (N70>=Quantile10 & N70<= Quantile90){
    	cat("N70 IS included in 80% credible interval")
    	} else {cat("N70 is NOT included in 80% credible interval")
}

	X11("",11,7)
	par(mar=c(6,5,4,5),bg = "white")

	hist(Ntotsims[,pos.year],50,xlab = "Population size",main = "Posterior distribution of total population in 2023",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
	lines(N70*rep(1,800),1:800,lwd = 3,col ="red")
	lines(Quantile10*rep(1,800),1:800,lwd = 3,lty = 2,col ="blue")
	lines(Quantile90*rep(1,800),1:800,lwd = 3,lty = 2,col ="blue")

	legend("topright",legend = c("80% Credible Interval",expression(N[70])),col = c("blue","red"),lty = c(2,1),lwd = c(3,3),bty = "n",cex = 1.4)

	X11("",11,8)
	par(mar=c(6,5,4,5),bg = "white",mfrow = c(4,1))
	plot(N0sims[,pos.year],xlab = "Realizations",ylab = "Posterior samples",main = "N0 - 2023",pch = ".")
	plot(N0sims[,pos.year-10],xlab = "Realizations",ylab = "Posterior samples",main = "N0 - 2013",pch = ".")
	plot(N0sims[,pos.year-20],xlab = "Realizations",ylab = "Posterior samples",main = "N0 - 2003",pch = ".")
	plot(N0sims[,pos.year-30],xlab = "Realizations",ylab = "Posterior samples",main = "N0 - 1993",pch = ".")


	X11("",11,8)
	par(mar=c(6,5,4,5),bg = "white",mfrow = c(4,1))
	plot(N1sims[,pos.year],xlab = "Realizations",ylab = "Posterior samples",main = "N1 - 2023",pch = ".")
	plot(N1sims[,pos.year-10],xlab = "Realizations",ylab = "Posterior samples",main = "N1 - 2013",pch = ".")
	plot(N1sims[,pos.year-20],xlab = "Realizations",ylab = "Posterior samples",main = "N1 - 2003",pch = ".")
	plot(N1sims[,pos.year-30],xlab = "Realizations",ylab = "Posterior samples",main = "N1 - 1993",pch = ".")


	X11("",11,8)
	par(mar=c(6,5,4,5),bg = "white",mfrow = c(2,1))
     plot(Ntotsims[,pos.year],xlab = "Realizations",ylab = "Posterior samples",main = "Samples from posterior distribution for the year 2023",pch = ".",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
     hist(Ntotsims[,pos.year],50,xlab = "Population size",main = "Posterior distribution of total population in 2023",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5)
	lines(N70*rep(1,800),1:800,lwd = 3,col ="red")
	lines(Quantile10*rep(1,800),1:800,lwd = 3,lty = 2,col ="blue")
	lines(Quantile90*rep(1,800),1:800,lwd = 3,lty = 2,col ="blue")

	legend("topright",legend = c("80% Credible Interval",expression(N[70])),col = c("blue","red"),lty = c(2,1),lwd = c(3,3),bty = "n",cex = 1.4)

}

#Means square Error
MSE.fit <- function(Xmat)
{
	Data <- read.survey.pup.production(".")

	MSEest = 0
	for(i in 1:length(Data)){
		i1 = which(Xmat[,1] == Data[i,1])
		MSEest = MSEest + (Xmat[i1,5]-Data[i,2])^2/Data[i,3]^2
		}
	return(MSEest)
	}

plot.mod.fit <- function(X,yr1 = 1995,yr2 = 2015,col = "red",lwd = 4,newplot = FALSE){

 scala = 100000

 est = read.survey.pup.production()

 Nest = dim(est)
 est$lower = est$pup*(1-1.96*est$cv)
 est$upper = est$pup*(1+1.96*est$cv)



if(newplot){
	X11("",10,7)
	par(mar=c(6,5,4,5),bg = "white",mfrow = c(1,1))
	plot(est$year,est$pup/scala,type = "p",col = "blue",pch=16,xlim = c(yr1,yr2),ylim = c(0,6),xlab = "Year",ylab = "Abundance (in 10 000)",cex.axis = 1.5,cex.lab = 1.5,cex.main = 1.5,bty = "l")
 }


  for (i in 1:Nest[1]){
  	  lines(rep(est$year[i],10),seq(est$pup[i]/scala+0.03,est$upper[i]/scala,length.out = 10),type = "l",col = "blue")
  	  lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$upper[i]/scala,10),type = "l",col = "blue")
	  lines(rep(est$year[i],10),seq(est$pup[i]/scala-0.03,est$lower[i]/scala,length.out = 10),type = "l",col = "blue")
	  lines(seq(est$year[i]-0.3,est$year[i]+0.3,length.out = 10),rep(est$lower[i]/scala,10),type = "l",col = "blue")
  	}
  	ind1 = which(X$year == yr1)
  	ind2 = which(X$year == yr2)
  	lines(X$year[ind1:ind2],X$N0[ind1:ind2]/scala,col = col,lwd = lwd)
}


buildPandF <- function(Fproj = Fproj)
{
  ##Prepare fecundity data
  Fdat <- read.table("wgharp.fdt",header = FALSE,sep = "")

  cdata <- scan("wgharp.cat",what="raw")
  yr1 = 1946
  yr2 = as.numeric(cdata[length(cdata)-2])
  nyr = yr2 - yr1 + 1

  Fvec = array(0,nyr)
  Fvec[1:(Fdat$V1[1] - yr1 + 1)] = Fdat$V2[1]

  for (i in 2:length(Fdat$V1)){
    if(Fdat$V1[i]-Fdat$V1[i-1]==1){
      Fvec[Fdat$V1[i]-yr1+1] = Fdat$V2[i]
    } else {
      i1 = Fdat$V1[i-1] - yr1 + 1
      i2 = Fdat$V1[i] - yr1 + 1
      Fvec[i1:i2] = seq(Fdat$V2[i-1],Fdat$V2[i],length.out = (i2-i1+1))
    }
  }
  Fvec[i2:nyr] = Fdat$V2[i]

  if(class(Fproj) == "character") Fvec[(yr2-yr1+1):nyr] = mean(Fdat$V2)

  if(class(Fproj)=="numeric") Fvec[(yr2-yr1+1):nyr] = Fproj

  write.table(Fvec,"wgharp.f",row.names = FALSE, col.names = FALSE)

##########################
                                        #Prepare birth ogive data
  Pdat <- read.table("wgharp.ogi",sep = "",header = TRUE)
  Pper <- read.table("wgharp.ogp",sep = "",header = TRUE)
  P = matrix(0,nyr,max(Pdat$Age))

  for(i in 1:(Pper$Pstop[1]-yr1+1)){
    P[i,] = Pdat[,2]
  }

  for(i in 2:length(Pper$Pstart)){
    i1 = Pper$Pstop[i-1]-yr1+1
    i2 = Pper$Pstart[i]-yr1+1
    i3 = Pper$Pstop[i]-yr1+1

    for(j in i2:i3){
      P[j,] = Pdat[,i+1]
    }

    for(age in 1:max(Pdat$Age)){
      P[i1:i2,age] = seq(P[i1,age],P[i2,age],length.out = (i2-i1+1))
    }

  }

  for(i in i3:nyr){
    P[i,] = Pdat[,length(Pper$Pstart)+1]
  }

  write.table(P,"wgharp.pma",row.names = FALSE, col.names = FALSE)

}
