#' Plot estimated population sizes from optimized population model for harp seals and hooded seals
#'
#' Plot population trajectories.
#' @param results Output from fitted model
#' @param dat Original data on estimated pup production (set to NA if these are not to be included).
#' @param component Which population component to plot. Can be either 'N0' for pups, 'N1' for adults, or a vector with both.
#' @param xLim Manually set the x axis extent, overrides the default which is the extent of the plotted data.
#' @param yLim Manually set the y axis extent, overrides the default which is the extent of the plotted data.
#' @param plot.legend Add legend to plot (TRUE/FALSE)
#' @param plot.Nlims True/False
#' @param projections Plot projections (TRUE/FALSE)
#' @param mean.proj True/False
#' @param width Width of the figure
#' @param height Height of the figure
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plot.N(res, data)

plot.N <- function(results=res,dat=data,component=c('N0', 'N1'),
                   xLim=NA,yLim=NA,plot.legend=TRUE,
                   plot.Nlims=TRUE, projections=TRUE, 
                   mean.proj=FALSE,width = 13, height = 7,
                   conf.int = TRUE)
{
  if(projections) {
    span <- c(1:length(results$indN0))
  } else {
    span <- c(1:nrow(dat$Cdata))
  }
  
  options(scipen=999)
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  windows("",width = width,height = height)
  
  #require(RColorBrewer)
  theCols <- RColorBrewer::brewer.pal(max(c(3, length(component))), 'Dark2')
  
  if(length(xLim)==1) xLim <- range(results$years)
  
  if('N1' %in% component) {
    indN1back <- results$indN1[c(1:match(max(dat$pupProductionData[,1]), results$years))]
    lCI <- results$rep.matrix[results$indN1,1]-(1.96*results$rep.matrix[results$indN1,2])
    uCI <- results$rep.matrix[results$indN1,1]+(1.96*results$rep.matrix[results$indN1,2])
    
    if(length(yLim)==1) yLim <- c(0, max(uCI))
    
    if(!projections) {
      if(length(xLim)==1) {
        xLim <- range(results$years[span])
      } else {
        xLim[2] <- min(c(xLim[2], max(dat$Cdata[,1])))
      }  
      if(length(yLim)==1) {
        yLim <- c(0, max(uCI[span]))
      } else {
        yLim[2] <- max(uCI[span])
      }  
    }
    par(mar=c(5,6,4,2)+0.1)
    
    plot(results$years, results$rep.matrix[results$indN1,1],
         type = "l",xlab = "Year",ylab = "Abundance",
         col=theCols[1], lwd=2, lty=2, xlim=xLim, ylim=yLim, axes=FALSE,cex.lab = 1.5)
    if(!projections | !mean.proj) lines(results$years, results$rep.matrix[results$indN1,1], lwd=2, lty=2, col='white')
    axis(1)
    axis(2, at=pretty(par('usr')[c(3,4)]), 
         labels=format(pretty(par('usr')[c(3,4)]), scientific=F))
    abline(h=par('usr')[3])
    abline(v=par('usr')[1])
    if(plot.Nlims) {
      Nlims <- c(0.3, 0.5, 0.7)*max(results$rep.matrix[results$indNTot,1]) 
      abline(h=Nlims, col='lightgrey')
      text(par('usr')[2], Nlims[1], expression(N[lim]), xpd=NA, adj=0, cex=0.9)
      text(par('usr')[2], Nlims[2], expression(N[50]), xpd=NA, adj=0, cex=0.9)
      text(par('usr')[2], Nlims[3], expression(N[70]), xpd=NA, adj=0, cex=0.9)
    }
    
    if(projections) {
      if(conf.int){
        lines(results$years, lCI, col=theCols[1], lty=2)
        lines(results$years, uCI, col=theCols[1], lty=2)
      }
      lines(results$years, results$rep.matrix[results$indN1,1],col=theCols[1], lwd = 2, lty=2)
    }
    
    if(conf.int){
      polygon(c(c(results$years[1]:max(dat$pupProductionData[,1])), 
              rev(c(results$years[1]:max(dat$pupProductionData[,1])))), 
            c(lCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))],
              rev(uCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))])),
            col=add.alpha(theCols[1], 0.5), border=NA)
    }
     
    lines(c(results$years[1]:max(dat$pupProductionData[,1])), results$rep.matrix[indN1back,1],
          col=theCols[1], lwd=2)
    
  } 
  
  if('N0' %in% component) {
    indN0back <- results$indN0[c(1:match(max(dat$pupProductionData[,1]), results$years))]
    plCI <- results$rep.matrix[results$indN0,1]-(1.96*results$rep.matrix[results$indN0,2])
    puCI <- results$rep.matrix[results$indN0,1]+(1.96*results$rep.matrix[results$indN0,2])
    if(!'N1' %in% component) {
      if(length(yLim)==1) yLim <- c(0, max(puCI))
      if(!projections) {
        if(length(xLim)==1) {
          xLim <- range(results$years[span])
        } else {
          xLim[2] <- min(c(xLim[2], max(dat$Cdata[,1])))
        }  
        if(length(yLim)==1) {
          yLim <- c(0, max(c(max(puCI[span]), max(dat$pupProductionData[,2]+(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2]))))))
        } else {
          yLim[2] <- max(c(max(puCI[span]), max(dat$pupProductionData[,2]+(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])))))
        }  
      }
      plot(results$years, results$rep.matrix[results$indN0,1],
           type = "l",xlab = "Year",ylab = "Abundance",
           col=theCols[2], lwd=2, lty=2, xlim=xLim, ylim=yLim, axes=FALSE,cex.lab = 1.5)
      if(!projections | !mean.proj) lines(results$years, results$rep.matrix[results$indN0,1], lwd=2, lty=2, col='white')
      axis(1)
      axis(2, at=pretty(par('usr')[c(3,4)]), 
           labels=format(pretty(par('usr')[c(3,4)]), scientific=FALSE))
      abline(h=par('usr')[3])
      abline(v=par('usr')[1])
    } else { 
      if(projections) lines(results$years, results$rep.matrix[results$indN0,1],
                            col=theCols[2], lwd=2,lty=2)
    }
    if(projections) {
      if(conf.int){
        lines(results$years, plCI, col=theCols[2], lty=2)
        lines(results$years, puCI, col=theCols[2], lty=2)
      }
    }  
    
    if(conf.int){
    polygon(c(c(results$years[1]:max(dat$pupProductionData[,1])), 
              rev(c(results$years[1]:max(dat$pupProductionData[,1])))), 
            c(plCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))],
              rev(puCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))])),
            col=add.alpha(theCols[2], 0.5), border=NA)
    }
    
    lines(c(results$years[1]:max(dat$pupProductionData[,1])), results$rep.matrix[indN0back,1],
          col=theCols[2], lwd=2)
    
  } 
  
  if(length(dat)>1 & 'N0' %in% component) {
    segments(dat$pupProductionData[,1], 
             dat$pupProductionData[,2]-(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])),
             dat$pupProductionData[,1],
             dat$pupProductionData[,2]+(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])),
             col=1)
    segments(dat$pupProductionData[,1]-0.5, 
             dat$pupProductionData[,2]-(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])),
             dat$pupProductionData[,1]+0.5,
             dat$pupProductionData[,2]-(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])),
             col=1)
    segments(dat$pupProductionData[,1]-0.5, 
             dat$pupProductionData[,2]+(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])),
             dat$pupProductionData[,1]+0.5,
             dat$pupProductionData[,2]+(1.96*(dat$pupProductionData[,3]*dat$pupProductionData[,2])),
             col=1)
    
    points(dat$pupProductionData[,1],dat$pupProductionData[,2],
           pch=21, bg=theCols[2], cex=1.5)
  }
  
  if(plot.legend) {
    if('N1' %in% component) {
      if(length(component)==1) {
        legend('topright', lwd=2, col=theCols[1], 
               pch=NA, 
               '1+ population size',
               bty='n')
      }  else {
        legend('topright', lwd=rep(2, 2), col=theCols[c(1,2)], 
               c('1+ population size', 
                 'Pup production'), pch=rep(NA, 2),
               bty='n')
        
        if(length(dat)>1) {
          legend('topright', lwd=rep(1, 2), lty=NA, 
                 pt.bg=c(NA, theCols[2]), pch=c(NA,21),
                 pt.cex=c(0, 1.5),
                 c('1+ population size', 
                   'Pup production'),
                 bty='n')
        }
      }
    } else {
      legend('topright', lwd=2, col=theCols[2], 
             'Pup production',
             bty='n')
      if(length(dat)>1) {
        legend('topright', lwd=1, lty=NA, 
               pt.bg=theCols[2], pch=21, 
               pt.cex=1.5, 'Pup production',
               bty='n')
      }
    }  
  }
  options(scipen=0)
}



#' Plot modelled population dynamics
#'
#' Plot population trajectories. This function does not work properly, use plot.N() instead.
#' @param results Output from fitted model
#' @param dat Original data on estimated pup production (set to NA if these are not to be included).
#' @param component Which population component to plot. Can be either 'N0' for pups, 'N1' for adults, or a vector with both.
#' @param xLim Manually set the x axis extent, overrides the default which is the extent of the plotted data.
#' @param yLim Manually set the y axis extent, overrides the default which is the extent of the plotted data.
#' @param plot.legend Add legend to plot (TRUE/FALSE)
#' @param plot.Nlims Plot horizontal lines indicating 30%, 50% and 70% of 1+ population True/False
#' @param col Specify colors used for (pups,1+,error bars)
#' @param projections Plot projections (TRUE/FALSE)
#' @param add.title Add title to figure (TRUE/FALSE)
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plot.res(res, data)

plot.res.ggplot <- function(res=res,dat=data,component=c('N0', 'N1'),
                   xLim=NA,yLim=NA,plot.legend=TRUE,
                   plot.Nlims=TRUE, 
                   col = c("royalblue","darkred","steelblue"),
                   projections=TRUE, 
                   conf.int = TRUE,
                   add.title = TRUE,
                   width = 13,
                   height = 7)
{
  require(ggplot2)
  
  if(projections){
  }
  
  dframe <- data.frame(Year = res$years,
                       N0 = as.vector(res$rep.matrix[res$indN0,1]),
                       N0sd = as.vector(res$rep.matrix[res$indN0,2]),
                       N1 = as.vector(res$rep.matrix[res$indN1,1]),
                       N1sd = as.vector(res$rep.matrix[res$indN1,2])
                       )
  
  
  dframe$Ntot = dframe$N0 + dframe$N1
  dframe$Ntotsd = sqrt(dframe$N0sd^2+dframe$N1sd^2)
  dframe$N0LL = dframe$N0 - 1.96*dframe$N0sd
  dframe$N0UL = dframe$N0 + 1.96*dframe$N0sd
  dframe$N1LL = dframe$N1 - 1.96*dframe$N1sd
  dframe$N1UL = dframe$N1 + 1.96*dframe$N1sd
  dframe$NtotLL = dframe$Ntot - 1.96*dframe$Ntotsd
  dframe$NtotUL = dframe$Ntot + 1.96*dframe$Ntotsd
  
  dfN0 = data.frame(Year = dframe$Year,Abun = dframe$N0,id = 'Pup population')
  dfN1 = data.frame(Year = dframe$Year,Abun = dframe$N1,id = '1+ population')
  dfNtot = data.frame(Year = dframe$Year,Abun = dframe$Ntot,id = "Total population")
  dfPupsbig = data.frame(Year = dframe$Year,Abun = NA,id = "Pup production estimates")
  dfPupsbig$Abun[which(dframe$Year %in% dfpups$Year)] = data$pupProductionData[,2]
  
  #Find out which components to plot
  if(length(component)>1){
    if(("N0" %in% component) & ("N1" %in% component))   dfall <- rbind.data.frame(dfN0,dfN1)
    if(("N0" %in% component) & ("Ntot" %in% component))   dfall <- rbind.data.frame(dfN0,dfNtot)
  }
  
  
  dfpups <- data.frame(Year = data$pupProductionData[,1],Pups = data$pupProductionData[,2],sd = data$pupProductionData[,2]*data$pupProductionData[,3])
  dfpups$LL = dfpups[,"Pups"]-1.96*dfpups[,"sd"]
  dfpups$UL = dfpups[,"Pups"]+1.96*dfpups[,"sd"]
  
  maxYear = max(dframe$Year)
  maxN = max(dframe$N1)
  
  pl <- ggplot(data = dframe,aes(x = Year))
  if(plot.Nlims){
    pl <- pl + geom_segment(aes(x=min(dframe$Year),xend=maxYear,y=(0.3*maxN),yend=(0.3*maxN)),color = "lightgrey",size = 0.5)
    pl <- pl + geom_segment(aes(x=min(dframe$Year),xend=maxYear,y=(0.5*maxN),yend=(0.5*maxN)),color = "lightgrey",size = 0.5)
    pl <- pl + geom_segment(aes(x=min(dframe$Year),xend=maxYear,y=(0.7*maxN),yend=(0.7*maxN)),color = "lightgrey",size = 0.5)
    
    #  geom_hline(yintercept = 0.3*max(dframe$N1),color = "lightgrey",size = 0.5)
    #pl <- pl + geom_hline(yintercept = 0.5*max(dframe$N1),color = "lightgrey",size = 0.5)
    #pl <- pl + geom_hline(yintercept = 0.7*max(dframe$N1),color = "lightgrey",size = 0.5)
    pl <- pl + geom_text(x = (max(dframe$Year)+10),y = 0.3*max(dframe$N1),label = (expression(N[lim])),size = 5)
    pl <- pl + geom_text(x = (max(dframe$Year)+10),y = 0.5*max(dframe$N1),label = (expression(N[50])), size = 5)
    pl <- pl + geom_text(x = (max(dframe$Year)+10),y = 0.7*max(dframe$N1),label = (cexpression(N[70])), size = 5)
  }
  
  # Add confidence intervals
  if(conf.int){
    if("N0" %in% component) pl <- pl  + geom_ribbon(data=dframe,aes(ymin=N0LL,ymax=N0UL,x=Year),fill = col[1],alpha = 0.3)
    if("N1" %in% component) pl <- pl  + geom_ribbon(data=dframe,aes(ymin=N1LL,ymax=N1UL,x=Year),fill = col[2],alpha = 0.3)
    if("Ntot" %in% component) pl <- pl  + geom_ribbon(data=dframe,aes(ymin=NtotLL,ymax=NtotUL,x=Year),fill = col[2],alpha = 0.3)
  }


  #If only one component is plottet, plot the right one
  if(length(component) == 1){
    if("N0" %in% component) pl  <- pl + geom_line(data = dframe,aes(x = Year,y = N0),color = col[1],size = 1.5)
    if("N0" %in% component) pl  <- pl + geom_point(data = dfpups,aes(x=Year,y = Pups,col = "red"),size = 3,color = col[3]) 
    if("N0" %in% component) pl  <- pl + geom_errorbar(data =dfpups,aes(x = Year,ymin = LL,ymax = UL), color = col[3],size = 1)                     
    if("N1" %in% component) pl  <- pl + geom_line(data = dframe,aes(x = Year,y = N1),color = col[2],size = 1.5) 
    if("Ntot" %in% component) pl  <- pl + geom_line(data = dframe,aes(x = Year,y = Ntot),color = col[2],size = 1.5) 
  } else {
        #Plot multiple lines
        pl <- pl + geom_line(data = dfall,aes(x = Year,y = Abun,colour = id),size = 1.5)
        if(("N0" %in% component) & ("N1" %in% component)) pl <- pl + scale_colour_manual(values = c('Pup population' = col[1],'1+ population' = col[2]))
        if(("N0" %in% component) & ("Ntot" %in% component)) pl <- pl + scale_colour_manual(values = c('Pup population' = col[1],'Total population' = col[2]))
        
  }
  
  #pl <- pl + scale_linetype_manual(labels = c("psavert", "uempmed"),values = c("psavert"="#00ba38", "uempmed"="#f8766d"))
  
  #pl <- pl + scale_colour_manual(values = c("purple", "green", "blue"),
  #                               guide = guide_legend(override.aes = list(
  #                                 linetype = c("solid","dashed"),
  #                                 shape = c(NA,NA))))
  
  #pl <- pl + expand_limits(x = (max(dframe$Year)+5))
  
  #pl <- pl + scale_fill_manual(labels = paste("long", c(5, 10, 15)),
  #                      guide = guide_legend(
  #                        direction = "horizontal",
  #                        title.position = "top",
  #                        label.position = "bottom",
  #                        label.hjust = 0.5,
  #                        label.vjust = 1,
  #                        label.theme = element_text(angle = 90)
  #                      ))
  
  pl <- pl + coord_cartesian(xlim = c(min(dframe$Year), (max(dframe$Year)+5)), clip = 'off')
  pl <- pl + theme_classic() 
  pl <- pl + theme(text = element_text(size=20),
                   plot.margin = unit(c(1,2,1,1), "cm"),
                   axis.text.y = element_text(angle = 90,margin = margin(t = 0, r = 20, b = 0, l = 30),vjust = 0.5),
                   axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),vjust = 1),
                   #legend.title = element_text(size = 20),
                   legend.title = element_blank(),
                   legend.position = "bottom",
                   legend.box = "vertical"
                   )
  pl <- pl + ylab("Population size")
  
  pl <- pl + labs(fill = "Dose (mg)")
  
  pl <- pl + scale_fill_manual(name = "Dose", labels = c("A", "B", "C"))
  
  #Add title if needed
  if(add.title)  pl <- pl + labs(title="Modelled population dynamics")
  

  #windows("",width = width,height = height)
  pl
  
  }


#' Plot the reported catch data
#'
#' Plot the reported catch data.
#' @param catch Reported catcj data
#' @param position Position of bars: If Pup catch and 1+ catch next to each other use position = "dodge" (default). On top of each other use position = "stack"
#' @param width Figure width
#' @param height Figure height
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plot.catch(data$Cdata)

plot.catch <- function(catch = cdata,width = 9,height = 7,position = "dodge")
{
  library(ggplot2)
  
  theCols <- RColorBrewer::brewer.pal(3, 'Dark2')
  
  dfpups = data.frame("Year" = catch[,1],"Catches" = catch[,2],id = "Pup catch")
  dfOnePluss = data.frame("Year" = catch[,1],"Catches" = catch[,3],id = "1+ catch")
  
  dfcatch = rbind.data.frame(dfpups,dfOnePluss)
  
  #pl <- ggplot(data = dfcatch,aes(x = Year,fill = id))
  
  #pl <- pl + geom_col(aes(y = Catches),position = "dodge")
  
  pl <- ggplot2::ggplot(data = dfcatch, aes(x = Year))
  pl <- pl + ggplot2::geom_col(aes(y = Catches/1000, fill = id),position = position)
  #pl <- pl + scale_colour_discrete(values = c('Pup catch' = "red",'1+ catch' = "blue"))
  pl <- pl + ggplot2::theme_classic() 
  pl <- pl + ggplot2::theme(text = element_text(size=20),
                   plot.margin = unit(c(1,2,1,1), "cm"),
                   axis.text.y = element_text(angle = 90,margin = margin(t = 0, r = 20, b = 0, l = 30),vjust = 0.5),
                   axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),vjust = 1),
                   #legend.title = element_text(size = 20),
                   legend.title = element_blank(),
                   #legend.position = "bottom",
                   #legend.box = "vertical"
  )
  pl <- pl + ggplot2::ylab("Catch level (in 1000)")
  pl <- pl + ggplot2::scale_fill_manual(values = c(theCols[1], theCols[2]))
  #pl <- pl + labs(fill = "Dose (mg)")
  
  #pl <- pl + scale_fill_manual(name = "Dose", labels = c("A", "B", "C"))
  
  
  #windows("",width = 15,height = 7)
  pl
  
}


#' Plot fecundity data
#'
#' Plot the fecundity data used in the model fit.
#' @param data Data object used in model fit
#' @param include.observations Plot the observed fecundity rates with 95% confidence intervals (default = TRUE)
#' @param population If include.observations is TRUE, define which population you want to use.

#' @return plot Returns a plot of predicted population size for different population components
#' @keywords fecundity data
#' @export
#' @examples
#' plot.fecundity(res, data)

plot.fecundity <- function(dat = data,include.observations = TRUE, population = "harpeast")
{


  par(mar=c(5,6,4,2)+0.1)
  plot(dat$Cdata[,1],dat$Ftmp,type = "l",lwd = 3, col = "royalblue", xlab='Year', ylab='Fecundity rate',ylim=c(0,1), axes=F,cex.lab=1.5)
  if(include.observations){
    fecundity <- read.table(paste("Data/",population,"/fecundity.dat",sep = ""),header = FALSE)
    segments(fecundity[,1], 
             fecundity[,2]-(1.96*(fecundity[,3]*fecundity[,2])),
             fecundity[,1],
             fecundity[,2]+(1.96*(fecundity[,3]*fecundity[,2])),
             col="steelblue")
    segments(fecundity[,1]-0.5, 
             fecundity[,2]-(1.96*(fecundity[,3]*fecundity[,2])),
             fecundity[,1]+0.5,
             fecundity[,2]-(1.96*(fecundity[,3]*fecundity[,2])),
             col="steelblue")
    segments(fecundity[,1]-0.5, 
             fecundity[,2]+(1.96*(fecundity[,3]*fecundity[,2])),
             fecundity[,1]+0.5,
             fecundity[,2]+(1.96*(fecundity[,3]*fecundity[,2])),
             col="steelblue")
    
    points(fecundity[,1],fecundity[,2],
           pch=21, bg="steelblue", cex=1.5)
    
    
  }
  
  axis(1)
  axis(2, at=pretty(par('usr')[c(3,4)]), 
       labels=format(pretty(par('usr')[c(3,4)]), scientific=FALSE))
  abline(h=par('usr')[3])
  abline(v=par('usr')[1])
  
  if(include.observations){
    legend('bottomright', legend = c("Fecundity used in model","Observed fecundity"),lty = c(1,1),pch = c(NA,19),lwd=c(3,2),bg = c(NA,"steelblue"),col = c("royalblue","steelblue"), bty='n', cex=1.3)
  }
    #box()
}

#' Create table with key parameters from fitted population model 
#'
#' Parameter table
#' @param results Output from fitted model
#' @param dat Original data on estimated pup production (set to NA if these are not to be included).
#' @return table Returns a table with key parameters from fitted population model.
#' @keywords population model
#' @export
#' @examples
#' par.table()

par.table <- function(results=res, dat=data, tab2flex=FALSE) {
  ## CHECK THIS FUNCTION!!
  means <- c(as.vector(results$Kest), 
             as.vector(results$M0est), 
             as.vector(results$Mest),
             results$N0Current, 
             results$N1Current,
             results$NTotCurrent,
             results$D1, results$DNmax)
  sds <- c(round(results$rep.matrix[results$indN1[1],2]), 
           results$rep.matrix[3,2],
           results$rep.matrix[2,2], 
           results$N0Current.sd, 
           results$N1Current.sd, 
           results$NTotCurrent.sd, 
           results$D1.sd, results$DNmax.sd)
  
  parNames <- c(paste0('N', results$years[1]),
                'M0', 'M1+', 
                paste0('N0,', max(dat$Cdata[,1])),
                paste0('N1+,', max(dat$Cdata[,1])),
                paste0('NTotal,', max(dat$Cdata[,1])),
                'D1+', paste0('NTotal,', max(results$years)))
                
  if(tab2flex) {
    means <- unlist(lapply(means, function(x) {
      ifelse(x<10, format(x, digits=2, scientific=F), 
             format(x, big.mark=' ', digits=0, scientific=F))
    }))
    sds <- unlist(lapply(sds, function(x) {
      ifelse(x<10, format(x, digits=2, scientific=F), 
             format(x, big.mark=' ', digits=0, scientific=F))
    }))
  } else {
    data.frame(par=c('Kest', 'M0est', 'M1est', 'N0Current', 'N1Current', 'NTotCurrent', 'D1', 'NTotprojected'), 
               parNames, Mean=means, SD=sds, stringsAsFactors = F)
  }
}

#' Plot birth ogive curves
#'
#' Plot the birth ogive curves for various time periods.
#' @param dat input data for the model
#' @param highlight.last Highlist last observed birth ogive curve (TRUE/FALSE)
#' @keywords birth ogive
#' @export
#' @examples
#' plot.Pmat()
plot.Pmat <- function(dat=data, highlight.last=TRUE) {
  opal <- palette()
  palette(RColorBrewer::brewer.pal(8, 'Dark2'))
  par(mar=c(5,6,4,2)+0.1)
  matplot(t(dat$Pmat), type='l', lty=1, col='grey',
          xlab='Age (years)',
          ylab='Proportion of mature females',
          ylim=c(0,1), axes=F,cex.lab=1.5)
  axis(1,at=seq(0, 20, by=2))
  axis(2, at=seq(0, 1, by=0.2))
  #box()
  ymat <- match(dat$Pper$Pstart, dat$Cdata[,1])
  matlines(t(dat$Pmat[ymat,]), lty=1, col=c(1:length(ymat)), lwd=2)
  if(highlight.last) {
    lines(c(1:dim(dat$Pmat)[2]), dat$Pmat[tail(ymat, 1),], 
                           lty=1, lwd=3, col='black')
  }  
  leglab <- apply(dat$Pper, 1, function(x) {
    if(x[1]==x[2]) {
      x[1]
    } else {
      paste(x, collapse='-')
    }  
  })
  
  if(highlight.last) {
    legend('bottomright', lwd=c(rep(2, length(ymat)-1), 3, 1), 
           col=c(c(1:(length(ymat)-1)), 'black', 'grey'), 
           c(leglab, 'Between periods'), bty='n', cex=1.1)
  } else {
    legend('bottomright', lwd=c(rep(2, length(ymat)), 1), 
           col=c(c(1:length(ymat)), 'grey'), 
           c(leglab, 'Between periods'), bty='n', cex=1.1)
  }
  palette(opal)       
}

