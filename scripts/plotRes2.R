plotRes2 <- function(results=res,
                    data=data,
                    component=c('N0', 'Ntot'),
                    plotNlims=TRUE, 
                    plotProjections=TRUE,
                    plotProjMean = TRUE,
                    width = 15, 
                    height = 10,
                    grDev = TRUE)
{
  
  
  
  # if(plotProjections) {
  #   span <- c(1:length(results$indN0))
  # } else {
  #   span <- c(1:nrow(data$Cdata))
  # }
  
  span = c(1:nrow(data$Cdata))
  
  
  NrowsDf = length(component)*length(span)
  
  df = data.frame(Year = rep(results$years[span],length(component)),
                  Abundance = NA,
                  LL=NA,
                  UL=NA,
                  group = NA)
  
  dfIndex = span
  
  
  #----------------------
  #Deal with projections
  #----------------------  
  if(plotProjections){
    spanProj = c(1:length(results$indN0))
    indProj = setdiff(spanProj,span)
    dfProjIndex = 1:length(indProj)
    dfProj = data.frame(Year=rep(results$years[indProj],length(component)),
                        Abundance = NA,
                        LL = NA,
                        UL = NA,
                        group = NA)
  }
  
  #--------------------------------
  #Adding the N0 (pups) abundance
  #--------------------------------
  if("N0" %in% component){
    df$Abundance[dfIndex] = results$rep.matrix[results$indN0[span],1]
    df$group[dfIndex] = "0 group"
    df$LL[dfIndex] = results$rep.matrix[results$indN0[span],1] - (1.96*results$rep.matrix[results$indN0[span],2]) 
    df$UL[dfIndex] = results$rep.matrix[results$indN0[span],1] + (1.96*results$rep.matrix[results$indN0[span],2]) 
    
    #Estimated pup production estimates
    dfPups = data.frame(Year = data$pupProductionData[,1],
                        pupEst = data$pupProductionData[,2],
                        LL = (data$pupProductionData[,2] - (1.96*(data$pupProductionData[,3]*data$pupProductionData[,2]))),
                        UL = (data$pupProductionData[,2] + (1.96*(data$pupProductionData[,3]*data$pupProductionData[,2]))))

    dfIndex = (tail(span,n=1)+1):(2*tail(span,n=1))
    
    if(plotProjections){
      dfProj$Abundance[dfProjIndex] = results$rep.matrix[results$indN0[indProj],1]
      dfProj$group[dfProjIndex] = "0 group"
      dfProj$LL[dfProjIndex] = results$rep.matrix[results$indN0[indProj],1] - (1.96*results$rep.matrix[results$indN0[indProj],2]) 
      dfProj$UL[dfProjIndex] = results$rep.matrix[results$indN0[indProj],1] + (1.96*results$rep.matrix[results$indN0[indProj],2])
      dfProjIndex = (tail(dfProjIndex,n=1)+1):(tail(dfProjIndex,n=1)+length(indProj))
    }
    
  }
  
  
  #---------------------------
  # Adding 1+ abundance
  if("N1" %in% component){
    
    df$Abundance[dfIndex] = results$rep.matrix[results$indN1[span],1]
    df$group[dfIndex] = "1+ group"
    df$LL[dfIndex] = results$rep.matrix[results$indN1[span],1] - (1.96*results$rep.matrix[results$indN1[span],2]) 
    df$UL[dfIndex] = results$rep.matrix[results$indN1[span],1] + (1.96*results$rep.matrix[results$indN1[span],2]) 
    
    dfIndex = (tail(dfIndex,n=1)+1):((tail(dfIndex,n=1)+length(span)))
    
    if(plotProjections){
      dfProj$Abundance[dfIndex] = results$rep.matrix[results$indN1[indProj],1]
      dfProj$group[dfIndex] = "1+ group"
      dfProj$LL[dfIndex] = results$rep.matrix[results$indN1[indProj],1] - (1.96*results$rep.matrix[results$indN1[indProj],2]) 
      dfProj$UL[dfIndex] = results$rep.matrix[results$indN1[indProj],1] + (1.96*results$rep.matrix[results$indN1[indProj],2]) 
      dfProjIndex = (tail(dfProjIndex,n=1)+1):(tail(dfProjIndex,n=1)+length(indProj))
    }
  }
  
  
  #-----------------------------
  # Adding total abundance
  #-----------------------------
  if("Ntot" %in% component){
    
    df$Abundance[dfIndex] = results$rep.matrix[results$indNTot[span],1]
    df$group[dfIndex] = "All age groups"
    df$LL[dfIndex] = results$rep.matrix[results$indNTot[span],1] - (1.96*results$rep.matrix[results$indNTot[span],2]) 
    df$UL[dfIndex] = results$rep.matrix[results$indNTot[span],1] + (1.96*results$rep.matrix[results$indNTot[span],2]) 
    
    if(plotProjections){
      dfProj$Abundance[dfProjIndex] = results$rep.matrix[results$indNTot[indProj],1]
      dfProj$group[dfProjIndex] = "All age groups"
      dfProj$LL[dfProjIndex] = results$rep.matrix[results$indNTot[indProj],1] - (1.96*results$rep.matrix[results$indNTot[indProj],2]) 
      dfProj$UL[dfProjIndex] = results$rep.matrix[results$indNTot[indProj],1] + (1.96*results$rep.matrix[results$indNTot[indProj],2]) 
    }
    
  }
  
  
  #-----------------------------------------------
  # Preparing to plot the N30, N50, and N70 lines
  #-----------------------------------------------
  
  # Do not add lines if only pup abundance is plotted, if so force plotNlims to be FALSE
  if((("N0" %in% component) & length(component) == 1)) plotNlims = FALSE
  
  if(plotNlims){
    Nlims <- c(0.3, 0.5, 0.7)*max(results$rep.matrix[results$indNTot,1]) 
    dfNlims = data.frame(Year = results$years[span],N30 = Nlims[1],N50 = Nlims[2], N70 = Nlims[3])
  }



theCols <- RColorBrewer::brewer.pal(max(c(3, length(component))), 'Dark2')
if(grDev) graphDev(width = width,height = height)

scalef = 10000

# Set the color of the pup production esitmates
if("N0" %in% component){
  if(length(component) == 1) colPupest = theCols[3]
  if(length(component) == 2) colPupest = theCols[3]
  if(length(component) == 3) colPupest = theCols[2]
}

p1 <- ggplot() + 
  geom_line(data = df,
            aes(x=Year,
                y=Abundance/scalef,
                group = group,
                color = group),
            size = 1.3,
            linetype = 1) +
  geom_ribbon(data=df,
              aes(x = Year, ymin=LL/scalef,ymax=UL/scalef, 
                  group = group,
                  color = NULL,
                  fill = group),
              alpha=0.3)

if(plotProjections){
  if(plotProjMean){
    p1 <- p1 + geom_line(data = dfProj,
                         aes(x=Year,
                             y=Abundance/scalef,
                             group = group,
                             color = group),
                         size = 0.8,
                         linetype = 2)
  }
    
  p1 <- p1 + geom_line(data = dfProj,
              aes(x=Year,
                  y=LL/scalef,
                  group = group,
                  color = group),
              size = 0.8,
              linetype = 2) +
    geom_line(data = dfProj,
              aes(x=Year,
                  y=UL/scalef,
                  group = group,
                  color = group),
              size = 0.8,
              linetype = 2) +
    xlim(NA,results$years[c(length(results$indN0))])
   # + geom_vline(xintercept = tail(data$Cdata[,1],n=1), 
   #              linetype="dotted", 
   #              color = "grey", 
   #              size=0.8) + 
   #   annotate("text", x = (tail(data$Cdata[,1],n=1)+10), y = max(df$UL)/scalef, 
   #            label = "Projections",
   #            size = 3,
   #            color = "darkgrey") +
   # )
  
}

if(plotNlims){
  p1 <- p1 + geom_line(data=dfNlims,
                       aes(x=Year,y=N30/scalef),
                       color="lightgrey",
                       size = 1.2,
                       alpha = 0.5) +
    geom_line(data=dfNlims,
              aes(x=Year,y=N50/scalef),
              color="lightgrey",
              size = 1.2,
              alpha = 0.5) + 
    geom_line(data=dfNlims,
              aes(x=Year,y=N70/scalef),
              color="lightgrey",
              size = 1.2,
              alpha = 0.5) +
    annotate("text", x = (max(df$Year)+4), y = Nlims[3]/scalef, label = "N[70]",parse = TRUE,size = 5) +
    annotate("text", x = (max(df$Year)+4), y = Nlims[2]/scalef, label = "N[50]",parse = TRUE,size = 5) +
    annotate("text", x = (max(df$Year)+4), y = Nlims[1]/scalef, label = "N[lim]",parse = TRUE,size = 5)
}

if("N0" %in% component){
  p1 <- p1 + geom_point(data = dfPups,
                        aes(x=Year,y = pupEst/scalef),
                        size = 2,
                        color = colPupest) + 
    geom_errorbar(data = dfPups,
                  aes(x = Year,ymin = LL/scalef,ymax = UL/scalef), 
                  width=1.0,
                  size = 0.5,
                  color = colPupest)
}

p1 <- p1 + theme_classic() +
  theme(text = element_text(size=20),
        plot.margin = unit(c(1,2,1,1), "cm"),
        axis.text.y = element_text(angle = 90,margin = margin(t = 0, r = 20, b = 0, l = 30),vjust = 0.5),
        axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),vjust = 1),
        #legend.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "top") +
  ylab("Abundance (in 10K)") + 
  scale_fill_manual(values = c(theCols[3], theCols[2], theCols[1])) +
  scale_colour_manual(values = c(theCols[3], theCols[2], theCols[1])) 


p1
}


plotRes <- function(results=res,
                    data=data,
                    component=c('N0', 'N1'),
                    xLim=NA,
                    yLim=NA,
                    plot.legend=TRUE,
                    plot.Nlims=TRUE, 
                    projections=TRUE, 
                    mean.proj=FALSE,
                    width = 13, 
                    height = 7,
                    conf.int = TRUE,
                    grDev = FALSE)
{
  if(projections) {
    span <- c(1:length(results$indN0))
  } else {
    span <- c(1:nrow(data$Cdata))
  }
  
  options(scipen=999)
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  if(grDev) graphDev(width = width,height = height)
  #windows("",width = width,height = height)
  
  #require(RColorBrewer)
  theCols <- RColorBrewer::brewer.pal(max(c(3, length(component))), 'Dark2')
  
  if(length(xLim)==1) xLim <- range(results$years)
  
  if('N1' %in% component) {
    indN1back <- results$indN1[c(1:match(max(data$pupProductionData[,1]), results$years))]
    lCI <- results$rep.matrix[results$indN1,1]-(1.96*results$rep.matrix[results$indN1,2])
    uCI <- results$rep.matrix[results$indN1,1]+(1.96*results$rep.matrix[results$indN1,2])
    
    if(length(yLim)==1) yLim <- c(0, max(uCI))
    
    if(!projections) {
      if(length(xLim)==1) {
        xLim <- range(results$years[span])
      } else {
        xLim[2] <- min(c(xLim[2], max(data$Cdata[,1])))
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
      polygon(c(c(results$years[1]:(max(data$Cdata[,1])+1)), 
                rev(c(results$years[1]:(max(data$Cdata[,1])+1)))), 
              c(lCI[c(1:((max(data$Cdata[,1])+1)-results$years[1]+1))],
                rev(uCI[c(1:((max(data$Cdata[,1])+1)-results$years[1]+1))])),
              col=add.alpha(theCols[1], 0.5), border=NA)
    }
    
    #lines(c(results$years[1]:max(dat$pupProductionData[,1])),results$rep.matrix[indN1back,1],
    #      col=theCols[1], lwd=2)
    lines((data$Cdata[1,1]:(max(data$Cdata[,1])+1)), results$rep.matrix[results$indN1[1:(length(data$Cdata[,1])+1)]],col=theCols[1], lwd=2)
    
  } 
  
  if('N0' %in% component) {
    indN0back <- results$indN0[c(1:match(max(data$pupProductionData[,1]), results$years))]
    plCI <- results$rep.matrix[results$indN0,1]-(1.96*results$rep.matrix[results$indN0,2])
    puCI <- results$rep.matrix[results$indN0,1]+(1.96*results$rep.matrix[results$indN0,2])
    if(!'N1' %in% component) {
      if(length(yLim)==1) yLim <- c(0, max(puCI))
      if(!projections) {
        if(length(xLim)==1) {
          xLim <- range(results$years[span])
        } else {
          xLim[2] <- min(c(xLim[2], max(data$Cdata[,1])))
        }  
        if(length(yLim)==1) {
          yLim <- c(0, max(c(max(puCI[span]), max(data$pupProductionData[,2]+(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2]))))))
        } else {
          yLim[2] <- max(c(max(puCI[span]), max(data$pupProductionData[,2]+(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])))))
        }  
      }
      plot(results$years, results$rep.matrix[results$indN0,1],
           type = "l",xlab = "Year",ylab = "Abundance",
           col=theCols[2], lwd=2, lty=2, xlim=xLim, ylim=yLim, axes=FALSE,cex.lab = 1.5)
      if(!projections | !mean.proj) lines(results$years, results$rep.matrix[results$indN0,1], lwd=2, lty=2, col=theCols[2])
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
      polygon(c(c(results$years[1]:(max(data$Cdata[,1])+0)), 
                rev(c(results$years[1]:(max(data$Cdata[,1])+0)))), 
              c(plCI[c(1:((max(data$Cdata[,1])+0)-results$years[1]+1))],
                rev(puCI[c(1:((max(data$Cdata[,1])+0)-results$years[1]+1))])),
              col=add.alpha(theCols[2], 0.5), border=NA)
    }
    
    #lines(c(results$years[1]:max(dat$pupProductionData[,1])), results$rep.matrix[indN0back,1],
    #      col=theCols[2], lwd=2)
    lines((data$Cdata[1,1]:(max(data$Cdata[,1])+0)), results$rep.matrix[results$indN0[1:(length(data$Cdata[,1])+0)]],col=theCols[2], lwd=2)
    
  } 
  
  if(length(data)>1 & 'N0' %in% component) {
    segments(data$pupProductionData[,1], 
             data$pupProductionData[,2]-(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])),
             data$pupProductionData[,1],
             data$pupProductionData[,2]+(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])),
             col=1)
    segments(data$pupProductionData[,1]-0.5, 
             data$pupProductionData[,2]-(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])),
             data$pupProductionData[,1]+0.5,
             data$pupProductionData[,2]-(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])),
             col=1)
    segments(data$pupProductionData[,1]-0.5, 
             data$pupProductionData[,2]+(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])),
             data$pupProductionData[,1]+0.5,
             data$pupProductionData[,2]+(1.96*(data$pupProductionData[,3]*data$pupProductionData[,2])),
             col=1)
    
    points(data$pupProductionData[,1],data$pupProductionData[,2],
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
        
        if(length(data)>1) {
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
      if(length(data)>1) {
        legend('topright', lwd=1, lty=NA, 
               pt.bg=theCols[2], pch=21, 
               pt.cex=1.5, 'Pup production',
               bty='n')
      }
    }  
  }
  options(scipen=0)
}
