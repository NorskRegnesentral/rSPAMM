#' Plot estimated population sizes from optimized population model for harp seals and hooded seals
#'
#' Plot population trajectories.
#' @param results Output from fitted model
#' @param dat Original data on estimated pup production (set to NA if these are not to be included).
#' @param component Which opulation component to plot. Can be either 'N0' for pups, 'N1' for adults, or a vector with both.
#' @param xLim Manually set the x axis extent, overrides the default which is the extent of the plotted data.
#' @param yLim Manually set the y axis extent, overrides the default which is the extent of the plotted data.
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plot.N(res, data)

plot.N <- function(results=res,dat=data,component=c('N0', 'N1'),
                   xLim=NA,yLim=NA,plot.legend=T,plot.Nlims=T, projections=T, mean.proj=F)
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
  
  require(RColorBrewer)
  theCols <- brewer.pal(max(c(3, length(component))), 'Dark2')
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
    plot(results$years, results$rep.matrix[results$indN1,1],
         type = "l",xlab = "",ylab = "Abundance",
         col=theCols[1], lwd=2, lty=2, xlim=xLim, ylim=yLim, axes=F)
    if(!projections | !mean.proj) lines(results$years, results$rep.matrix[results$indN1,1], lwd=2, lty=2, col='white')
    axis(1)
    axis(2, at=pretty(par('usr')[c(3,4)]), 
         labels=format(pretty(par('usr')[c(3,4)]), scientific=F))
    abline(h=par('usr')[3])
    abline(v=par('usr')[1])
    if(plot.Nlims) {
      Nlims <- c(0.3, 0.5, 0.7)*max(results$rep.matrix[results$indN1,1]) 
      abline(h=Nlims, col='lightgrey')
      text(par('usr')[2], Nlims[1], expression(N[lim]), xpd=NA, adj=0, cex=0.7)
      text(par('usr')[2], Nlims[2], expression(N[50]), xpd=NA, adj=0, cex=0.7)
      text(par('usr')[2], Nlims[3], expression(N[70]), xpd=NA, adj=0, cex=0.7)
    }
    
    if(projections) {
      lines(results$years, lCI, col=theCols[1], lty=2)
      lines(results$years, uCI, col=theCols[1], lty=2)
    }
    
    polygon(c(c(results$years[1]:max(dat$pupProductionData[,1])), 
              rev(c(results$years[1]:max(dat$pupProductionData[,1])))), 
            c(lCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))],
              rev(uCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))])),
            col=add.alpha(theCols[1], 0.5), border=NA)
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
           type = "l",xlab = "",ylab = "Abundance",
           col=theCols[2], lwd=2, lty=2, xlim=xLim, ylim=yLim, axes=F)
      if(!projections | !mean.proj) lines(results$years, results$rep.matrix[results$indN0,1], lwd=2, lty=2, col='white')
      axis(1)
      axis(2, at=pretty(par('usr')[c(3,4)]), 
           labels=format(pretty(par('usr')[c(3,4)]), scientific=F))
      abline(h=par('usr')[3])
      abline(v=par('usr')[1])
    } else { 
      if(projections) lines(results$years, results$rep.matrix[results$indN0,1],
            col=theCols[2], lwd=2,lty=2)
    }
    if(projections) {
      lines(results$years, plCI, col=theCols[2], lty=2)
      lines(results$years, puCI, col=theCols[2], lty=2)
    }  
    
    polygon(c(c(results$years[1]:max(dat$pupProductionData[,1])), 
              rev(c(results$years[1]:max(dat$pupProductionData[,1])))), 
            c(plCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))],
              rev(puCI[c(1:(max(dat$pupProductionData[,1])-results$years[1]+1))])),
            col=add.alpha(theCols[2], 0.5), border=NA)
    
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
               'Total population size',
               bty='n')
        }  else {
        legend('topright', lwd=rep(2, 2), col=theCols[c(1,2)], 
               c('Total population size', 
                 'Pup production'), pch=rep(NA, 2),
               bty='n')
      
        if(length(dat)>1) {
          legend('topright', lwd=rep(1, 2), lty=NA, 
                 pt.bg=c(NA, theCols[2]), pch=c(NA,21),
                 pt.cex=c(0, 1.5),
                 c('Total population size', 
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
             results$D1, results$D1New)
  sds <- c(round(results$rep.matrix[results$indN1[1],2]), 
           results$rep.matrix[3,2],
           results$rep.matrix[2,2], 
           results$N0Current.sd, 
           results$N1Current.sd, 
           results$NTotCurrent.sd, 
           results$D1.sd, results$D1New.sd)
  
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
  matplot(t(dat$Pmat), type='l', lty=1, col='grey',
          xlab='Age (years)',
          ylab='Proportion of mature females',
          ylim=c(0,1), axes=F)
  axis(1)
  axis(2, at=seq(0, 1, by=0.2))
  box()
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
           c(leglab, 'Between periods'), bty='n', cex=0.7)
  } else {
    legend('bottomright', lwd=c(rep(2, length(ymat)), 1), 
           col=c(c(1:length(ymat)), 'grey'), 
           c(leglab, 'Between periods'), bty='n', cex=0.7)
  }
  palette(opal)       
}

