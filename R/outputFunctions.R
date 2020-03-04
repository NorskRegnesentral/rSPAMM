#' Plot estimated population sizes from optimized population model for harp seals and hooded seals
#'
#' Plot population trajectories.
#' @param results Output from fitted model
#' @param data Original data on estimated pup production (set to NA if these are not to be included).
#' @param component Which population component to plot. Can be either 'N0' for pups, 'N1' for adults, or a vector with both.
#' @param xLim Manually set the x axis extent, overrides the default which is the extent of the plotted data.
#' @param yLim Manually set the y axis extent, overrides the default which is the extent of the plotted data.
#' @param plot.legend Add legend to plot (TRUE/FALSE)
#' @param plot.Nlims True/False
#' @param projections Plot projections (TRUE/FALSE)
#' @param mean.proj True/False
#' @param width Width of the figure
#' @param height Height of the figure
#' @param conf.int Logical parameter to decide wether to plot 95 percent confidence interval or not
#' @param grDev Logical parameter to decide wether to open a OS independent graphical window
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plotRes(res, data)

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

#' Plot the model results and the reported catch data
#'
#' Plot the model results and the reported catch data in the same figure.
#' @param results Output from fitted model
#' @param data Original data on estimated pup production (set to NA if these are not to be included).
#' @param component Which population component to plot. Can be either 'N0' for pups, 'N1' for adults, or a vector with both.
#' @param xLim Manually set the x axis extent, overrides the default which is the extent of the plotted data.
#' @param yLim Manually set the y axis extent, overrides the default which is the extent of the plotted data.
#' @param plot.legend Add legend to plot (TRUE/FALSE)
#' @param plot.Nlims True/False
#' @param projections Plot projections (TRUE/FALSE)
#' @param mean.proj True/False
#' @param width Width of the figure
#' @param height Height of the figure
#' @param conf.int Logical parameter to decide wether to plot 95 percent confidence interval or not
#' @param grDev Logical parameter to decide wether to open a OS independent graphical window
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plotCResatch(data$Cdata)

plotResCatch <- function (results = res, 
                          data = data, 
                          component = c("N0","N1"), 
                          xLim = NA, 
                          yLim = NA, 
                          plot.legend = TRUE, 
                          plot.Nlims = TRUE, 
                          projections = TRUE, 
                          mean.proj = FALSE, 
                          width = 13, 
                          height = 13, 
                          conf.int = TRUE, 
                          grDev = FALSE, 
                          labels='English') 
{
  if (projections) {
    span <- c(1:length(results$indN0))
  } else {
    span <- c(1:nrow(data$Cdata))
  }
  
  options(scipen = 999)
  add.alpha <- function(col, alpha = 1) {
    if (missing(col)) 
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], 
                                                       x[2], x[3], alpha = alpha))
  }
  if (grDev) graphDev(width = width, height = height)
  theCols <- RColorBrewer::brewer.pal(max(c(3, length(component))), 
                                      "Dark2")
  if (length(xLim) == 1) 
    xLim <- range(results$years)
  if ("N1" %in% component) {
    indN1back <- results$indN1[c(1:match(max(data$pupProductionData[, 
                                                                    1]), results$years))]
    lCI <- results$rep.matrix[results$indN1, 1] - (1.96 * 
                                                     results$rep.matrix[results$indN1, 2])
    uCI <- results$rep.matrix[results$indN1, 1] + (1.96 * 
                                                     results$rep.matrix[results$indN1, 2])
    if (length(yLim) == 1) 
      yLim <- c(0, max(uCI))
    if (!projections) {
      if (length(xLim) == 1) {
        xLim <- range(results$years[span])
      }
      else {
        xLim[2] <- min(c(xLim[2], max(data$Cdata[, 1])))
      }
      if (length(yLim) == 1) {
        yLim <- c(0, max(uCI[span]))
      }
      else {
        yLim[2] <- max(uCI[span])
      }
    }
    layout(matrix(c(1,2), ncol=1), heights=c(0.7, 0.3))
    par(mar = c(1, 8, 1, 2) + 0.1, las=1)
    plot(results$years, results$rep.matrix[results$indN1, 
                                           1], type = "l", xlab = "Year", ylab = "", 
         col = theCols[1], lwd = 2, lty = 2, xlim = xLim, 
         ylim = yLim, axes = FALSE, cex.lab = 1.5)
    if(labels=='Norwegian') {
      mtext(side=2, line=6, 'Bestandsst?rrelse', las=0, cex=1.5)
    } else {
      mtext(side=2, line=6, 'Abundance', las=0, cex=1.5)
    }
    if (!projections | !mean.proj) 
      lines(results$years, results$rep.matrix[results$indN1, 
                                              1], lwd = 2, lty = 2, col = "white")
    axis(1, labels=F)
    axis(2, at = pretty(par("usr")[c(3, 4)]), 
         labels = format(pretty(par("usr")[c(3, 4)]), scientific = F, big.mark=' '))
    abline(h = par("usr")[3])
    abline(v = par("usr")[1])
    if (plot.Nlims) {
      Nlims <- c(0.3, 0.5, 0.7) * max(results$rep.matrix[results$indNTot, 
                                                         1])
      abline(h = Nlims, col = "lightgrey")
      text(par("usr")[2], Nlims[1], expression(N[lim]), 
           xpd = NA, adj = 0, cex = 0.9)
      text(par("usr")[2], Nlims[2], expression(N[50]), 
           xpd = NA, adj = 0, cex = 0.9)
      text(par("usr")[2], Nlims[3], expression(N[70]), 
           xpd = NA, adj = 0, cex = 0.9)
    }
    if (projections) {
      if (conf.int) {
        lines(results$years, lCI, col = theCols[1], lty = 2)
        lines(results$years, uCI, col = theCols[1], lty = 2)
      }
      lines(results$years, results$rep.matrix[results$indN1, 
                                              1], col = theCols[1], lwd = 2, lty = 2)
    }
    if (conf.int) {
      polygon(c(c(results$years[1]:(max(data$Cdata[, 1]) + 
                                      1)), rev(c(results$years[1]:(max(data$Cdata[, 
                                                                                  1]) + 1)))), c(lCI[c(1:((max(data$Cdata[, 1]) + 
                                                                                                             1) - results$years[1] + 1))], rev(uCI[c(1:((max(data$Cdata[, 
                                                                                                                                                                        1]) + 1) - results$years[1] + 1))])), col = add.alpha(theCols[1], 
                                                                                                                                                                                                                              0.5), border = NA)
    }
    lines((data$Cdata[1, 1]:(max(data$Cdata[, 1]) + 1)), 
          results$rep.matrix[results$indN1[1:(length(data$Cdata[, 
                                                                1]) + 1)]], col = theCols[1], lwd = 2)
  }
  if ("N0" %in% component) {
    indN0back <- results$indN0[c(1:match(max(data$pupProductionData[, 
                                                                    1]), results$years))]
    plCI <- results$rep.matrix[results$indN0, 1] - (1.96 * 
                                                      results$rep.matrix[results$indN0, 2])
    puCI <- results$rep.matrix[results$indN0, 1] + (1.96 * 
                                                      results$rep.matrix[results$indN0, 2])
    if (!"N1" %in% component) {
      if (length(yLim) == 1) 
        yLim <- c(0, max(puCI))
      if (!projections) {
        if (length(xLim) == 1) {
          xLim <- range(results$years[span])
        }
        else {
          xLim[2] <- min(c(xLim[2], max(data$Cdata[, 
                                                   1])))
        }
        if (length(yLim) == 1) {
          yLim <- c(0, max(c(max(puCI[span]), max(data$pupProductionData[, 
                                                                         2] + (1.96 * (data$pupProductionData[, 3] * 
                                                                                         data$pupProductionData[, 2]))))))
        }
        else {
          yLim[2] <- max(c(max(puCI[span]), max(data$pupProductionData[, 
                                                                       2] + (1.96 * (data$pupProductionData[, 3] * 
                                                                                       data$pupProductionData[, 2])))))
        }
      }
      plot(results$years, results$rep.matrix[results$indN0, 
                                             1], type = "l", xlab = "", ylab = "", 
           col = theCols[2], lwd = 2, lty = 2, xlim = xLim, 
           ylim = yLim, axes = FALSE, cex.lab = 1.5)
      if(labels=='Norwegian') {
        mtext(side=2, line=6, 'Bestandsst?rrelse', las=0, cex=1.5)
      } else {
        mtext(side=2, line=6, 'Abundance', las=0, cex=1.5)
      }
      
      if (!projections | !mean.proj) 
        lines(results$years, results$rep.matrix[results$indN0, 
                                                1], lwd = 2, lty = 2, col = theCols[2])
      axis(1)
      axis(2, at = pretty(par("usr")[c(3, 4)]), 
           labels = format(pretty(par("usr")[c(3, 4)]), scientific = FALSE, big.mark=' '))
      abline(h = par("usr")[3])
      abline(v = par("usr")[1])
    } else {
      if (projections) 
        lines(results$years, results$rep.matrix[results$indN0, 
                                                1], col = theCols[2], lwd = 2, lty = 2)
    }
    if (projections) {
      if (conf.int) {
        lines(results$years, plCI, col = theCols[2], 
              lty = 2)
        lines(results$years, puCI, col = theCols[2], 
              lty = 2)
      }
    }
    if (conf.int) {
      polygon(c(c(results$years[1]:(max(data$Cdata[, 1]) + 
                                      0)), rev(c(results$years[1]:(max(data$Cdata[, 
                                                                                  1]) + 0)))), c(plCI[c(1:((max(data$Cdata[, 1]) + 
                                                                                                              0) - results$years[1] + 1))], rev(puCI[c(1:((max(data$Cdata[, 
                                                                                                                                                                          1]) + 0) - results$years[1] + 1))])), col = add.alpha(theCols[2], 
                                                                                                                                                                                                                                0.5), border = NA)
    }
    lines((data$Cdata[1, 1]:(max(data$Cdata[, 1]) + 0)), 
          results$rep.matrix[results$indN0[1:(length(data$Cdata[, 
                                                                1]) + 0)]], col = theCols[2], lwd = 2)
  }
  if (length(data) > 1 & "N0" %in% component) {
    segments(data$pupProductionData[, 1], data$pupProductionData[, 
                                                                 2] - (1.96 * (data$pupProductionData[, 3] * data$pupProductionData[, 
                                                                                                                                    2])), data$pupProductionData[, 1], data$pupProductionData[, 
                                                                                                                                                                                              2] + (1.96 * (data$pupProductionData[, 3] * data$pupProductionData[, 
                                                                                                                                                                                                                                                                 2])), col = 1)
    segments(data$pupProductionData[, 1] - 0.5, data$pupProductionData[, 
                                                                       2] - (1.96 * (data$pupProductionData[, 3] * data$pupProductionData[, 
                                                                                                                                          2])), data$pupProductionData[, 1] + 0.5, data$pupProductionData[, 
                                                                                                                                                                                                          2] - (1.96 * (data$pupProductionData[, 3] * data$pupProductionData[, 
                                                                                                                                                                                                                                                                             2])), col = 1)
    segments(data$pupProductionData[, 1] - 0.5, data$pupProductionData[, 
                                                                       2] + (1.96 * (data$pupProductionData[, 3] * data$pupProductionData[, 
                                                                                                                                          2])), data$pupProductionData[, 1] + 0.5, data$pupProductionData[, 
                                                                                                                                                                                                          2] + (1.96 * (data$pupProductionData[, 3] * data$pupProductionData[, 
                                                                                                                                                                                                                                                                             2])), col = 1)
    points(data$pupProductionData[, 1], data$pupProductionData[, 
                                                               2], pch = 21, bg = theCols[2], cex = 1.5)
  }
  
  par(mar = c(3, 8, 1, 2) + 0.1, las=1)
  plot(par('usr')[c(1,2)], c(0, max(data$Cdata[,-1])), 
       type='n', axes=F, xlab='', ylab='', xlim=xLim, 
       cex.lab = 1.5)
  axis(1)
  axis(2, at=pretty(c(0, max(data$Cdata[,-1]))), 
       labels=format(pretty(c(0, max(data$Cdata[,-1]))), big.mark = ' '))
  abline(h=par('usr')[3])
  abline(v=par('usr')[1])
  
  points(data$Cdata[,1]-0.2, data$Cdata[,2], 
         type='h', col=theCols[2], lwd=3,
         lend=2)
  points(data$Cdata[,1]+0.2, data$Cdata[,3], 
         type='h', col=theCols[1], lwd=3,
         lend=2)
  if(labels=='Norwegian') {
    mtext(side=2, line=6, 'Fangst', las=0, cex=1.5)
    legend('topright', pch=15, col=theCols[c(2,1)], c('Unger (0 gruppe)', 'Voksne (1+ gruppe)'), bty='n')
  } else {
    mtext(side=2, line=6, 'Catch size', las=0, cex=1.5)
    legend('topright', pch=15, col=theCols[c(2,1)], c('Pups (0 group)', 'Adults (1+ group)'), bty='n')
  }
  
  options(scipen = 0)
}



#' Plot the reported catch data
#'
#' Plot the reported catch data.
#' @param catch Reported catcj data
#' @param position Position of bars: If Pup catch and 1+ catch next to each other use position = "dodge" (default). On top of each other use position = "stack"
#' @param width Figure width
#' @param height Figure height
#' @param grDev Logical parameter to decide wether to open a OS independent graphical window
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords population model
#' @export
#' @examples
#' plotCatch(data$Cdata)

plotCatch <- function(catch = cdata,
                      width = 9,
                      height = 7,
                      position = "dodge",
                      grDev = FALSE)
{

  #theCols <- RColorBrewer::brewer.pal(3, 'Dark2')
  theCols <- RColorBrewer::brewer.pal(3, "Dark2")[c(2,1,3)]
  
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
  if(grDev) graphDev(width = width,height = height)
  pl
  
}


#' Plot fecundity data
#'
#' Plot the fecundity data used in the model fit.
#' @param data Data object used in model fit
#' @param include.observations Plot the observed fecundity rates with 95 percent confidence intervals (default = TRUE)
#' @return plot Returns a plot of predicted population size for different population components
#' @keywords fecundity data
#' @export
#' @examples
#' plotFecundity(res, data)

plotFecundity <- function(dat = data,
                          include.observations = TRUE, 
                          grDev = FALSE)
{

  if(grDev) graphDev(width = width,height = height)
  par(mar=c(5,6,4,2)+0.1)
  plot(dat$Cdata[,1],dat$Ftmp,type = "l",lwd = 3, col = "royalblue", xlab='Year', ylab='Fecundity rate',ylim=c(0,1), axes=F,cex.lab=1.5)
  if(include.observations){
    #fecundity <- read.table(paste("vignettes/data/",population,"/fecundity.dat",sep = ""),header = FALSE)
    segments(dat$fecundity[,1], 
             dat$fecundity[,2]-(1.96*(dat$fecundity[,3]*dat$fecundity[,2])),
             dat$fecundity[,1],
             dat$fecundity[,2]+(1.96*(dat$fecundity[,3]*dat$fecundity[,2])),
             col="steelblue")
    segments(dat$fecundity[,1]-0.5, 
             dat$fecundity[,2]-(1.96*(dat$fecundity[,3]*dat$fecundity[,2])),
             dat$fecundity[,1]+0.5,
             dat$fecundity[,2]-(1.96*(dat$fecundity[,3]*dat$fecundity[,2])),
             col="steelblue")
    segments(dat$fecundity[,1]-0.5, 
             dat$fecundity[,2]+(1.96*(dat$fecundity[,3]*dat$fecundity[,2])),
             dat$fecundity[,1]+0.5,
             dat$fecundity[,2]+(1.96*(dat$fecundity[,3]*dat$fecundity[,2])),
             col="steelblue")
    
    points(dat$fecundity[,1],dat$fecundity[,2],
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
#' @param grDev Logical parameter to decide wether to open a OS independent graphical window
#' @keywords birth ogive
#' @export
#' @examples
#' plotOgive()

plotOgive <- function(dat=data, 
                      highlight.last=TRUE,
                      grDev = FALSE) {
  opal <- palette()
  palette(RColorBrewer::brewer.pal(8, 'Dark2'))
  
  if(grDev) graphDev(width = width,height = height)
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

# Prephare graphical device for a given OS
# @param width Width of the graphical window
# @param height Height of the grapchical window
# @return The an OS dependent graphical window
graphDev = function(width = 7,height = 5) {
  
  system = Sys.info()
  if(system['sysname']=="Windows"){
    windows(width = 7,height = 5)
  }
  
  if(system['sysname']=="Linux"){
    X11(width = 7,height = 5)
  }
  
  if(system['sysname']=="Darwin"){
    quartz("",width = 7,height = 5)
  }
}

