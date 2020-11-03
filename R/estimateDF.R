#' \code{makeDF} Create data frame from model output
#' @param res Model output object
#' @param data Model input data object
#' @param include.catch Logical. Should catch data be included in data frame.
#' @param include.fec Logical. Should fecundity data be included in data frame.
#' @return Data frame of modelled population size per age class and (optionally) catch and fecundity data
#' @keywords output, data
#' @export
#' @examples
#' df <- makeDF(res, data)
#' 

makeDF <- function(res, data, include.catch=T, include.fec=T) {
  repnames <- dimnames(res$rep.matrix)[[1]]
  n0 <- res$rep.matrix[which(repnames=='N0'),]
  n1 <- res$rep.matrix[which(repnames=='N1'),]
  ntot <- res$rep.matrix[which(repnames=='NTot'),]
  
  df <- data.frame(Year=res$years,muN0=n0[-1,1], sigN0=n0[-1,2],
                   muN1=n1[-1,1], sigN1=n1[-1,2],
                   muNTot=ntot[-1,1], sigNTot=ntot[-1,2])
  
  df$lwrN0 <- df$muN0-(1.96 * df$sigN0)
  df$uprN0 <- df$muN0+(1.96 * df$sigN0)
  df$lwrN1 <- df$muN1-(1.96 * df$sigN1)
  df$uprN1 <- df$muN1+(1.96 * df$sigN1)
  df$lwrNTot <- df$muNTot-(1.96 * df$sigNTot)
  df$uprNTot <- df$muNTot+(1.96 * df$sigNTot)
  df$projections <- rep(T, nrow(df))
  df$projections[match(data$Cdata[,1], df$Year)] <- F
  
  pupProd <- as.data.frame(data$pupProductionData)
  names(pupProd) <- c('Year', 'mu', 'cv')
  pupProd$sd <- pupProd$cv * pupProd$mu
  pupProd$lwr <- pupProd$mu - (1.96*pupProd$sd)
  pupProd$upr <- pupProd$mu + (1.96*pupProd$sd)
  
  df$muCount <-rep(NA, nrow(df))
  df$muCount[match(pupProd$Year, df$Year)] <- pupProd$mu
  
  df$sdCount <-rep(NA, nrow(df))
  df$sdCount[match(pupProd$Year, df$Year)] <- pupProd$sd
  
  df$lwrCount <-rep(NA, nrow(df))
  df$lwrCount[match(pupProd$Year, df$Year)] <- pupProd$lwr
  
  df$uprCount <-rep(NA, nrow(df))
  df$uprCount[match(pupProd$Year, df$Year)] <- pupProd$upr
  
  if(include.catch) {
    Cdata <- as.data.frame(data$Cdata)
    names(Cdata) <- c('Year', 'C0', 'C1')
    df <- merge(df, Cdata, by='Year', all.x=T)
  }
  
  if(include.fec) {
    Fdata <- as.data.frame(data$fecundity)
    names(Fdata) <- c('Year', 'mu', 'sd')
    df <- merge(df, Fdata, by='Year', all.x=T)
    names(df)[length(df)-c(1,0)] <- c('muF', 'sdF')
  }
  
  df
}



#' \code{plotN} Plot modelled abundance and pup production estimates from data frame object
#' @param df Output data frame from \code{makeDF}
#' @param component Population components to include (N0, N1 and NTot)
#' @param plotNlims Logical. Should N-limits be plotted (currently not implemented)
#' @param plotProjections Logical. Should projections be included in trajectory.
#' @param plotMean Logical. Should mean projection be included in trajectory.
#' @param width Plot width (if \code{grDev} is set to TRUE)
#' @param height  Plot height (if \code{grDev} is set to TRUE)
#' @param grDev Set up gr device
#' @param language Set language for labels. Currently only Norwegian ('NO') and English ('EN') implemented.
#' @param brewer.palette Select which palette to use from the \code{RColorBrewer} package.
#' @return ggplot of selected population components
#' @keywords output, plotting
#' @export
#' @examples
#' plotN(df)
#' plotN(df, component=c('N0', 'NTot'), language='EN')
#' 

plotN <- function(df, 
                  component = c("N0","N1", "NTot"), 
                  plotNlims=TRUE, 
                  plotProjections=TRUE,
                  plotProjMean = TRUE,
                  width = 15, 
                  height = 10,
                  grDev = FALSE,
                  language='NO',
                  brewer.palette='Dark2') {
  
  require(ggplot2)
  require(tidyr)

  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Set extended projection back to last yea of catch data 
  ## (to avoid gap at the start of in projection lines)
  
  df$ext.proj <- df$projections
  df$ext.proj[which(df$projections)[1]-1] <- T

  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Get variable numbers and pivot data to long format:
  
  mus <- match(paste0('mu', component), names(df))
  lwrs <- match(paste0('lwr', component), names(df))
  uprs <- match(paste0('upr', component), names(df))
  columns <- c(1, match(c('projections', 'ext.proj'), names(df)))
  
  muDF <- df[,c(columns, mus)] %>% pivot_longer(cols=starts_with('mu'),
                                 names_to='group',
                                 values_to = 'mu')
  
  muDF$group <- gsub('mu', '', muDF$group)
  
  lwrDF <- df[,c(columns, lwrs)] %>% pivot_longer(cols=starts_with('lwr'),
                                          names_to='group',
                                          values_to = 'lwr')
  lwrDF$group <- gsub('lwr', '', lwrDF$group)
  
  uprDF <- df[,c(columns, uprs)] %>% pivot_longer(cols=starts_with('upr'),
                                            names_to='group',
                                            values_to = 'upr')
  uprDF$group <- gsub('upr', '', uprDF$group)
  
  DF <- merge(muDF, lwrDF, by=c('Year', 'group', 'projections', 'ext.proj'))
  DF <- merge(DF, uprDF, by=c('Year', 'group', 'projections', 'ext.proj'))

  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Get axis limits
  
  if(plotProjections) {
    projDF <- DF[which(DF$ext.proj),]
    xLims <- range(DF$Year)
    yLims <- c(0, max(DF$upr))
  }
  
  DF <- DF[which(!DF$projections),]
  
  if(!plotProjections) {
    xLims <- range(DF$Year)
    yLims <- c(0, max(DF$upr))
  }
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Set label and legend language:    
  
  if(language=='NO') {
    labs <- c('År', 'Antall')
    legs <- component
    legs <- gsub('N0', 'Unger (0-gruppe)', legs)
    legs <- gsub('N1', 'Voksne (+1-gruppe)', legs)
    legs <- gsub('NTot', 'Total', legs)
    leg.title <- 'Aldersklasse'
  } else {
    labs <- c('Year', 'Abundance')
    legs <- component
    legs <- gsub('N0', 'Pups  (0 group)', legs)
    legs <- gsub('N1', 'Adults (1+ group)', legs)
    legs <- gsub('NTot', 'Total', legs)
    leg.title <- 'Age class'
  }
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Set up graphics:    
  
  if(grDev) graphDev(width = width,height = height)
  
  y.format <- function(y) format(y, big.mark = " ")
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Make ggplot:    
  
  pt <- ggplot() +
    geom_line(data = DF,
              aes(x=Year,y=mu,
                  group = group,
                  color = group),
                  size = 1.3,
                  linetype = 1) +
    geom_ribbon(data=DF, show.legend=F, aes(x = Year, ymin=lwr,ymax=upr,
                               group = group,
                               color = NULL,
                               fill = group),
                  alpha=0.3) +
    theme_bw() +
    labs(x=labs[1], y=labs[2]) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = y.format, limits=yLims) +
    scale_x_continuous(breaks=breaks_extended(10), limits=xLims) +
    scale_color_brewer(palette=brewer.palette) + 
    scale_fill_brewer(palette=brewer.palette) +
    scale_color_discrete(name = leg.title, labels = legs)
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Add projections:
  
  if(plotProjections) {
  pt <- pt +  
    geom_line(data = projDF,
              aes(x=Year,y=lwr,
                  group = group,
                  color = group),
              size = 0.5,
              linetype = 3) +
    geom_line(data = projDF,
              aes(x=Year,y=upr,
                  group = group,
                  color = group),
              size = 0.5,
              linetype = 3)
  }
  
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Add pup count estimates:
  
  cnt <- df[which(!is.na(df$muCount)), 
            match(c('Year', 'muCount', 'lwrCount', 'uprCount'), names(df))]
  
  names(cnt) <- gsub('Count', '', names(cnt))
  
  pt <- pt +
    geom_segment(data=cnt, aes(x=Year, y=lwr, xend=Year, yend=upr),
                 color=RColorBrewer::brewer.pal(9, 'Reds')[9]) +
    geom_point(data=cnt, aes(x=Year, y=mu),
               color=RColorBrewer::brewer.pal(9, 'Reds')[9])
  
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Return plot:
  
  pt
  
} 
