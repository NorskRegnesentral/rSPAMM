#' Plot modelled population dynamics
#'
#' Plot population trajectories. This function does not work properly, use plot.N() instead.
#' @param results Output from fitted model
#' @param data Original data on estimated pup production (set to NA if these are not to be included).
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

plot.res.ggplot <- function(res=res,data=data,component=c('N0', 'N1'),
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



#CODE BELOW DOES NOT WORK ANYMORE. NEW VERSIONS OF FUNCTIONS
#ARE CREATED
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

