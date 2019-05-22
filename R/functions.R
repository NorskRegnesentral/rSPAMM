ilogit1 <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(exp(x)/(1+exp(x)))
}


buildPandF <- function(Fdat = fecundity,cdata = catch_data,Fproj = Fproj,Pdat = Pdat,Pper = Pper)
{
  ##Prepare fecundity data
  #Fdat <- read.table("wgharp.fdt",header = FALSE,sep = "")

  #cdata <- scan("wgharp.cat",what="raw")
  yr1 = cdata[1,1]              #First year of catch data
  yr2 = cdata[dim(cdata)[1],1]  #Last year of catch data
  nyr = yr2 - yr1 + 1           #Number of years with catch data

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

  #write.table(Fvec,"wgharp.f",row.names = FALSE, col.names = FALSE)

  ##########################
  #Prepare birth ogive data
  #Pdat <- read.table("wgharp.ogi",sep = "",header = TRUE)
  #Pper <- read.table("wgharp.ogp",sep = "",header = TRUE)
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

  #write.table(P,"wgharp.pma",row.names = FALSE, col.names = FALSE)
 return(list(Fdt = Fvec,Pmatrix = P))
}
