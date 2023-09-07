#***************************************************************************
# Code to run contour plots over various targets on brood stock as specified 
# ppns of returns to river that determine hatchery size, against 2 other 
# management levers: percent marking and percent selective 
# removal of marked fish
# Date last revised 7 July 2023
#***************************************************************************

library(akima)
library(RColorBrewer)

# Source population dynamics models
source(here::here("R", "run.lever.model.r"))
source(here::here("R", "FunDefs.r"))

scenario <- "BSppnRR_thetaHatch60"# "BSppnRR"
c<-400000#Beverton-Holt capacity parameter (Rmax)
# Wither et al. 2018 used 400000
# Bradford et al. 2023 used 100000 smolts
percent.hatch<-0# 1-percent of hatchery fish that return to natural spawning 
HR<-0.4#Harvest rate
# Withler et al. 2018 used 0.4
# Bradford et al. 2023 use ?

w<-sqrt(100)
# Withler et al. 2018 used sqrt(100) with sqrt(1000) and sqrt(40) as sensitivity analyses
# Bradford et al. 2023 used sqrt(50) with sqrt(30) as sensitivity 

h<-sqrt(0.25)
# Withler et al. 2018 used sqrt(0.25) with sqrt(0.05) and sqrt(0.5) as sensitivity analyses
# Bradford et al. 2023 used sqrt(0.25)
Theta.hatch <- 60
# Withler et al. used 80
# Bradford et al. used 60
RS<-0.8
# Withler et al. use 0.8 with 0.2 in sensitivity analysis
# Bradford et al. 2023 used 0.8
sex.ratio <- 0.5
p <- 175#Beverton-Holt productivity parameter
# Withler et al. 2018 used 175
# Bradford et al. 2023 use 138 smolts/spawner
mar.surv.const<-TRUE#Is the marine suvival assumed to be 0.02 when estimating 
# Seq, or derived from inputed marine survival?
mar.surv <- 0.02
# Withler et al. 2018 used 0.02 (0.01 and 0.05 in sensitivity analyses)
# Bradford et al. 2023 used 0.03 
mar.surv.hatch <- 0.0024
# Withler et al. 2018 used 0.0024 (0.001 and 0.005 in sensitivity analyses)
# Bradford et al. 2023 used 0.005

# Fecundity (hard coded in run.lever.model). Withler used 4900
# Bradford et al. 2023 use 4000

# Survival from capture to brood stock (hard coded in run.lever.model)
# Withler et al. 2018 used 0.8
# Bradford et al. 2023 used 0.8

# Survival from egg to release in hatchery (hard coded in run.lever.model)
# Withler et al. 2018 used 0.88
# Bradford et al . 2023 used 0.8


# First, set up empty matrices for contour plot
PNImh<-matrix(NA,100,100) 
pHOSmh<-matrix(NA,100,100) 
pNOBmh<-matrix(NA,100,100)   
RperSmh<-matrix(NA,100,100) 
Ret.natmh<-matrix(NA,100,100) 
Ret.hatmh<-matrix(NA,100,100) 
Catchmh<-matrix(NA,100,100) 
fitmh<-matrix(NA,100,100) 
BSmarkmh<-matrix(NA,100,100)
HLmarkmh<-matrix(NA,100,100)
PNIsh<-matrix(NA,100,100) 
pHOSsh<-matrix(NA,100,100)   
pNOBsh<-matrix(NA,100,100) 
RperSsh<-matrix(NA,100,100) 
Ret.natsh<-matrix(NA,100,100) 
Ret.hatsh<-matrix(NA,100,100) 
Catchsh<-matrix(NA,100,100) 
fitsh<-matrix(NA,100,100) 
BSmarksh<-matrix(NA,100,100)
HLmarksh<-matrix(NA,100,100)
PNIms<-matrix(NA,100,100) 
pHOSms<-matrix(NA,100,100) 
pNOBms<-matrix(NA,100,100)  
RperSms<-matrix(NA,100,100) 
Ret.natms<-matrix(NA,100,100) 
Ret.hatms<-matrix(NA,100,100) 
Catchms<-matrix(NA,100,100) 
fitms<-matrix(NA,100,100) 
BSmarkms<-matrix(NA,100,100)
HLmarkms<-matrix(NA,100,100)

#Repeat selective harvest vs hatchery size relationship for 100% marking case
PNIsh100<-matrix(NA,100,100) 
pHOSsh100<-matrix(NA,100,100) 
pNOBsh100<-matrix(NA,100,100)  
RperSsh100<-matrix(NA,100,100) 
Ret.natsh100<-matrix(NA,100,100) 
Ret.hatsh100<-matrix(NA,100,100) 
Catchsh100<-matrix(NA,100,100) 
fitsh100<-matrix(NA,100,100) 
BSmarksh100<-matrix(NA,100,100)
HLmarksh100<-matrix(NA,100,100)

mh.hatchery.size<-NA 
mh.sel<-NA 
mh.BS<-NA
sh.hatchery.size<-NA 
sh.per.mark<-NA 
sh.BS<-NA
ms.hatchery.size<-NA 
ms.BS<-NA
sh100.hatchery.size<-NA 
sh100.per.mark<-NA 
sh100.BS<-NA
  
for(i in 1:100){
  print(i)
  for (j in 1:100){
    mh<-run.lever.model (per.mark = i*0.01, 
                         hatchery.size = 0.15, # Not used as BS.ppnRR=TRUE
                         sel = 0, 
                         Theta.hatch = Theta.hatch, # changed from 80 as per 
                         # Bradford et al. 2023 
                         c = c,  
                         percent.hatch = percent.hatch, 
                         HR = HR, 
                         h = h, 
                         w = w, 
                         mar.surv = mar.surv, 
                         RS = RS, 
                         mar.surv.hatch = mar.surv.hatch, 
                         sex.ratio = sex.ratio, 
                         ppn.RR = j*0.005,
                         BS.ppnRR = TRUE)
    mh.sel<-mh$sel
    if ( i == 1){ 
      mh.hatchery.size[j] <- mh$hatchery.size
      mh.BS[j]<-mh$BS[100]
      }
    PNImh[i,j]<-mh$PNI[100]
    pHOSmh[i,j]<-mh$pHOSeff[100]
    pNOBmh[i,j]<-mh$pNOB[100]
    RperSmh[i,j]<-mh$RperS[100]
    Ret.natmh[i,j]<-mh$ret.nat.preharvest[100]
    Ret.hatmh[i,j]<-mh$ret.hatch.preharvest[100]
    Catchmh[i,j]<-mh$catch[100]
    fitmh[i,j]<-mh$fit.adult[100]
    BSmarkmh[i,j]<-mh$BS.mark
    HLmarkmh[i,j]<-mh$handling.limit.mark
    
    sh <- run.lever.model (per.mark = 0.5, 
                           hatchery.size = 0.15,# Not used as BS.ppnRR=TRUE 
                           sel = i*0.01, 
                           Theta.hatch = Theta.hatch, # changed from 80 as per 
                           # Bradford et al. 2023 
                           c = c, 
                           percent.hatch = percent.hatch, 
                           HR = HR, 
                           h = h, 
                           w = w, 
                           mar.surv = mar.surv, 
                           RS = RS, 
                           mar.surv.hatch = mar.surv.hatch, 
                           sex.ratio = sex.ratio, 
                           ppn.RR = j*0.005,
                           BS.ppnRR = TRUE)
    
    sh.per.mark < -sh$per.mark
    if (i==1){
      sh.hatchery.size[j] <- sh$hatchery.size
      sh.BS[j] <- sh$BS[100]
      }
    PNIsh[i,j]<-sh$PNI[100]
    pHOSsh[i,j]<-sh$pHOSeff[100]
    pNOBsh[i,j]<-sh$pNOB[100]
    RperSsh[i,j]<-sh$RperS[100]  
    Ret.natsh[i,j]<-sh$ret.nat.preharvest[100]  
    Ret.hatsh[i,j]<-sh$ret.hatch.preharvest[100]  
    Catchsh[i,j]<-sh$catch[100]  
    fitsh[i,j]<-sh$fit.adult[100]
    BSmarksh[i,j]<-sh$BS.mark
    HLmarksh[i,j]<-sh$handling.limit.mark
    
      
    ms <- run.lever.model (per.mark = i*0.01, 
                           hatchery.size = 0.15, # Not used as BS.ppnRR=TRUE 
                           sel = j*0.01, 
                           Theta.hatch = Theta.hatch, # changed from 80 as per  
                           # Bradford et al. 2023 
                           c = c, 
                           percent.hatch = percent.hatch, 
                           HR = HR, 
                           h = h, 
                           w = w, 
                           mar.surv = mar.surv, 
                           RS = RS, 
                           mar.surv.hatch = mar.surv.hatch, 
                           sex.ratio = sex.ratio, 
                           ppn.RR = 0.30,
                           BS.ppnRR = TRUE)
      
    ms.hatchery.size <- ms$hatchery.size
    ms.BS<ms$BS[100]
    PNIms[i,j]<-ms$PNI[100]
    pHOSms[i,j]<-ms$pHOSeff[100]
    pNOBms[i,j]<-ms$pNOB[100]
    RperSms[i,j]<-ms$RperS[100]
    Ret.natms[i,j]<-ms$ret.nat.preharvest[100]
    Ret.hatms[i,j]<-ms$ret.hatch.preharvest[100]
    Catchms[i,j]<-ms$catch[100]
    fitms[i,j]<-ms$fit.adult[100]
    BSmarkms[i,j]<-ms$BS.mark
    HLmarkms[i,j]<-ms$handling.limit.mark
    
      
    #Repeat selective harvest vs hatchery size relationship for 100% marking case
    sh100 <- run.lever.model (per.mark = 1.0, 
                              hatchery.size = 0.15, # Not used as BS.ppnRR=TRUE 
                              sel = i * 0.01, 
                              Theta.hatch = Theta.hatch, # changed from 80 as 
                              # per Bradford et al. 2023 
                              c = c, 
                              percent.hatch = percent.hatch, 
                              HR = HR, 
                              h = h, 
                              w = w, 
                              mar.surv = mar.surv,
                              RS = RS, 
                              mar.surv.hatch = mar.surv.hatch, 
                              sex.ratio = sex.ratio, 
                              ppn.RR = j * 0.005,
                              BS.ppnRR = TRUE)
    
      sh100.per.mark <- sh100$per.mark
      if (i == 1){
        sh100.hatchery.size[j]<-sh100$hatchery.size
        sh100.BS[j]<-sh100$BS[100]
        }
      PNIsh100[i,j]<-sh100$PNI[100]
      pHOSsh100[i,j]<-sh100$pHOSeff[100]
      pNOBsh100[i,j]<-sh100$pNOB[100]
      RperSsh100[i,j]<-sh100$RperS[100]  
      Ret.natsh100[i,j]<-sh100$ret.nat.preharvest[100]  
      Ret.hatsh100[i,j]<-sh100$ret.hatch.preharvest[100]  
      Catchsh100[i,j]<-sh100$catch[100]  
      fitsh100[i,j]<-sh100$fit.adult[100]
      BSmarksh100[i,j]<-sh100$BS.mark
      HLmarksh100[i,j]<-sh100$handling.limit.mark
      
    }#End of  for (j in 1:100){
  }#End of for(i in 1:100){
  
  #Contour plots of management levers
  



plot.BSppnRR.contours <- function(output.label, 
                                  outputmh,  
                                  outputsh, 
                                  outputms, 
                                  outputsh100,
                                  outputNAmh,  
                                  outputNAsh, 
                                  outputNAms, 
                                  outputNAsh100){
  par(mar=c(3,4,3,2), mfrow=c(2,2), oma=c(2,1,2,0.5))
  #png("PNIonlyw4029Aug2017.png", width=4, height=4, units="in", res=1000)
  #par(mar=c(3,4,3,2), mfrow=c(1,1), oma=c(2,1,2,0.5))
  
  
  #Where are NAs?
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputNAsh100)
  InterpListHLNAsh100 <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  

  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputNAsh)
  InterpListHLNAsh <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  

  Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputNAmh)
  InterpListHLNAmh <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  
  Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Sel <- rep(1:100*0.01, each = 100)
  Data_DF$z <- as.vector(outputNAms)
  InterpListHLNAms <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=200, ny=200)
  
  # Put data in vector form
  Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputmh)
  #write.csv(Data_DF, "ContourData_DF.csv")
  
  InterpList <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  # InterpList[[3]][which(InterpListHLNAmh[[3]]>0)]<-NA
  
  # get x and y ranges for plotting
  xrange<-range(InterpList[[1]])
  yrange<-range(InterpList[[2]])
  zrange<-range(InterpList[[3]], na.rm=T)
  zlevels<-pretty(zrange,20)
  nlevels<-(length(zlevels))
  
  # make color palette
  cp<-colorRampPalette(brewer.pal(9,'Blues'))
  cp4<-rev(cp(nlevels+2)[-c(1,2)])
  if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}#(output.label=="pHOSeff") 
  cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
  
  plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
  .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange, nlevels) ,cp4) 
  if(output.label=="PNI"|output.label=="pNOB")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1) ,cp4[c(1,nlevels/2, nlevels)])} 
  if(output.label=="pHOSeff")#"pHOSeff")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=5, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.1), add=TRUE, lty=c("solid"), lwd=c(1))}
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  usr<-par("usr")
  axis(1); axis(2)
  scaledBS <- mh.BS/mh$Seq
  if(scenario == "BSppnRR") {int <- c(0,0.1, 0.25, 0.5, 1, 10)}
  if(scenario != "BSppnRR") {
    # int <- pretty(quantile(scaledBS, c(0.05,0.2,0.4,0.6,0.8, 0.95)))
    int <- signif(quantile(scaledBS, c(0.05,0.5,0.95)), 1)
    # int <- int[which(int>0.00999)]
  }
  
  axis(side=4, at=c(length(which(scaledBS<int[1]))*0.005,
                    length(which(scaledBS<int[2]))*0.005, 
                    length(which(scaledBS<int[3]))*0.005,                  
                    length(which(scaledBS<int[4]))*0.005,
                    length(which(scaledBS<int[5]))*0.005,
                    length(which(scaledBS<int[6]))*0.005),
       labels=NA, tcl=-0.3, cex.axis=0.8)
  mtext(side=4, text=int[1], line=0.3, at=length(which(scaledBS<int[1]))*0.005, cex=0.8)
  mtext(side=4, text=int[2], line=0.3, at=length(which(scaledBS<int[2]))*0.005, cex=0.8)
  mtext(side=4, text=int[3], line=0.3, at=length(which(scaledBS<int[3]))*0.005, cex=0.8)
  mtext(side=4, text=int[4], line=0.3, at=length(which(scaledBS<int[4]))*0.005, cex=0.8)
  mtext(side=4, text=int[5], line=0.3, at=length(which(scaledBS<int[5]))*0.005, cex=0.8)
  mtext(side=4, text=int[6], line=0.3, at=length(which(scaledBS<int[6]))*0.005, cex=0.8)
  if(scenario != "BSppnRR") {# Add tick at 0.3 of Seq for sensitivity analyses
    axis(side=4, at=(length(which(scaledBS<0.3))*0.005), 
         labels=NA, tcl=-0.3, cex.axis=0.8)
    mtext(side=4, text=0.3, line=0.3, at=length(which(scaledBS<0.3))*0.005, 
          cex=0.8)
    if(int[2]>0.3){ #Add tick mark at 1, if there is space
      axis(side=4, at=(length(which(scaledBS<1))*0.005), 
           labels=NA, tcl=-0.3, cex.axis=0.8)
      mtext(side=4, text=1, line=0.3, at=length(which(scaledBS<1))*0.005, cex=0.8)
    }
  }
  mtext(side=4, text="Hatchery size (ppn of Seq)", line=1, cex=0.7)
  abline(h=usr[3], v=usr[1])
  abline(v=usr[2])
  mtext(side=1, text="Proportion marked", line=2.5, cex=0.8)
  mtext(side=2, text="Brood stock(ppn of returns)\n", line=1.1, cex=0.8)
  mtext(side=3, text="(No selective removals of marked fish)", line=0.1, cex=0.6, at=1, adj=1)
  mtext(side=3, text="(a)", line=0.1, cex=0.8, at=0, adj=0)
  #mtext(3, text=paste(output.label, sep=""), line=2, outer=T, cex=1.5)
  
  # Make vector of grids for hatching/stippling where NAs are
  incl <- which(InterpListHLNAmh[[3]]>0)
  
  # Make polygons for each grid for hatching/stippling
  # See FuncDefs.r and 
  # https://stackoverflow.com/questions/11736996/adding-stippling-to-image-contour-plot
  polys <- matrix.poly(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], n=incl)
  # Hatching where NAs are
  for(i in seq(polys)){
    polygon(polys[[i]], density=25, angle=45, border=NA, col=grey(0.33), lwd=0.5)
    polygon(polys[[i]], density=25, angle=-45, border=NA, col=grey(0.33), lwd=0.5)
  }
  
  # Variables (management levers) to assess: Sel vs Hatchery size
  # Put data in vector form
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputsh)
  #write.csv(Data_DF, "ContourData_DF.csv")
  
  InterpList <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  # InterpList[[3]][which(InterpListHLNAsh[[3]]>0)]<-NA
  
  # get x and y ranges for plotting
  xrange<-range(InterpList[[1]])
  yrange<-range(InterpList[[2]])
  zrange<-range(InterpList[[3]], na.rm=T)
  zlevels<-pretty(zrange,20)
  nlevels<-(length(zlevels))
  
  # make color pallette
  cp<-colorRampPalette(brewer.pal(9,'Blues'))
  cp4<-rev(cp(nlevels+2)[-c(1,2)])
  if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}#(output.label=="pHOSeff") 
  cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
  
  plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
  .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange, nlevels) ,cp4) 
  if(output.label=="PNI"|output.label=="pNOB")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1) ,cp4[c(1,nlevels/2, nlevels)])} 
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=5, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.1), add=TRUE, lty=c("solid"), lwd=c(1))}
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  usr<-par("usr")
  axis(1); axis(2)
  
  scaledBS <- sh.BS/sh$Seq
  if(scenario == "BSppnRR") {int <- c(0,0.1, 0.25, 0.5, 1, 10)}
  if(scenario != "BSppnRR") {
    # int <- pretty(quantile(scaledBS, c(0.05,0.2,0.4,0.6,0.8, 0.95)))
    int <- signif(quantile(scaledBS, c(0.05,0.5,0.95)), 1)
    # int <- int[which(int>0.00999)]
  }
  axis(side=4, at=c(length(which(scaledBS<int[1]))*0.005, 
                    length(which(scaledBS<int[2]))*0.005, 
                    length(which(scaledBS<int[3]))*0.005, 
                    length(which(scaledBS<int[4]))*0.005,
                    length(which(scaledBS<int[5]))*0.005,
                    length(which(scaledBS<int[6]))*0.005),
       labels=NA, tcl=-0.3, cex.axis=0.8)
  mtext(side=4, text=int[1], line=0.3, at=length(which(scaledBS<int[1]))*0.005, cex=0.8)
  mtext(side=4, text=int[2], line=0.3, at=length(which(scaledBS<int[2]))*0.005, cex=0.8)
  mtext(side=4, text=int[3], line=0.3, at=length(which(scaledBS<int[3]))*0.005, cex=0.8)
  mtext(side=4, text=int[4], line=0.3, at=length(which(scaledBS<int[4]))*0.005, cex=0.8)
  mtext(side=4, text=int[5], line=0.3, at=length(which(scaledBS<int[5]))*0.005, cex=0.8)
  mtext(side=4, text=int[6], line=0.3, at=length(which(scaledBS<int[6]))*0.005, cex=0.8)
  if(scenario != "BSppnRR") {# Add tick at 0.3 of Seq for sensitivity analyses
    axis(side=4, at=(length(which(scaledBS<0.3))*0.005), 
         labels=NA, tcl=-0.3, cex.axis=0.8)
    mtext(side=4, text=0.3, line=0.3, at=length(which(scaledBS<0.3))*0.005, cex=0.8)
    if(int[2]>0.3){ #Add tick mark at 1, if there is space
      axis(side=4, at=(length(which(scaledBS<1))*0.005), 
           labels=NA, tcl=-0.3, cex.axis=0.8)
      mtext(side=4, text=1, line=0.3, at=length(which(scaledBS<1))*0.005, cex=0.8)
    }
  }
  mtext(side=4, text="Hatchery size (ppn of Seq)", line=1, cex=0.7)
  abline(h=usr[3], v=usr[1])
  abline(v=usr[2])
  mtext(side=1, text="Proportion of marked fish\nselectively removed", line=2.8, cex=0.8)
  mtext(side=2, text="Brood stock(ppn of returns)\n", line=1.1, cex=0.8)
  mtext(side=3, text="(50% marking)", line=0.1, cex=0.6, at=1, adj=1)
  mtext(side=3, text="(b)", line=0.1, cex=0.8, at=0, adj=0)
  # if(output.label=="pHOSeff") mtext(side=3, text=expression(paste("pHOSeff: Low ", gamma)), line=2, outer=F, cex=1.5, at=-0.2)#(output.label=="pHOSeff") #text=expression('pHOS'[eff])
  # if(output.label=="PNI") mtext(side=3, text=expression(paste("PNI: Low ", gamma)), line=2, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="pNOB") mtext(side=3, text=expression(paste("pNOB: Low ", gamma)), line=2, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="Recruits from river spawners (HOS+NOS)") mtext(side=3, text=expression(paste("Recruits from river spawners (HOS+NOS):Low ", gamma)), line=1.5, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="Recruits from hatchery production") mtext(side=3, text=expression(paste("Recruits from hatchery production:Low ", gamma)), line=1.5, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="Catch") mtext(side=3, text=expression(paste("Catch: Low ", gamma)), line=2, outer=F, cex=1.5, at=-0.2)
  if(RS!=0.8){
    if(output.label=="pHOSeff") mtext(side=3, text=expression(paste("pHOSeff: low", gamma)), line=2, outer=F, cex=1.5, at=-0.2)#(output.label=="pHOSeff") #text=expression('pHOS'[eff])
    if(output.label=="PNI") mtext(side=3, text=expression(paste("PNI: low", gamma)), line=2, outer=F, cex=1.5, at=-0.2)
    if(output.label=="pNOB") mtext(side=3, text=expression(paste("pNOB: low", gamma)), line=2, outer=F, cex=1.5, at=-0.2)
    if(output.label=="Recruits from river spawners (HOS+NOS)") mtext(side=3, text=expression(paste("Recruits from river spawners (HOS+NOS): low", gamma)), line=1.5, outer=F, cex=1.5, at=-0.2)
    if(output.label=="Recruits from hatchery production") mtext(side=3, text=expression(paste("Recruits from hatchery production: low", gamma)), line=1.5, outer=F, cex=1.5, at=-0.2)
    if(output.label=="Catch") mtext(side=3, text=expression(paste("Catch: low", gamma)), line=2, outer=F, cex=1.5, at=-0.2)
    
  }
  if(RS==0.8){
    if(output.label=="pHOSeff") mtext(side=3, text=paste("pHOSeff", title, sep=""), line=1.5, outer=F, cex=1.3, at=-0.2)#(output.label=="pHOSeff") #text=expression('pHOS'[eff])
    if(output.label=="PNI") mtext(side=3, text=paste("PNI", title, sep=""), line=1.5, outer=F, cex=1.3, at=-0.2)
    if(output.label=="pNOB") mtext(side=3, text=paste("pNOB", title, sep=""), line=1.5, outer=F, cex=1.3, at=-0.2)
    if(output.label=="Recruits from river spawners (HOS+NOS)") mtext(side=3, text=paste("Recruits from river spawners (HOS+NOS)", title, sep=""), line=1.5, outer=F, cex=1.3, at=-0.2)
    if(output.label=="Recruits from hatchery production") mtext(side=3, text=paste("Recruits from hatchery production", title, sep=""), line=1.5, outer=F, cex=1.3, at=-0.2)
    if(output.label=="Catch") mtext(side=3, text=paste("Catch", title, sep=""), line=1.5, outer=F, cex=1.3, at=-0.2)
    
  }
  # if(output.label=="pHOSeff") mtext(side=3, text=expression(paste("pHOSeff\nWeak selection")), line=2, outer=F, cex=1.5, at=-0.2)#(output.label=="pHOSeff") #text=expression('pHOS'[eff])
  # if(output.label=="PNI") mtext(side=3, text=expression(paste("PNI\nWeak selection")), line=2, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="pNOB") mtext(side=3, text=expression(paste("pNOB\nWeak selection")), line=2, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="Recruits from river spawners (HOS+NOS)") mtext(side=3, text=expression(paste("Recruits from river spawners (HOS+NOS)\nWeak selection")), line=1.5, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="Recruits from hatchery production") mtext(side=3, text=expression(paste("Recruits from hatchery production\nWeak selection")), line=1.5, outer=F, cex=1.5, at=-0.2)
  # if(output.label=="Catch") mtext(side=3, text=expression(paste("Catch\nWeak selection")), line=2, outer=F, cex=1.5, at=-0.2)
  #if(output.label=="Catch") mtext(side=3, text=paste(output.label, ": Low cf", sep=""), line=2, outer=F, cex=1.5, at=-0.2)
  # Make vector of grids for hatching/stippling where NAs are
  incl <- which(InterpListHLNAsh[[3]]>0)
  
  # Make polygons for each grid for hatching/stippling
  # See FuncDefs.r and 
  # https://stackoverflow.com/questions/11736996/adding-stippling-to-image-contour-plot
  polys <- matrix.poly(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], n=incl)
  # Hatching where NAs are
  for(i in seq(polys)){
    polygon(polys[[i]], density=25, angle=45, border=NA, col=grey(0.33), lwd=0.5)
    polygon(polys[[i]], density=25, angle=-45, border=NA, col=grey(0.33), lwd=0.5)
  }
  
  
  # Variables (management levers) to assess: %Marking vs selectivity
  # Put data in vector form
  Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Sel <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(outputms)
  
  # InterpList <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=15, ny=15)
  InterpList <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=200, ny=200)
  # InterpList[[3]][which(InterpListHLNAms[[3]]>0)]<-NA
  
  # If all NAs, then skip, where 200 x 200 is 
  # interpolation matrix size
  if(sum(is.na(InterpList[[3]])) == (200*200)){
    plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", 
         ylab="")
  }
  
  if(sum(is.na(InterpList[[3]])) != (200*200)){ #If not all NAs, then plot
    # get x and y ranges for plotting
    xrange<-range(InterpList[[1]])
    yrange<-range(InterpList[[2]])
    zrange<-range(InterpList[[3]], na.rm=T)
    zlevels<-pretty(zrange,20)
    nlevels<-(length(zlevels))
    
    # make color pallette
    cp<-colorRampPalette(brewer.pal(9,'Blues'))
    cp4<-rev(cp(nlevels+2)[-c(1,2)])
    if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}#(output.label=="pHOSeff") 
    cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
    cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
    
    plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
    .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange, nlevels) ,cp4) 
    if(output.label=="PNI"|output.label=="pNOB")
    {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1), cp4[c(1,nlevels/2, nlevels)])} 
    if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
    {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
    contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=6, add=TRUE)
    if(output.label=="PNI"|output.label=="pNOB")
    {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
    if(output.label=="PNI"|output.label=="pNOB")
    {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.1), add=TRUE, lty=c("solid"), lwd=c(1))}
    if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
    {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
    usr<-par("usr")
    axis(1); axis(2)
    abline(h=usr[3], v=usr[1])
    mtext(side=1, text="Proportion marked", line=2.5, cex=0.8)
    mtext(side=2, text="Proportion marked fish\nselectively removed", line=2, cex=0.8)
    mtext(side=3, text="(BS = 0.3 of returns)", line=0.1, cex=0.6, at=1, adj=1)
    mtext(side=3, text="(c)", line=0.1, cex=0.8, at=0, adj=0)
    #mtext(side=3, text=paste(output.label, sep=""), line=2, outer=F, cex=1.5, at=0.5)
    if(output.label=="Recruits from hatchery production"){text(x=0.1, y=0.95, labels=round(Ret.hatms[1],0), cex=0.8)}
  }
  
  # Make vector of grids for hatching/stippling where NAs are
  incl <- which(InterpListHLNAms[[3]]>0)
  if(length(incl)>0){
    # Make polygons for each grid for hatching/stippling
    # See FuncDefs.r and 
    # https://stackoverflow.com/questions/11736996/adding-stippling-to-image-contour-plot
    polys <- matrix.poly(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], n=incl)
    # Hatching where NAs are
    for(i in seq(polys)){
      polygon(polys[[i]], density=25, angle=45, border=NA, col=grey(0.33), lwd=0.5)
      polygon(polys[[i]], density=25, angle=-45, border=NA, col=grey(0.33), lwd=0.5)
    }
  }
    
    
  
  # Variables (management levers) to assess: Sel vs Hatchery size for 100% marking scenario
  # Put data in vector form
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputsh100)
  #write.csv(Data_DF, "ContourData_DF.csv")
  
  InterpList <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  # InterpList[[3]][which(InterpListHLNAsh100[[3]]>0)]<-NA
  
  
  # get x and y ranges for plotting
  xrange<-range(InterpList[[1]])
  yrange<-range(InterpList[[2]])
  zrange<-range(InterpList[[3]], na.rm=T)
  zlevels<-pretty(zrange,20)
  nlevels<-(length(zlevels))
  
  # make color pallette
  cp<-colorRampPalette(brewer.pal(9,'Blues'))
  cp4<-rev(cp(nlevels+2)[-c(1,2)])
  if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}#(output.label=="pHOSeff") 
  cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
  
  plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
  .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange, nlevels) ,cp4) 
  if(output.label=="PNI"|output.label=="pNOB")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1) ,cp4[c(1,nlevels/2, nlevels)])} 
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=5, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.1), add=TRUE, lty=c("solid"), lwd=c(1))}
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  usr<-par("usr")
  axis(1); axis(2)
  
  scaledBS <- sh100.BS/sh100$Seq
  if(scenario == "BSppnRR") {int <- c(0,0.1, 0.25, 0.5, 1, 10)}
  if(scenario != "BSppnRR") {
    # int <- pretty(quantile(scaledBS, c(0.05,0.2,0.4,0.6,0.8, 0.95)))
    int <- signif(quantile(scaledBS, c(0.05,0.5,0.95)), 1)
    # int <- int[which(int>0.00999)]
  }
  axis(side=4, at=c(length(which(scaledBS<int[1]))*0.005, 
                    length(which(scaledBS<int[2]))*0.005, 
                    length(which(scaledBS<int[3]))*0.005, 
                    length(which(scaledBS<int[4]))*0.005,
                    length(which(scaledBS<int[5]))*0.005,
                    length(which(scaledBS<int[6]))*0.005),
       labels=NA, tcl=-0.3, cex.axis=0.8)
  
  mtext(side=4, text=int[1], line=0.3, at=length(which(scaledBS<int[1]))*0.005, cex=0.8)
  mtext(side=4, text=int[2], line=0.3, at=length(which(scaledBS<int[2]))*0.005, cex=0.8)
  mtext(side=4, text=int[3], line=0.3, at=length(which(scaledBS<int[3]))*0.005, cex=0.8)
  mtext(side=4, text=int[4], line=0.3, at=length(which(scaledBS<int[4]))*0.005, cex=0.8)
  mtext(side=4, text=int[5], line=0.3, at=length(which(scaledBS<int[5]))*0.005, cex=0.8)
  mtext(side=4, text=int[6], line=0.3, at=length(which(scaledBS<int[6]))*0.005, cex=0.8)
  if(scenario != "BSppnRR") {# Add tick at 0.3 of Seq for sensitivity analyses
    axis(side=4, at=(length(which(scaledBS<0.3))*0.005), 
         labels=NA, tcl=-0.3, cex.axis=0.8)
    mtext(side=4, text=0.3, line=0.3, at=length(which(scaledBS<0.3))*0.005, cex=0.8)
    if(int[2]>0.3){ #Add tick mark at 1, if there is space
      axis(side=4, at=(length(which(scaledBS<1))*0.005), 
           labels=NA, tcl=-0.3, cex.axis=0.8)
      mtext(side=4, text=1, line=0.3, at=length(which(scaledBS<1))*0.005, cex=0.8)
    }
  }
  mtext(side=4, text="Hatchery size (ppn of Seq)", line=1, cex=0.7)
  abline(h=usr[3], v=usr[1])
  abline(v=usr[2])
  mtext(side=1, text="Proportion marked fish\nselectively removed", line=2.8, cex=0.8)
  mtext(side=2, text="Brood stock(ppn of returns)\n", line=1.1, cex=0.8)
  mtext(side=3, text="(100% marking)", line=0.1, cex=0.6, at=1, adj=1)
  mtext(side=3, text="(d)", line=0.1, cex=0.8, at=0, adj=0)
  if(output.label=="pNOB"){if(pNOBsh100[1,100]==1){text(x=0.1, y=0.48, labels=c("1.0"), cex=0.7)}}#round(pNOBsh100[1],1)
  # Make vector of grids for hatching/stippling where NAs are
  incl <- which(InterpListHLNAsh100[[3]]>0)
  
  # Make polygons for each grid for hatching/stippling
  # See FuncDefs.r and 
  # https://stackoverflow.com/questions/11736996/adding-stippling-to-image-contour-plot
  polys <- matrix.poly(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], n=incl)
  # Hatching where NAs are
  for(i in seq(polys)){
    polygon(polys[[i]], density=25, angle=45, border=NA, col=grey(0.33), lwd=0.5)
    polygon(polys[[i]], density=25, angle=-45, border=NA, col=grey(0.33), lwd=0.5)
  }
  
  
}# End of function plot.BSppnRR.contours()


title <- ""
dir <- here::here("Results", "PpnReturnsRiver", "MainResults")
output.label <-  "PNI"
# outputmh<-PNImh
# outputsh<-PNIsh
# outputms<-PNIms
# outputsh100<-PNIsh100
# outputNAsh100<-HLmarksh100
# outputNAsh<-HLmarksh
# outputNAmh<-HLmarkmh
# outputNAms<-HLmarkms
# output.label<-"PNI"


png(paste(dir, "/", output.label, scenario, ".png", sep=""), width=6, height=6, 
    units="in", res=1000)

plot.BSppnRR.contours (output.label = output.label, 
                       outputmh = PNImh,  
                       outputsh = PNIsh, 
                       outputms = PNIms, 
                       outputsh100 = PNIsh100,
                       outputNAmh = HLmarkmh,  
                       outputNAsh = HLmarksh, 
                       outputNAms = HLmarkms, 
                       outputNAsh100 = HLmarksh100)
dev.off()
