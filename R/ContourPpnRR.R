#***************************************************************************
# Code to run contour plots over various limits on escapement removals for 
# brood take as specified ppns of returns to river, against 3 management levers
# (a) target hartchery size, (b) percent marking, and (c) percent selective 
# removal of marked fish
# Date last revised 26 Feb 2023
#***************************************************************************


# Source population dynamics models
source(here::here("R", "run.lever.model.r"))


# Run population dynamics model over various removal limits and management 
# levers

# First, set up empty matrices for contour plot
PNIph <- matrix(NA,100,100)
pHOSph <- matrix(NA,100,100)
pNOBph <- matrix(NA,100,100)
RperSph <- matrix(NA,100,100)
Ret.natph <- matrix(NA,100,100)
Ret.hatph <- matrix(NA,100,100)
Catchph <- matrix(NA,100,100) 
fitph <- matrix(NA,100,100)
BSmarkph <- matrix(NA,100,100)

ph.hatchery.size <- NA
ph.sel <- NA
ph.BS <- NA

PNIpm <- matrix(NA,100,100)
pHOSpm <- matrix(NA,100,100)
pNOBpm <- matrix(NA,100,100)
RperSpm <- matrix(NA,100,100)
Ret.natpm <- matrix(NA,100,100)
Ret.hatpm <- matrix(NA,100,100)
Catchpm <- matrix(NA,100,100) 
fitpm <- matrix(NA,100,100)
BSmarkpm <- matrix(NA,100,100)

pm.hatchery.size <- NA
pm.sel <- NA
pm.BS <- NA

PNIps <- matrix(NA,100,100)
pHOSps <- matrix(NA,100,100)
pNOBps <- matrix(NA,100,100)
RperSps <- matrix(NA,100,100)
Ret.natps <- matrix(NA,100,100)
Ret.hatps <- matrix(NA,100,100)
Catchps <- matrix(NA,100,100) 
fitps <- matrix(NA,100,100)
BSmarkps <- matrix(NA,100,100)

ps.hatchery.size <- NA
ps.sel <- NA
ps.BS <- NA

# Run model over values of 100 values of ppn.RR and 3 management levers: 
# hatchery.size,per.mark, and sel
for(i in 1:100){
  print(i)
  for (j in 1:100){
    ph<-run.lever.model(ppn.RR = i*0.005, 
                        per.mark = 1, 
                        hatchery.size = j*0.005, 
                        sel = 0, 
                        Theta.hatch = 80, 
                        c = c, 
                        percent.hatch = percent.hatch, 
                        HR = HR, 
                        h = h, 
                        w = w, 
                        mar.surv = 0.02, 
                        RS = RS, 
                        mar.surv.hatch = 0.0024, 
                        sex.ratio = sex.ratio)
    ph.sel <- ph$sel
    
    if (i==1) { 
      ph.hatchery.size[j] <- ph$hatchery.size
      ph.BS[j]<-ph$BS[100]
      }
    
    PNIph[i,j]<-ph$PNI[100]
    pHOSph[i,j]<-ph$pHOSeff[100]
    pNOBph[i,j]<-ph$pNOB[100]
    RperSph[i,j]<-ph$RperS[100]
    Ret.natph[i,j]<-ph$ret.nat.preharvest[100]
    Ret.hatph[i,j]<-ph$ret.hatch.preharvest[100]
    Catchph[i,j]<-ph$catch[100]
    fitph[i,j]<-ph$fit.adult[100]
    BSmarkph[i,j]<-ph$BS.mark
    
    
    pm<-run.lever.model(ppn.RR = i*0.005, 
                        per.mark = j*0.01, 
                        hatchery.size = 0.15, 
                        sel = 0, 
                        Theta.hatch = 80, 
                        c = c, 
                        percent.hatch = percent.hatch, 
                        HR = HR, 
                        h = h, 
                        w = w, 
                        mar.surv = 0.02, 
                        RS = RS, 
                        mar.surv.hatch = 0.0024, 
                        sex.ratio = sex.ratio)
    pm.sel <- pm$sel
    
    if (i==1) { 
      pm.hatchery.size[j] <- pm$hatchery.size
      pm.BS[j]<-ph$BS[100]
    }
    
    PNIpm[i,j]<-pm$PNI[100]
    pHOSpm[i,j]<-pm$pHOSeff[100]
    pNOBpm[i,j]<-pm$pNOB[100]
    RperSpm[i,j]<-pm$RperS[100]
    Ret.natpm[i,j]<-pm$ret.nat.preharvest[100]
    Ret.hatpm[i,j]<-pm$ret.hatch.preharvest[100]
    Catchpm[i,j]<-pm$catch[100]
    fitpm[i,j]<-pm$fit.adult[100]
    BSmarkpm[i,j]<-pm$BS.mark
    
    ps<-run.lever.model(ppn.RR = i*0.005, 
                        per.mark = 1, 
                        hatchery.size = 0.15, 
                        sel = j*0.01, 
                        Theta.hatch = 80, 
                        c = c, 
                        percent.hatch = percent.hatch, 
                        HR = HR, 
                        h = h, 
                        w = w, 
                        mar.surv = 0.02, 
                        RS = RS, 
                        mar.surv.hatch = 0.0024, 
                        sex.ratio = sex.ratio)
    ps.per.mark<-ps$per.mark
    ps.hatchery.size <- ps$hatchery.size
    if (i==1) { 
      ps.BS[j]<-ps$BS[100]
    }
    
    PNIps[i,j] <- ps$PNI[100]
    pHOSps[i,j] <- ps$pHOSeff[100]
    pNOBps[i,j] <- ps$pNOB[100]
    RperSps[i,j] <- ps$RperS[100]
    Ret.natps[i,j] <- ps$ret.nat.preharvest[100]
    Ret.hatps[i,j] <- ps$ret.hatch.preharvest[100]
    Catchps[i,j] <- ps$catch[100]
    fitps[i,j] <- ps$fit.adult[100]
    BSmarkps[i,j] <- ps$BS.mark
  }#End of  for (j in 1:100){
}#End of for(i in 1:100){


# Function to plot contours plots
plot.ppnRR.contours <- function(output.label, outputph, outputNAph, outputpm, 
                                outputNApm, outputps, outputNAps){
  #Contour plots of management levers
  par(mar=c(3,4,3,2), mfrow=c(2,2), oma=c(2,1,2,0.5))
  
  #Where are NAs?
  
  BS.mark.NAph <- outputNAph
  Data_DF <- data.frame (Ppn.RR = numeric(10000), Hatchery.size = numeric(10000), 
                         z = numeric(10000) )
  Data_DF$Ppn.RR <- rep(1:100*0.005, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(BS.mark.NAph)
  InterpListBSNAph <- interp(Data_DF$Ppn.RR, Data_DF$Hatchery.size, Data_DF$z, 
                             nx=200, ny=200)
  
  BS.mark.NApm <- outputNApm
  Data_DF <- data.frame (Ppn.RR = numeric(10000), Per.mark = numeric(10000), 
                         z = numeric(10000) )
  Data_DF$Ppn.RR <- rep(1:100*0.005, 100)
  Data_DF$Per.mark <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(BS.mark.NApm)
  InterpListBSNApm <- interp(Data_DF$Ppn.RR, Data_DF$Per.mark, Data_DF$z, 
                             nx=200, ny=200)
  
  BS.mark.NAps <- outputNAps
  Data_DF <- data.frame (Ppn.RR = numeric(10000), Sel = numeric(10000), 
                         z = numeric(10000) )
  Data_DF$Ppn.RR <- rep(1:100*0.005, 100)
  Data_DF$Sel <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(BS.mark.NAps)
  InterpListBSNAps <- interp(Data_DF$Ppn.RR, Data_DF$Sel, Data_DF$z, 
                             nx=200, ny=200)
  
  
  # For plot(a), Put data in vector form
  Data_DF <- data.frame (Ppn.RR=numeric(10000), Hatchery.size=numeric(10000), 
                         z=numeric(10000) )
  Data_DF$Ppn.RR <- rep(1:100*0.005, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputph)
  #write.csv(Data_DF, "ContourData_DF.csv")
  
  # Mask contours where target BS exceeds limit set by Ppn.RR
  InterpList <- interp (Data_DF$Ppn.RR, Data_DF$Hatchery.size, Data_DF$z, 
                        nx=200, ny=200)
  InterpList[[3]][which (InterpListBSNAph[[3]] > 0) ] <- NA
  
  # get x and y ranges for plotting
  xrange <- range(InterpList[[1]])
  yrange <- range(InterpList[[2]])
  zrange <- range(InterpList[[3]], na.rm=T)
  zlevels <- pretty(zrange,20)
  nlevels <- (length(zlevels))
  
  # make color palette
  cp <- colorRampPalette(brewer.pal(9,'Blues'))
  cp4 <- rev( cp( nlevels + 2)[-c(1, 2)])
  
  if (output.label=="pHOSeff") {
    cp4 <- (cp(nlevels + 2)[-c(1, 2)])
  } 
  cp_grey <- colorRampPalette(brewer.pal(9, 'Greys')) # Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)]) # Not implemented
  
  plot(NA, xlim=xrange, ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", 
       ylab="")
  .filled.contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
                   levels = pretty(zrange, nlevels) ,cp4) 
  if (output.label == "PNI" | output.label == "pNOB") {
    .filled.contour (x = InterpList[[1]], y = InterpList[[2]], 
                     z = InterpList[[3]], levels = c(0,0.5, 0.72,1), 
                     cp4[c(1, nlevels / 2, nlevels)])
  } 
  if (output.label=="pHOSeff") {
    .filled.contour(x=InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
                    levels = c(0, 0.28,0.5,1), cp4[ c( round( nlevels / 3), 
                                                       nlevels / 1.5, nlevels)])
  } 
  contour( x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
           labcex = 0.6, xlab = "", ylab = "", nlevels = 5, add = TRUE)
  if (output.label=="PNI"|output.label=="pNOB") {
    contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
             labcex = 0.6, xlab = "", ylab = "", level = c(0.72, 0.5), 
             add = TRUE, 
             lty = c("dashed", "dotted"), lwd = c(2, 2) ) }
  if (output.label == "pHOSeff") {
    contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
             labcex = 0.6, xlab = "", ylab = "", level = c(0.28, 0.5), 
             add = TRUE, 
             lty = c("dashed", "dotted"), lwd = c(2,2))
  }
  usr <- par("usr")
  axis(1); axis(2)
  int <- pretty(ph.BS, n=3)
  axis( side = 4, at = c(0, length (which (ph.BS < int[2]) ) * 0.005, 
                         length (which (ph.BS < int[3]) ) * 0.005, 
                         length (which (ph.BS < int[4]) ) * 0.005), 
        labels = NA, tcl = -0.3, cex.axis = 0.8)
  mtext (side = 4, text = int[1], line = 0.3, at = 0, cex = 0.8)
  mtext (side = 4, text = int[2], line = 0.3, 
         at = length (which (ph.BS < int[2]) ) * 0.005, cex = 0.8)
  mtext (side = 4, text = int[3], line = 0.3, 
         at = length (which (ph.BS < int[3]) ) * 0.005, cex = 0.8)
  mtext( side = 4, text = int[4], line = 0.3, 
         at = length (which (ph.BS < int[4]) ) * 0.005, cex = 0.8)
  mtext(side = 4, text = "Brood stock", line = 1, cex = 0.7)
  abline(h = usr[3], v = usr[1])
  abline(v = usr[2])
  mtext(side = 1, text = "Max escapement removal (ppn returns)", line = 2.5, 
        cex = 0.8)
  mtext(side = 2, text = "Hatchery size\n", line = 1.1, cex = 0.8)
  mtext(side = 3, text = "(100% marking,no selective removals)", line = 0.1, 
        cex = 0.6, at = 0.5, adj = 1)
  mtext(side = 3, text = "(a)", line = 0.1, cex = 0.8, at = 0, adj = 0)
  
  # For plot(b), Put data in vector form
  Data_DF <- data.frame (Ppn.RR=numeric(10000), Per.mark=numeric(10000), 
                         z=numeric(10000) )
  Data_DF$Ppn.RR <- rep(1:100*0.005, 100)
  Data_DF$Per.mark <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(outputpm)
  #write.csv(Data_DF, "ContourData_DF.csv")
  
  # Mask contours where target BS exceeds limit set by Ppn.RR
  InterpList <- interp (Data_DF$Ppn.RR, Data_DF$Per.mark, Data_DF$z, 
                        nx=200, ny=200)
  InterpList[[3]][which (InterpListBSNApm[[3]] > 0) ] <- NA
  
  # get x and y ranges for plotting
  xrange <- range(InterpList[[1]])
  yrange <- range(InterpList[[2]])
  zrange <- range(InterpList[[3]], na.rm=T)
  zlevels <- pretty(zrange,20)
  nlevels <- (length(zlevels))
  
  # make color palette
  cp <- colorRampPalette(brewer.pal(9,'Blues'))
  cp4 <- rev( cp( nlevels + 2)[-c(1, 2)])
  
  if (output.label=="pHOSeff") {
    cp4 <- (cp(nlevels + 2)[-c(1, 2)])
  } 
  cp_grey <- colorRampPalette(brewer.pal(9, 'Greys')) # Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)]) # Not implemented
  
  plot(NA, xlim=xrange, ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", 
       ylab="")
  .filled.contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
                   levels = pretty(zrange, nlevels) ,cp4) 
  if (output.label == "PNI" | output.label == "pNOB") {
    .filled.contour (x = InterpList[[1]], y = InterpList[[2]], 
                     z = InterpList[[3]], levels = c(0,0.5, 0.72,1), 
                     cp4[c(1, nlevels / 2, nlevels)])
  } 
  if (output.label=="pHOSeff") {
    .filled.contour(x=InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
                    levels = c(0, 0.28,0.5,1), cp4[ c( round( nlevels / 3), 
                                                       nlevels / 1.5, nlevels)])
  } 
  contour( x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
           labcex = 0.6, xlab = "", ylab = "", nlevels = 5, add = TRUE)
  if (output.label=="PNI"|output.label=="pNOB") {
    contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
             labcex = 0.6, xlab = "", ylab = "", level = c(0.72, 0.5), 
             add = TRUE, 
             lty = c("dashed", "dotted"), lwd = c(2, 2) ) }
  if (output.label == "pHOSeff") {
    contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
             labcex = 0.6, xlab = "", ylab = "", level = c(0.28, 0.5), 
             add = TRUE, 
             lty = c("dashed", "dotted"), lwd = c(2,2))
  }
  usr <- par("usr")
  axis(1); axis(2)
  int <- pretty(ph.BS, n=3)
  abline(h = usr[3], v = usr[1])
  abline(v = usr[2])
  mtext(side = 1, text = "Max escapement removal (ppn returns)", line = 2.5, 
        cex = 0.8)
  mtext(side = 2, text = "Percent marked\n", line = 1.1, cex = 0.8)
  mtext(side = 3, text = "(Hatchery size = 0.15\nno selective removals)", 
        line = 0.1, cex = 0.6, at = 0.5, adj = 1)
  mtext(side = 3, text = "(b)", line = 0.1, cex = 0.8, at = 0, adj = 0)
  if (output.label=="PNI")  {
    mtext( side=3, text=paste("PNI", title, sep=""), line = 1.5, outer = F, 
           cex = 1.3, at = -0.2)
  }
  
  # For plot(c), Put data in vector form
  Data_DF <- data.frame (Ppn.RR=numeric(10000), Sel=numeric(10000), 
                         z=numeric(10000) )
  Data_DF$Ppn.RR <- rep(1:100*0.005, 100)
  Data_DF$Sel <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(outputps)
  #write.csv(Data_DF, "ContourData_DF.csv")
  
  # Mask contours where target BS exceeds limit set by Ppn.RR
  InterpList <- interp (Data_DF$Ppn.RR, Data_DF$Sel, Data_DF$z, 
                        nx=200, ny=200)
  InterpList[[3]][which (InterpListBSNAps[[3]] > 0) ] <- NA
  
  # get x and y ranges for plotting
  xrange <- range(InterpList[[1]])
  yrange <- range(InterpList[[2]])
  zrange <- range(InterpList[[3]], na.rm=T)
  zlevels <- pretty(zrange,20)
  nlevels <- (length(zlevels))
  
  # make color palette
  cp <- colorRampPalette(brewer.pal(9,'Blues'))
  cp4 <- rev( cp( nlevels + 2)[-c(1, 2)])
  
  if (output.label=="pHOSeff") {
    cp4 <- (cp(nlevels + 2)[-c(1, 2)])
  } 
  cp_grey <- colorRampPalette(brewer.pal(9, 'Greys')) # Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)]) # Not implemented
  
  plot(NA, xlim=xrange, ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", 
       ylab="")
  .filled.contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
                   levels = pretty(zrange, nlevels) ,cp4) 
  if (output.label == "PNI" | output.label == "pNOB") {
    .filled.contour (x = InterpList[[1]], y = InterpList[[2]], 
                     z = InterpList[[3]], levels = c(0,0.5, 0.72,1), 
                     cp4[c(1, nlevels / 2, nlevels)])
  } 
  if (output.label=="pHOSeff") {
    .filled.contour(x=InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
                    levels = c(0, 0.28,0.5,1), cp4[ c( round( nlevels / 3), 
                                                       nlevels / 1.5, nlevels)])
  } 
  contour( x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
           labcex = 0.6, xlab = "", ylab = "", nlevels = 5, add = TRUE)
  if (output.label=="PNI"|output.label=="pNOB") {
    contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
             labcex = 0.6, xlab = "", ylab = "", level = c(0.72, 0.5), 
             add = TRUE, 
             lty = c("dashed", "dotted"), lwd = c(2, 2) ) }
  if (output.label == "pHOSeff") {
    contour (x = InterpList[[1]], y = InterpList[[2]], z = InterpList[[3]], 
             labcex = 0.6, xlab = "", ylab = "", level = c(0.28, 0.5), 
             add = TRUE, 
             lty = c("dashed", "dotted"), lwd = c(2,2))
  }
  usr <- par("usr")
  axis(1); axis(2)
  int <- pretty(ph.BS, n=3)
  abline(h = usr[3], v = usr[1])
  abline(v = usr[2])
  mtext(side = 1, text = "Max escapement removal (ppn returns)", line = 2.5, 
        cex = 0.8)
  mtext(side = 2, text = "Percent selectively havested\n", line = 1.1, 
        cex = 0.8)
  mtext(side = 3, text = "(Hatchery size = 0.15, 100% marking)", line = 0.1, 
        cex = 0.6, at = 0.5, adj = 1)
  mtext(side = 3, text = "(c)", line = 0.1, cex = 0.8, at = 0, adj = 0)
  
} #end of plot.ppnRR.contours

output.label<-"PNI"  #PNI, pHOS, pNOB, RperS, Ret.nat
outputph <- PNIph
outputNAph <- BSmarkph
outputpm <- PNIpm
outputNApm <- BSmarkpm
outputps <- PNIps
outputNAps <- BSmarkps
scenario <- ""
dir <- here::here("Results", "PpnReturnsRiver", "MainResults")

png(paste(dir, "/", output.label, scenario, ".png", sep=""), width=6, height=6, 
    units="in", res=1000)
plot.ppnRR.contours( output.label, outputph, outputNAph, outputpm, outputNApm, 
                     outputps, outputNAps)
dev.off()
