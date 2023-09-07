#***************************************************************************
#Contour plots for sensitivity analyses on management levers
#Date last revised  6 July 2023
#Uses  data.frames of results for 3 combinations of management levers (3 rows), and 4 combinations of marine survival (high and low for natural/hatchery pop)
# To create for 3x4 plot
# Requires output from ContourBSPpnRR_ms.r
#***************************************************************************

#Source("ContourBSPpnRRSA_ms.r")


# Need to install older versin of akima for contour plots
# detach("package:akima", unload=TRUE)
# install.packages("C:/Users/HoltC/AppData/Local/Programs/R/R-4.2.1/library/akima_0.6-2.tar.gz", repos = NULL, type="source")


plot.3x4<-function(metric){#metric=="PNI", "Ret.nat", "Ret.hat"
  
par(mfcol=c(3,4), oma=c(2,1,2,0.5))
#pdf("test.pdf"); par(mfcol=c(3,4), oma=c(2,1,2,0.5))
for (ms.scenario in 1:4){
  if(ms.scenario==1){par(mar=c(1.5,5,2.5,1))}
  if(ms.scenario==1|ms.scenario==2|ms.scenario==3){par(mar=c(1.5,3,2.5,1))}
  if(ms.scenario==4){par(mar=c(1.5,3,2.5,2))}
  
  #Contour plots of management levers
  if(ms.scenario==1){input<-lowmshatch}#lowms}#lowms.lowmshatch}
  if(ms.scenario==2){input<-highmshatch}#highms}#lowms.highmshatch}
  if(ms.scenario==3){input<-lowms}#lowmshatch}#highms.lowmshatch}
  if(ms.scenario==4){input<-highms}#highmshatch}#highms.highmshatch}
  
    
  if(metric=="PNI"){outputmh<-input$PNImh; outputsh<-input$PNIsh; outputms<-input$PNIms; output.label<-"PNI"}  #PNI, pHOS, pNOB, RperS, Ret.nat
  if(metric=="Ret.nat"){outputmh<-input$Ret.natmh; outputsh<-input$Ret.natsh; outputms<-input$Ret.natms; output.label<-"Recruits from river spawners (HOS+NOS)"}  #PNI, pHOS, pNOB, RperS, Ret.nat
  if(metric=="Ret.hat"){outputmh<-input$Ret.hatmh; outputsh<-input$Ret.hatsh; outputms<-input$Ret.hatms; output.label<-"Recruits from hatchery production"}  #PNI, pHOS, pNOB, RperS, Ret.nat
  if(metric=="pHOSeff"){outputmh<-input$pHOSmh; outputsh<-input$pHOSsh; outputms<-input$pHOSms; output.label<-"pHOSeff"}  #PNI, pHOS, pNOB, RperS, Ret.nat
  if(metric=="pNOB"){outputmh<-input$pNOBmh; outputsh<-input$pNOBsh; outputms<-input$pNOBms; output.label<-"pNOB"}  #PNI, pHOS, pNOB, RperS, Ret.nat
  
#Keep Ret.nat for all management levers for plotting NA region (where population is extinct) on plots of other metrics
  Ret.nat.NAmh<-input$Ret.natmh
  Ret.nat.NAsh<-input$Ret.natsh
  Ret.nat.NAms<-input$Ret.natms
  BS.mark.NAsh<-input$BSmarksh
  BS.mark.NAmh<-input$BSmarkmh
  BS.mark.NAms<-input$BSmarkms
  HL.mark.NAsh<-input$HLmarksh
  HL.mark.NAmh<-input$HLmarkmh
  HL.mark.NAms<-input$HLmarkms
  
  indmh<-NA;indsh<-NA; indms<-NA
# Variables (management levers) to assess: %Marked vs Hatchery size
  Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(Ret.nat.NAmh)
  InterpListNAmh <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=40, ny=40)
# Variables (management levers) to assess: Sel vs Hatchery size
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(Ret.nat.NAsh)
  InterpListNAsh <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=40, ny=40)
# Variables (management levers) to assess: %Marking vs selectivity
  Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Sel <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(Ret.nat.NAms)
  InterpListNAms <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=40, ny=40)

 # NAs where brook take exceeds 33% of returns to river. Now masking these values (Feb 2023)
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(HL.mark.NAsh)
  InterpListHLNAsh <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=40, ny=40)

  Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(HL.mark.NAmh)
  InterpListHLNAmh <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=40, ny=40)
  
  Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Sel <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(HL.mark.NAms)
  InterpListHLNAms <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=40, ny=40)
  
  # Variables (management levers) to assess: %Marked vs Hatchery size
  # Put data in vector form
  Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputmh)
  
  InterpList <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=40, ny=40)
  # Mask where extinct and  where bs exceeds 33% of returns to river
  # InterpList[[3]][which(InterpListNAmh[[3]]<2|InterpListBSNAmh[[3]]>0)]<-NA  

  # get x and y ranges for plotting
  xrange<-range(InterpList[[1]])
  yrange<-range(InterpList[[2]])
  zrange<-range(InterpList[[3]], na.rm=T)
  zlevels<-pretty(zrange,20)
  nlevels<-(length(zlevels))
  
  # make color pallette
  cp<-colorRampPalette(brewer.pal(9,'Blues'))
  cp4<-rev(cp(nlevels+2)[-c(1,2)])
  if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}
  cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
  
  plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
  .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange, nlevels) ,cp4) 
  if(metric=="Ret.nat")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,2,pretty(zrange,nlevels)[2:nlevels]), c("black", cp4))}
  if(output.label=="PNI"|output.label=="pNOB")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1) ,cp4[c(1,nlevels/3, nlevels/1.5)])}#cp4[c(1,nlevels/2, nlevels)])} 
  if(output.label=="pHOSeff")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", nlevels=5, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="pHOSeff")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  usr<-par("usr")
  axis(1, cex.axis=0.9); axis(2, cex.axis=0.9)
  
  # Make vector of grids for hatching/stippling where NAs are
  incl <- which(InterpListHLNAmh[[3]]>0)
  
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
  

  # int<-pretty(input$mh.BS,n=3)
  # if(int[length(int)]>max(input$mh.BS)){int<-int[1:(length(int)-1)]}
  # if(length(int)==4){axis(side=4, at=c(0, length(which(input$mh.BS<int[2]))*0.005, length(which(input$mh.BS<int[3]))*0.005, length(which(input$mh.BS<int[4]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.6)}
  # if(length(int)==3){axis(side=4, at=c(0, length(which(input$mh.BS<int[2]))*0.005, length(which(input$mh.BS<int[3]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.6)}
  # mtext(side=4, text=int[1], line=0.3, at=0, cex=0.6)
  # mtext(side=4, text=int[2], line=0.3, at=length(which(input$mh.BS<int[2]))*0.005, cex=0.6)
  # mtext(side=4, text=int[3], line=0.3, at=length(which(input$mh.BS<int[3]))*0.005, cex=0.6)
  # if(length(int)==4)mtext(side=4, text=int[4], line=0.3, at=length(which(input$mh.BS<int[4]))*0.005, cex=0.6)
  scaledBS <- input$mh.BS/input$mh.Seq
  int <- signif(quantile(scaledBS, c(0.05,0.5,0.95)), 1)
  axis(side=4, at=c(length(which(scaledBS<int[1]))*0.005,
                    length(which(scaledBS<int[2]))*0.005, 
                    length(which(scaledBS<int[3]))*0.005),
       labels=NA, tcl=-0.3, cex.axis=0.8)
  mtext(side=4, text=int[1], line=0.3, at=length(which(scaledBS<int[1]))*0.005, cex=0.5)
  mtext(side=4, text=int[2], line=0.3, at=length(which(scaledBS<int[2]))*0.005, cex=0.5)
  mtext(side=4, text=int[3], line=0.3, at=length(which(scaledBS<int[3]))*0.005, cex=0.5)
  axis(side=4, at=(length(which(scaledBS<0.3))*0.005), 
       labels=NA, tcl=-0.3, cex.axis=0.8)
  mtext(side=4, text=0.3, line=0.3, at=length(which(scaledBS<0.3))*0.005, 
        cex=0.5)
  if(int[2]>0.3){ #Add tick mark at 1, if there is space
    axis(side=4, at=(length(which(scaledBS<1))*0.005), 
         labels=NA, tcl=-0.3, cex.axis=0.8)
    mtext(side=4, text=1, line=0.3, at=length(which(scaledBS<1))*0.005, cex=0.5)
    }
  if(ms.scenario==4){mtext(side=4, text="Hatchery size (ppn of Seq)", line=1, cex=0.6)}
  abline(h=usr[3], v=usr[1])
  abline(v=usr[2], h=usr[4])
  mtext(side=1, text="Proportion marked", line=2, cex=0.6)
  if(ms.scenario==1){mtext(side=2, text="Brood stock(ppn of returns)\n", line=2, cex=0.6)}
  if(ms.scenario==1)mtext(side=3,line=0.3, cex=0.5, at=0.6,  text="nat. surv.=0.02\nhatch. surv.=0.001")# text="(No selective removals of marked fish)", )
  if(ms.scenario==2)mtext(side=3,line=0.3, cex=0.5, at=0.6,  text="nat. surv.=0.02\nhatch. surv.=0.005")# text="(No selective removals of marked fish)", )
  if(ms.scenario==3)mtext(side=3,line=0.3, cex=0.5, at=0.6,  text="nat. surv.=0.01\nhatch. surv.=0.0024")# text="(No selective removals of marked fish)", )
  if(ms.scenario==4)mtext(side=3,line=0.3, cex=0.5, at=0.6,  text="nat. surv.=0.05\nhatch. surv.=0.0024")# text="(No selective removals of marked fish)", )
  if(ms.scenario==1){mtext(side=3, text="(a)", line=0.1, cex=0.6, at=0, adj=0)}
  if(ms.scenario==2){mtext(side=3, text="(d)", line=0.1, cex=0.6, at=0, adj=0)}
  if(ms.scenario==3){mtext(side=3, text="(g)", line=0.1, cex=0.6, at=0, adj=0)}
  if(ms.scenario==4){mtext(side=3, text="(j)", line=0.1, cex=0.6, at=0, adj=0)}
  #mtext(3, text=paste(output.label, sep=""), line=2, outer=T, cex=1.5)
  if(ms.scenario==2){
    if(output.label!="pHOSeff"){mtext(side=3, text=paste(output.label, sep=""), line=2.5, outer=F, cex=1, at=1.3)}
    if(output.label=="pHOSeff"){mtext(side=3, text=expression('pHOS'[eff]), line=2.5, outer=F, cex=1, at=1.3)}
  }
 
  
  # Variables (management levers) to assess: Sel vs Hatchery size
  # Put data in vector form
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(outputsh)
  
  InterpList <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=40, ny=40)
  # Mask where extinct and  where bs exceeds 33% of returns to river
  # InterpList[[3]][which(InterpListNAsh[[3]]<2|InterpListBSNAsh[[3]]>0)]<-NA

  # get x and y ranges for plotting
  xrange<-range(InterpList[[1]])
  yrange<-range(InterpList[[2]])
  zrange<-range(InterpList[[3]], na.rm=T)
  zlevels<-pretty(zrange,20)
  nlevels<-(length(zlevels))
  
  # make color pallette
  cp<-colorRampPalette(brewer.pal(9,'Blues'))
  cp4<-rev(cp(nlevels+2)[-c(1,2)])
  if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}
  cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
  
  plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
  .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange,nlevels) ,cp4) 
  if(metric=="Ret.nat")
    {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,2,pretty(zrange,nlevels)[2:nlevels]), c("black", cp4))}
  if(output.label=="PNI"|output.label=="pNOB")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1) ,cp4[c(1,nlevels/3, nlevels/1.5)])} 
  if(output.label=="pHOSeff")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", nlevels=5, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="pHOSeff")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))} 
  usr<-par("usr")
  axis(1, cex.axis=0.9); axis(2, cex.axis=0.9)
  
  
  # Make vector of grids for hatching/stippling where NAs are
  incl <- which(InterpListHLNAsh[[3]]>0)

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
  
  scaledBS <- input$sh.BS/input$sh.Seq
  int <- signif(quantile(scaledBS, c(0.05,0.5,0.95)), 1)
  axis(side=4, at=c(length(which(scaledBS<int[1]))*0.005,
                    length(which(scaledBS<int[2]))*0.005, 
                    length(which(scaledBS<int[3]))*0.005),
       labels=NA, tcl=-0.3, cex.axis=0.8)
  mtext(side=4, text=int[1], line=0.3, at=length(which(scaledBS<int[1]))*0.005, cex=0.5)
  mtext(side=4, text=int[2], line=0.3, at=length(which(scaledBS<int[2]))*0.005, cex=0.5)
  mtext(side=4, text=int[3], line=0.3, at=length(which(scaledBS<int[3]))*0.005, cex=0.5)
  axis(side=4, at=(length(which(scaledBS<0.3))*0.005), 
       labels=NA, tcl=-0.3, cex.axis=0.8)
  mtext(side=4, text=0.3, line=0.3, at=length(which(scaledBS<0.3))*0.005, 
        cex=0.5)
  if(int[2]>0.3){ #Add tick mark at 1, if there is space
    axis(side=4, at=(length(which(scaledBS<1))*0.005), 
         labels=NA, tcl=-0.3, cex.axis=0.8)
    mtext(side=4, text=1, line=0.3, at=length(which(scaledBS<1))*0.005, cex=0.5)
  }
  if(ms.scenario==4){mtext(side=4, text="Hatchery size (ppn of Seq)", line=1, cex=0.6)}
  
  # 
  # int<-pretty(input$sh.BS,n=3)
  # if(int[length(int)]>max(input$sh.BS)){int<-int[1:(length(int)-1)]}
  # if(length(int)==4){axis(side=4, at=c(0, length(which(input$sh.BS<int[2]))*0.005, length(which(input$sh.BS<int[3]))*0.005, length(which(input$sh.BS<int[4]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.6)}
  # if(length(int)==3){axis(side=4, at=c(0, length(which(input$sh.BS<int[2]))*0.005, length(which(input$sh.BS<int[3]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.6)}
  # #axis(side=4, at=c(0, length(which(input$sh.BS<int[2]))*0.005, length(which(input$sh.BS<int[3]))*0.005, length(which(input$sh.BS<int[4]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.8)
  # mtext(side=4, text=int[1], line=0.3, at=0, cex=0.6)
  # mtext(side=4, text=int[2], line=0.3, at=length(which(input$sh.BS<int[2]))*0.005, cex=0.6)
  # mtext(side=4, text=int[3], line=0.3, at=length(which(input$sh.BS<int[3]))*0.005, cex=0.6)
  # if(length(int)==4)mtext(side=4, text=int[4], line=0.3, at=length(which(input$sh.BS<int[4]))*0.005, cex=0.6)
  # if(ms.scenario==4){mtext(side=4, text="Brood stock(ppn of returns)", line=1, cex=0.6)}
  abline(h=usr[3], v=usr[1])
  abline(v=usr[2], h=usr[4])
  mtext(side=1, text="Proportion marked fish\nselectively harvested", line=2.7, cex=0.6)
  if(ms.scenario==1){mtext(side=2, text="Brood stock(ppn of returns)", line=2, cex=0.6)}
  #mtext(side=3, text="(50% marking)", line=0.1, cex=0.6, at=1, adj=1)
  if(ms.scenario==1)mtext(side=3, text="(b)", line=0.1, cex=0.6, at=0, adj=0)
  if(ms.scenario==2)mtext(side=3, text="(e)", line=0.1, cex=0.6, at=0, adj=0)
  if(ms.scenario==3)mtext(side=3, text="(h)", line=0.1, cex=0.6, at=0, adj=0)
  if(ms.scenario==4)mtext(side=3, text="(k)", line=0.1, cex=0.6, at=0, adj=0)
  #mtext(side=3, text=paste(output.label, sep=""), line=2, outer=F, cex=1.5, at=-0.2)
  
  
  
  
  # Variables (management levers) to assess: %Marking vs selectivity
  # Put data in vector form
  Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Sel <- rep(1:100*0.01, each=100)
  Data_DF$z <- as.vector(outputms)
  
  InterpList <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=40, ny=40)
  # Mask where extinct and  where bs exceeds 33% of returns to river
  # InterpList[[3]][which(InterpListNAms[[3]]<2|InterpListBSNAms[[3]]>0)]<-NA

  # get x and y ranges for plotting
  xrange<-range(InterpList[[1]])
  yrange<-range(InterpList[[2]])
  zrange<-range(InterpList[[3]], na.rm=T)
  zlevels<-pretty(zrange,20)
  nlevels<-(length(zlevels))
  
  # make color pallette
  cp<-colorRampPalette(brewer.pal(9,'Blues'))
  cp4<-rev(cp(nlevels+2)[-c(1,2)])
  if (output.label=="pHOSeff") {cp4<-(cp(nlevels+2)[-c(1,2)])}
  cp_grey <- colorRampPalette(brewer.pal(9,'Greys'))#Not implemented
  cp4_grey<-rev(cp_grey(nlevels+2)[-c(1,2)])#Not implemented
  
  plot(NA, xlim=xrange,ylim=yrange, frame=F, axes=F,xaxs="i", yaxs="i", xlab="", ylab="")
  .filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = pretty(zrange, nlevels) ,cp4) 
  if(metric=="Ret.nat")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,2,pretty(zrange,nlevels)[2:nlevels]), c("black", cp4))}
  if(output.label=="PNI"|output.label=="pNOB")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.72,1), cp4[c(1,nlevels/3, nlevels/1.5)])} 
  if(output.label=="pHOSeff")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.28,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", nlevels=6, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", level=c(0.72, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="pHOSeff")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.5, xlab="", ylab="", level=c(0.28, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  usr<-par("usr")
  axis(1, cex.axis=0.9); axis(2, cex.axis=0.9); 
  
  
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
    
  
  abline(h=usr[3], v=usr[1]); abline(h=usr[4], v=usr[2])
  mtext(side=1, text="Proportion marked", line=2, cex=0.6)
  if(ms.scenario==1){mtext(side=2, text="Proportion marked fish\nselectively harvested", line=2.0, cex=0.7)}
  #mtext(side=3, text="(Hatchery size = 0.1)", line=0.1, cex=0.6, at=1, adj=1)
  if(ms.scenario==1)mtext(side=3, text="(c)", line=0.1, cex=0.6, at=-0.15, adj=0)
  if(ms.scenario==2)mtext(side=3, text="(f)", line=0.1, cex=0.6, at=-0.15, adj=0)
  if(ms.scenario==3)mtext(side=3, text="(i)", line=0.1, cex=0.6, at=-0.15, adj=0)
  if(ms.scenario==4)mtext(side=3, text="(l)", line=0.1, cex=0.6, at=-0.2, adj=0)
  if(output.label=="Recruits from hatchery production"){text(x=0.1, y=0.95, labels=round(input$Ret.hatms[1],0), cex=0.8)}
  
  }#End of 4 ms.scenarios
#dev.off()#test pdf

}#End of plot.3x4()

#Run Contour Plots
png(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "PNI.png"), width=6, height=6, units="in", res=1000)
plot.3x4(metric="PNI")
dev.off()

#Repeat for "pHOSeff", "pNOB", "Ret.nat" "Ret.hat"
