#Contour plots for sensitivity analyses on management levers
#Date last revised 5 Apr 2017


Run.contours<-function(c, percent.hatch, sex.ratio, ppn.RR, HR, h, w, RS, scenario=NA, title=NA, dir=NA){  

PNImh<-matrix(NA,100,100); pHOSmh<-matrix(NA,100,100); pNOBmh<-matrix(NA,100,100);  RperSmh<-matrix(NA,100,100); Ret.natmh<-matrix(NA,100,100); Ret.hatmh<-matrix(NA,100,100); Catchmh<-matrix(NA,100,100); fitmh<-matrix(NA,100,100); BSmarkmh<-matrix(NA,100,100);
PNIsh<-matrix(NA,100,100); pHOSsh<-matrix(NA,100,100); pNOBsh<-matrix(NA,100,100);  RperSsh<-matrix(NA,100,100); Ret.natsh<-matrix(NA,100,100); Ret.hatsh<-matrix(NA,100,100); Catchsh<-matrix(NA,100,100); fitsh<-matrix(NA,100,100); BSmarksh<-matrix(NA,100,100);
PNIms<-matrix(NA,100,100); pHOSms<-matrix(NA,100,100); pNOBms<-matrix(NA,100,100);  RperSms<-matrix(NA,100,100); Ret.natms<-matrix(NA,100,100); Ret.hatms<-matrix(NA,100,100); Catchms<-matrix(NA,100,100); fitms<-matrix(NA,100,100); BSmarkms<-matrix(NA,100,100);
#Repeat selective harvest vs hatchery size relationship for 100% marking case
PNIsh100<-matrix(NA,100,100); pHOSsh100<-matrix(NA,100,100); pNOBsh100<-matrix(NA,100,100);  RperSsh100<-matrix(NA,100,100); Ret.natsh100<-matrix(NA,100,100); Ret.hatsh100<-matrix(NA,100,100); Catchsh100<-matrix(NA,100,100); fitsh100<-matrix(NA,100,100); BSmarksh100<-matrix(NA,100,100);
mh.hatchery.size<-NA; mh.sel<-NA; mh.BS<-NA
sh.hatchery.size<-NA; sh.per.mark<-NA; sh.BS<-NA
ms.hatchery.size<-NA; ms.BS<-NA
sh100.hatchery.size<-NA; sh100.per.mark<-NA; sh100.BS<-NA

for(i in 1:100){
  print(i)
  for (j in 1:100){
  mh<-run.lever.model(per.mark=i*0.01, hatchery.size=j*0.005, sel=0, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=RS, mar.surv.hatch=0.0024, sex.ratio=sex.ratio, ppn.RR=ppn.RR)
  mh.sel<-mh$sel
  if (i==1){mh.hatchery.size[j]<-mh$hatchery.size; mh.BS[j]<-mh$BS[100]}
  PNImh[i,j]<-mh$PNI[100]
  pHOSmh[i,j]<-mh$pHOSeff[100]
  pNOBmh[i,j]<-mh$pNOB[100]
  RperSmh[i,j]<-mh$RperS[100]
  Ret.natmh[i,j]<-mh$ret.nat.preharvest[100]
  Ret.hatmh[i,j]<-mh$ret.hatch.preharvest[100]
  Catchmh[i,j]<-mh$catch[100]
  fitmh[i,j]<-mh$fit.adult[100]
  BSmarkmh[i,j]<-mh$BS.mark
  
  sh<-run.lever.model(per.mark=0.5, hatchery.size=j*0.005, sel=i*0.01, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=RS, mar.surv.hatch=0.0024, sex.ratio=sex.ratio, ppn.RR=ppn.RR)
  sh.per.mark<-sh$per.mark
  if (i==1){sh.hatchery.size[j]<-sh$hatchery.size;  sh.BS[j]<-sh$BS[100]}
  PNIsh[i,j]<-sh$PNI[100]
  pHOSsh[i,j]<-sh$pHOSeff[100]
  pNOBsh[i,j]<-sh$pNOB[100]
  RperSsh[i,j]<-sh$RperS[100]  
  Ret.natsh[i,j]<-sh$ret.nat.preharvest[100]  
  Ret.hatsh[i,j]<-sh$ret.hatch.preharvest[100]  
  Catchsh[i,j]<-sh$catch[100]  
  fitsh[i,j]<-sh$fit.adult[100]
  BSmarksh[i,j]<-sh$BS.mark
  
  ms<-run.lever.model(per.mark=i*0.01, hatchery.size=0.15, sel=j*0.01, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=RS, mar.surv.hatch=0.0024, sex.ratio=sex.ratio, ppn.RR=ppn.RR)
  ms.hatchery.size<-ms$hatchery.size
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
  
  #Repeat selective harvest vs hatchery size relationship for 100% marking case
  sh100<-run.lever.model(per.mark=1.0, hatchery.size=j*0.005, sel=i*0.01, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=RS, mar.surv.hatch=0.0024, sex.ratio=sex.ratio, ppn.RR=ppn.RR)
  sh100.per.mark<-sh100$per.mark
  if (i==1){sh100.hatchery.size[j]<-sh100$hatchery.size;  sh100.BS[j]<-sh100$BS[100]}
  PNIsh100[i,j]<-sh100$PNI[100]
  pHOSsh100[i,j]<-sh100$pHOSeff[100]
  pNOBsh100[i,j]<-sh100$pNOB[100]
  RperSsh100[i,j]<-sh100$RperS[100]  
  Ret.natsh100[i,j]<-sh100$ret.nat.preharvest[100]  
  Ret.hatsh100[i,j]<-sh100$ret.hatch.preharvest[100]  
  Catchsh100[i,j]<-sh100$catch[100]  
  fitsh100[i,j]<-sh100$fit.adult[100]
  BSmarksh100[i,j]<-sh100$BS.mark
  }#End of  for (j in 1:100){
}#End of for(i in 1:100){

#Contour plots of management levers


Do.contours<-function(){
par(mar=c(3,4,3,2), mfrow=c(2,2), oma=c(2,1,2,0.5))
#png("PNIonlyw4029Aug2017.png", width=4, height=4, units="in", res=1000)
#par(mar=c(3,4,3,2), mfrow=c(1,1), oma=c(2,1,2,0.5))

#Where are NAs in sh100?
  BS.mark.NAsh100<-outputNAsh100
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(BS.mark.NAsh100)
  InterpListBSNAsh100 <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)

  BS.mark.NAsh<-outputNAsh
  Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Sel <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(BS.mark.NAsh)
  InterpListBSNAsh <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  
  BS.mark.NAmh<-outputNAmh
  Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
  Data_DF$z <- as.vector(BS.mark.NAmh)
  InterpListBSNAmh <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
  
  BS.mark.NAms<-outputNAms
  Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
  Data_DF$Per.mark <- rep(1:100*0.01, 100)
  Data_DF$Sel <- rep(1:100*0.01, each = 100)
  Data_DF$z <- as.vector(BS.mark.NAms)
  InterpListBSNAms <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=200, ny=200)
  
  # Put data in vector form
Data_DF <- data.frame(Per.mark=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
Data_DF$Per.mark <- rep(1:100*0.01, 100)
Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
Data_DF$z <- as.vector(outputmh)
#write.csv(Data_DF, "ContourData_DF.csv")

InterpList <- interp(Data_DF$Per.mark, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
InterpList[[3]][which(InterpListBSNAmh[[3]]>0)]<-NA

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
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.8,1) ,cp4[c(1,nlevels/2, nlevels)])} 
if(output.label=="pHOSeff")#"pHOSeff")
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.19,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=5, add=TRUE)
if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.8, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.19, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
usr<-par("usr")
axis(1); axis(2)
int<-pretty(mh.BS,n=3)
axis(side=4, at=c(0, length(which(mh.BS<int[2]))*0.005, length(which(mh.BS<int[3]))*0.005, length(which(mh.BS<int[4]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.8)
mtext(side=4, text=int[1], line=0.3, at=0, cex=0.8)
mtext(side=4, text=int[2], line=0.3, at=length(which(mh.BS<int[2]))*0.005, cex=0.8)
mtext(side=4, text=int[3], line=0.3, at=length(which(mh.BS<int[3]))*0.005, cex=0.8)
mtext(side=4, text=int[4], line=0.3, at=length(which(mh.BS<int[4]))*0.005, cex=0.8)
mtext(side=4, text="Brood stock", line=1, cex=0.7)
abline(h=usr[3], v=usr[1])
abline(v=usr[2])
mtext(side=1, text="Proportion marked", line=2.5, cex=0.8)
mtext(side=2, text="Hatchery size\n", line=1.1, cex=0.8)
mtext(side=3, text="(No selective removals of marked fish)", line=0.1, cex=0.6, at=1, adj=1)
mtext(side=3, text="(a)", line=0.1, cex=0.8, at=0, adj=0)
#mtext(3, text=paste(output.label, sep=""), line=2, outer=T, cex=1.5)



# Variables (management levers) to assess: Sel vs Hatchery size
# Put data in vector form
Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
Data_DF$Sel <- rep(1:100*0.01, 100)
Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
Data_DF$z <- as.vector(outputsh)
#write.csv(Data_DF, "ContourData_DF.csv")

InterpList <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
InterpList[[3]][which(InterpListBSNAsh[[3]]>0)]<-NA

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
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.8,1) ,cp4[c(1,nlevels/2, nlevels)])} 
if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.19,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=5, add=TRUE)
if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.8, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.19, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
usr<-par("usr")
axis(1); axis(2)
int<-pretty(mh.BS,n=3)
axis(side=4, at=c(0, length(which(sh.BS<int[2]))*0.005, length(which(sh.BS<int[3]))*0.005, length(which(sh.BS<int[4]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.8)
mtext(side=4, text=int[1], line=0.3, at=0, cex=0.8)
mtext(side=4, text=int[2], line=0.3, at=length(which(mh.BS<int[2]))*0.005, cex=0.8)
mtext(side=4, text=int[3], line=0.3, at=length(which(mh.BS<int[3]))*0.005, cex=0.8)
mtext(side=4, text=int[4], line=0.3, at=length(which(mh.BS<int[4]))*0.005, cex=0.8)
mtext(side=4, text="Brood stock", line=1, cex=0.7)
abline(h=usr[3], v=usr[1])
abline(v=usr[2])
mtext(side=1, text="Proportion of marked fish\nselectively removed", line=2.8, cex=0.8)
mtext(side=2, text="Hatchery size\n", line=1.1, cex=0.8)
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

# Variables (management levers) to assess: %Marking vs selectivity
# Put data in vector form
Data_DF <- data.frame(Per.mark=numeric(10000), Sel=numeric(10000), z=numeric(10000) )
Data_DF$Per.mark <- rep(1:100*0.01, 100)
Data_DF$Sel <- rep(1:100*0.01, each=100)
Data_DF$z <- as.vector(outputms)

# InterpList <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=15, ny=15)
InterpList <- interp(Data_DF$Per.mark, Data_DF$Sel, Data_DF$z, nx=200, ny=200)
InterpList[[3]][which(InterpListBSNAms[[3]]>0)]<-NA

# If all NAs, then skip (Hatchery size of 15% too big), where 200 x 200 is 
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
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.8,1), cp4[c(1,nlevels/2, nlevels)])} 
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.19,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
  contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=6, add=TRUE)
  if(output.label=="PNI"|output.label=="pNOB")
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.8, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
  {contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.19, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
  usr<-par("usr")
  axis(1); axis(2)
  abline(h=usr[3], v=usr[1])
  mtext(side=1, text="Proportion marked", line=2.5, cex=0.8)
  mtext(side=2, text="Proportion marked fish\nselectively removed", line=2, cex=0.8)
  mtext(side=3, text="(Hatchery size = 0.15)", line=0.1, cex=0.6, at=1, adj=1)
  mtext(side=3, text="(c)", line=0.1, cex=0.8, at=0, adj=0)
  #mtext(side=3, text=paste(output.label, sep=""), line=2, outer=F, cex=1.5, at=0.5)
  if(output.label=="Recruits from hatchery production"){text(x=0.1, y=0.95, labels=round(Ret.hatms[1],0), cex=0.8)}
  
  }  

# Variables (management levers) to assess: Sel vs Hatchery size for 100% marking scenario
# Put data in vector form
Data_DF <- data.frame(Sel=numeric(10000), Hatchery.size=numeric(10000), z=numeric(10000) )
Data_DF$Sel <- rep(1:100*0.01, 100)
Data_DF$Hatchery.size <- rep(1:100*0.005, each=100)
Data_DF$z <- as.vector(outputsh100)
#write.csv(Data_DF, "ContourData_DF.csv")

InterpList <- interp(Data_DF$Sel, Data_DF$Hatchery.size, Data_DF$z, nx=200, ny=200)
InterpList[[3]][which(InterpListBSNAsh100[[3]]>0)]<-NA


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
{.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0,0.5, 0.8,1) ,cp4[c(1,nlevels/2, nlevels)])} 
if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
{.filled.contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], levels = c(0, 0.19,0.5,1) ,cp4[c(round(nlevels/3), nlevels/1.5, nlevels)])} 
contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", nlevels=5, add=TRUE)
if(output.label=="PNI"|output.label=="pNOB")
{contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.8, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
if(output.label=="pHOSeff")#(output.label=="pHOSeff") 
{contour(x=InterpList[[1]],y=InterpList[[2]], z=InterpList[[3]], labcex=0.6, xlab="", ylab="", level=c(0.19, 0.5), add=TRUE, lty=c("dashed", "dotted"), lwd=c(2,2))}
usr<-par("usr")
axis(1); axis(2)
int<-pretty(sh100.BS,n=3)
axis(side=4, at=c(0, length(which(sh100.BS<int[2]))*0.005, length(which(sh100.BS<int[3]))*0.005, length(which(sh100.BS<int[4]))*0.005), labels=NA, tcl=-0.3, cex.axis=0.8)
mtext(side=4, text=int[1], line=0.3, at=0, cex=0.8)
mtext(side=4, text=int[2], line=0.3, at=length(which(sh100.BS<int[2]))*0.005, cex=0.8)
mtext(side=4, text=int[3], line=0.3, at=length(which(sh100.BS<int[3]))*0.005, cex=0.8)
mtext(side=4, text=int[4], line=0.3, at=length(which(sh100.BS<int[4]))*0.005, cex=0.8)
mtext(side=4, text="Brood stock", line=1, cex=0.7)
abline(h=usr[3], v=usr[1])
abline(v=usr[2])
mtext(side=1, text="Proportion marked fish\nselectively removed", line=2.8, cex=0.8)
mtext(side=2, text="Hatchery size\n", line=1.1, cex=0.8)
mtext(side=3, text="(100% marking)", line=0.1, cex=0.6, at=1, adj=1)
mtext(side=3, text="(d)", line=0.1, cex=0.8, at=0, adj=0)
if(output.label=="pNOB"){if(pNOBsh100[1,100]==1){text(x=0.1, y=0.48, labels=c("1.0"), cex=0.7)}}#round(pNOBsh100[1],1)

}#End of Do.contours()


#Run Contour Plots
#png("PNIbasecase6April2017.png", width=6, height=6, units="in", res=1000)
#png("PNIhigheritability6April2017.png", width=6, height=6, units="in", res=1000)
#png("PNIweakselection6April2017.png", width=6, height=6, units="in", res=1000)
# png("PNIrs0.228Nov2017.png", width=6, height=6, units="in", res=1000)
png(paste(dir, "PNI", scenario, ".png", sep=""), width=6, height=6, units="in", res=1000)
outputmh<-PNImh; outputsh<-PNIsh; outputms<-PNIms; outputsh100<-PNIsh100; outputNAsh100<-BSmarksh100; outputNAsh<-BSmarksh;  outputNAmh<-BSmarkmh;  outputNAms<-BSmarkms; output.label<-"PNI"  #PNI, pHOS, pNOB, RperS, Ret.nat
Do.contours()
dev.off()
# png("test.png"); plot(x=1:10, y=1:10); dev.off()

#png("pHOSbasecase6April2017.png", width=6, height=6, units="in", res=1000)
#png("pHOShigheritability6April2017.png", width=6, height=6, units="in", res=1000)
#png("pHOSweakselection6April2017.png", width=6, height=6, units="in", res=1000)
# png("pHOSrs0.228Nov2017.png", width=6, height=6, units="in", res=1000)
png(paste(dir, "pHOS", scenario, ".png", sep=""), width=6, height=6, units="in", res=1000)
outputmh<-pHOSmh; outputsh<-pHOSsh; outputms<-pHOSms; outputsh100<-pHOSsh100; outputNAsh100<-BSmarksh100; outputNAsh<-BSmarksh; outputNAmh<-BSmarkmh; outputNAms<-BSmarkms; output.label<-"pHOSeff"  #pHOSeff #PNI, pHOS, pNOB, RperS, Ret.nat
Do.contours()
dev.off()

#png("pNOBbasecase6April2017.png", width=6, height=6, units="in", res=1000)
#png("pNOBhigheritability6April2017.png", width=6, height=6, units="in", res=1000)
#png("pNOBweakselection6April2017.png", width=6, height=6, units="in", res=1000)
# png("pNOBrs0.228Nov2017.png", width=6, height=6, units="in", res=1000)
png(paste(dir, "pNOB", scenario, ".png", sep=""), width=6, height=6, units="in", res=1000)
outputmh<-pNOBmh; outputsh<-pNOBsh; outputms<-pNOBms; outputsh100<-pNOBsh100; outputNAsh100<-BSmarksh100; outputNAsh<-BSmarksh; outputNAmh<-BSmarkmh; outputNAms<-BSmarkms; output.label<-"pNOB"  #PNI, pHOS, pNOB, RperS, Ret.nat
Do.contours()
dev.off()

#png("Retnatbasecase6April2017.png", width=6, height=6, units="in", res=1000)
#png("Retnathigheritability6April2017.png", width=6, height=6, units="in", res=1000)
#png("Retnatweakselection6April2017.png", width=6, height=6, units="in", res=1000)
# png("Retnatrs0.228Nov2017.png", width=6, height=6, units="in", res=200)
png(paste(dir, "Retnat", scenario, ".png", sep=""), width=6, height=6, units="in", res=200)
outputmh<-Ret.natmh; outputsh<-Ret.natsh; outputms<-Ret.natms; outputsh100<-Ret.natsh100; outputNAsh100<-BSmarksh100; outputNAsh<-BSmarksh; outputNAmh<-BSmarkmh; outputNAms<-BSmarkms; output.label<-"Recruits from river spawners (HOS+NOS)"  #PNI, pHOS, pNOB, RperS, Ret.nat
Do.contours()
dev.off()

#png("Rethatbasecase6April2017.png", width=6, height=6, units="in", res=1000)
#png("Rethathigheritability6April2017.png", width=6, height=6, units="in", res=1000)
#png("Rethatweakselection6April2017.png", width=6, height=6, units="in", res=1000)
# png("Rethatrs0.228Nov2017.png", width=6, height=6, units="in", res=1000)
png(paste(dir, "Rethat", scenario, ".png", sep=""), width=6, height=6, units="in", res=1000)
outputmh<-Ret.hatmh; outputsh<-Ret.hatsh; outputms<-Ret.hatms; outputsh100<-Ret.hatsh100; outputNAsh100<-BSmarksh100; outputNAsh<-BSmarksh; outputNAmh<-BSmarkmh; output.label<-"Recruits from hatchery production"  #PNI, pHOS, pNOB, RperS, Ret.nat
Do.contours()
dev.off()

#png("Catchbasecase6April2017.png", width=6, height=6, units="in", res=1000)
#png("Catchhigheritability6April2017.png", width=6, height=6, units="in", res=1000)
#png("Catchweakselection6April2017.png", width=6, height=6, units="in", res=1000)
# png("Catchrs0.228Nov2017.png", width=6, height=6, units="in", res=1000)
png(paste(dir, "Catch", scenario, ".png", sep=""), width=6, height=6, units="in", res=1000)
outputmh<-Catchmh; outputsh<-Catchsh; outputms<-Catchms; outputsh100<-Catchsh100; outputNAsh100<-BSmarksh100; outputNAsh<-BSmarksh; outputNAmh<-BSmarkmh; output.label<-"Catch"  #PNI, pHOS, pNOB, RperS, Ret.nat
Do.contours()
dev.off()

#outputmh<-RperSmh; outputsh<-RperSsh; outputms<-RperSms; outputsh100<-RperSsh100; output.label<-"Recruits/spawner"  #PNI, pHOS, pNOB, RperS, Ret.nat
#Do.contours()



}#End of Run.contours<-function(c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w){
