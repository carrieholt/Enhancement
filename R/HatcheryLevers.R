#Code to estimate effect of 3 hatchery levers: harchery size, marking, and selective harvest
#Last updated 7 Dec. 2016
rm(list = ls())

source("R/FuncDefs.r")
source("R/SeqFormulationBH.r")
source("R/run.lever.model.r")

panel.plots.lever<-function(res=res, res.nogenetics=res.nogenetics){

par(mfrow=c(3,2), mar=c(4,4,2,1))
plot(1:100,res$fit.smolt^2, type="l", ylim=c(0,1), ylab="Fitness", xlab="Generation")
#points(x=1:100,y=res$P.nat/100, col="red",type="l")
#points(x=1:100,y=res$P.hatch/100, col="blue", type="l")
#abline(h=res$Theta.nat/100, col="red", lty="dotted")
#abline(h=res$Theta.hatch/100, col="blue", lty="dotted")
#legend(x=5,y=0.40, legend=c("Population fitness", "Phenotype: natural origin fish", "Phenotype: hatchery origin fish"), col=c("black", "red", "blue"), cex=0.7, lty=c("solid", "solid"), bty="n")
#text(x=2, y=0.93, labels="(a)")
title<-paste("Mark rate=",res$per.mark*100,"%;  Brood stock=",round(res$BS[1],0), " (",res$hatchery.size*100, "% ave natural recruitment);  Removal of marked fish=",res$sel*100,"%", sep="")
mtext(title, side=3, line=0.5, at=110, cex=0.75)

y.upper.lim<-max(res$Sp.nat, na.rm=T)#res$Seq*res$mar.surv*40
if(max(res$Sp.nat)>res$Seq*1.2){y.upper.lim<-max(res$Sp.nat)}

plot(1:100, res$Sp.nat, type="l", col="red", ylim=c(0,y.upper.lim), ylab="Spawners (HOS+NOS)", xlab="Generation")
points(x=1:100, y=res$NOS, col="red", type="l", lty="dashed")
points(x=1:100, y=res$HOS, col="red", type="l",lty="dotted")
points(x=1:100, y=res$BS, col="blue", type="l")
points(x=1:100, y=res$NOB, col="blue", type="l",lty="dashed")
points(x=1:100, y=res$HOB, col="blue", type="l",lty="dotted")
legend(x=5,y=y.upper.lim*0.83, legend=c("Natural spawners", "NOS", "HOS", "Brood stock removals", "NOB (accounting for mortality prior to spawning)", "HOB"), col=c(rep("red",3),rep("blue",3)), cex=0.7, lty=rep(c("solid", "dashed","dotted"),2), bty="n")
text(x=2, y=y.upper.lim*0.95, labels="(b)")


plot(1:100, res$pHOSeff, type="l", col="red", ylim=c(0,1), ylab="Proportion (pHOSeff or pNOB)", xlab="Generation")
lines(1:100, res$pNOB, col="blue")
legend(x=5,y=0.5, legend=c("pHOSeff","pNOB"), cex=0.7, col=c("red", "blue"), lty=c("solid", "solid"), bty="n")
text(x=2, y=0.93, labels="(c)")

plot(1:100, res$PNI, type="l", col="black", ylim=c(0,1), ylab="PNI",  xlab="Generation")
text(x=2, y=0.95, labels="(d)")

plot(1:100, res$RperS, ylim=c(min(res$RperS)*0.8, max(res$RperS)*1.2), type="l", col="black",  ylab="Recruits/spawner on spawning grounds", xlab="Generation")
text(x=2, y=max(res$RperS)*1.18, labels="(e)")

y.upper.lim<-max(res$ret.nat.preharvest, na.rm=T)#res$c*res$mar.surv
#if(max(res$ret.nat)>res$Seq*1.2){y.upper.lim<-max(res$ret.nat)*1.1}
plot(1:100, res$ret.nat.preharvest, type="l", col="red", ylim=c(0,y.upper.lim), ylab="Recruitment from HOS+NOS", xlab="Generation")
points(x=1:100, res.nogenetics$ret.nat.preharvest, type="l", col="blue")
legend(x=5,y=y.upper.lim*0.8, legend=c("Without genetic impacts",  "With genetic impacts"), cex=0.7, col=c("blue", "red"), lty=c("solid", "solid"), bty="n")
text(x=2, y=y.upper.lim*0.95, labels="(f)")

}

c<-40000; percent.hatch<-0; HR<-0.4; h=sqrt(0.25); w<-sqrt(100)
pdf("TimeseriesLevers18Jan2016.pdf")
res<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0, Theta.hatch=80, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
res.nogenetics<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0, Theta.hatch=100, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever(res, res.nogenetics)

res<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0.5, Theta.hatch=80, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
res.nogenetics<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0.5, Theta.hatch=100, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever(res, res.nogenetics)

res<-run.lever.model(per.mark=0.5,hatchery.size=0.3, sel=0.5, Theta.hatch=80, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
res.nogenetics<-run.lever.model(per.mark=0.5,hatchery.size=0.3, sel=0.5, Theta.hatch=100, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever(res, res.nogenetics)
dev.off()

res<-run.lever.model(per.mark=0.5,hatchery.size=0.3, sel=0, Theta.hatch=80, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
res.nogenetics<-run.lever.model(per.mark=1,hatchery.size=0.03, sel=1, Theta.hatch=100, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever(res, res.nogenetics)
#per.mark=1;hatchery.size=0.1; sel=1; Theta.hatch=80; c=c; percent.hatch=percent.hatch; HR=HR; h=h; w=w

PNIm<-NA; pHOSm<-NA; pNOBm<-NA; RperSm<-NA; Ret.natm<-NA; PNIs<-NA; pHOSs<-NA; pNOBs<-NA; RperSs<-NA; Ret.nats<-NA; PNIh<-NA; pHOSh<-NA; pNOBh<-NA; BSh<-NA; RperSh<-NA; Ret.nath<-NA;  
for(i in 1:100){
  m<-run.lever.model(per.mark=i*0.01, hatchery.size=0.3, sel=0, Theta.hatch=80, c=c, sex.ratio=0.5, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
  m.BS<-m$BS
  m.hatchery.size<-m$hatchery.size
  m.sel<-m$sel
  PNIm[i]<-m$PNI[100]
  pHOSm[i]<-m$pHOSeff[100]
  pNOBm[i]<-m$pNOB[100]
  RperSm[i]<-m$RperS[100]
  Ret.natm[i]<-m$ret.nat.preharvest[100]
  s<-run.lever.model(per.mark=0.1, hatchery.size=0.3, sel=i*0.01, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
  s.BS<-s$BS
  s.hatchery.size<-s$hatchery.size
  s.per.mark<-s$per.mark
  PNIs[i]<-s$PNI[100]
  pHOSs[i]<-s$pHOSeff[100]
  pNOBs[i]<-s$pNOB[100]
  RperSs[i]<-s$RperS[100]  
  Ret.nats[i]<-s$ret.nat.preharvest[100]  
  h<-run.lever.model(per.mark=0.1, hatchery.size=i*0.003, sel=0, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
  h.per.mark<-h$per.mark
  h.sel<-h$sel
  BSh[i]<-h$BS[100]
  PNIh[i]<-h$PNI[100]
  pHOSh[i]<-h$pHOSeff[100]
  pNOBh[i]<-h$pNOB[100]
  RperSh[i]<-h$RperS[100]
  Ret.nath[i]<-h$ret.nat.preharvest[100]
  }
 
pdf("EquilibriumLevers29Dec2016.pdf")
par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(x=(1:100)*0.01, y=pHOSm, ylim=c(0,1), ylab="pHOSeff or pNOB at equilibrium", xlab="", type="l")
mtext("% Marked",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=0.97, label="(a)")
lines(x=(1:100)*0.01, y=pNOBm, lty="dashed")
legend(x=0.1, y=0.5, legend=c("pHOSeff", "pNOB"), lty=c("solid", "dashed"), bty="n")
plot(x=(1:100)*0.01, y=PNIm, ylim=c(0,1), ylab="PNI at equilbrium", xlab="", type="l")
mtext("% Marked",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=0.97, label="(b)")
plot(x=(1:100)*0.01, y=RperSm, ylim=c(0,2.5), ylab="Recruits/spawner at equilbrium", xlab="", type="l")
mtext("% Marked",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=2.5*0.97, label="(c)")
plot(x=(1:100)*0.01, y=Ret.natm, ylim=c(0,max(Ret.natm)), ylab="Recruits from HOS+HOS at equilbrium", xlab="", type="l")
mtext("% Marked",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=max(Ret.natm)*0.97, label="(d)")
title<-paste("Levers held constant: Brood stock=", round(m.BS[100],0)," (",m.hatchery.size*100,"% of ave natural recruitment); Removal of marked fish=",m.sel*100,"%", sep="")
mtext(title, side=1, line=2.8, at=-0.2, cex=0.8)
#mtext("test", side=1, line=2.8, at=-0.2, cex=0.8)


par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(x=(1:100)*0.01, y=pHOSs, ylim=c(0,1), ylab="pHOSeff or pNOB at equilibrium", xlab="", type="l")
mtext("% Removal of marked fish",1,line=1.8,at=0.5, cex=0.8)
lines(x=(1:100)*0.01, y=pNOBs, lty="dashed")
text(x=0.02, y=0.97, label="(a)")
legend(x=0.1, y=0.5, legend=c("pHOSeff", "pNOB"), lty=c("solid", "dashed"), bty="n")
plot(x=(1:100)*0.01, y=PNIs, ylim=c(0,1), ylab="PNI at equilbrium", xlab="", type="l")
mtext("% Removal of marked fish",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=0.97, label="(b)")
plot(x=(1:100)*0.01, y=RperSs, ylim=c(0,2.5), ylab="Recruits/spawner at equilbrium", xlab="", type="l")
mtext("% Removal of marked fish",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=2.5*0.97, label="(c)")
title<-paste("Levers held constant: Mark rate=", s.per.mark*100, "%;  Brood stock=", round(s.BS[100],0)," (",s.hatchery.size*100,"% of ave natural recruitment)", sep="")
mtext(title, side=1, line=2.8, at=1.1, cex=0.8)
plot(x=(1:100)*0.01, y=Ret.nats, ylim=c(0,max(Ret.nats)), ylab="Recruits from HOS+NOS at equilbrium", xlab="", type="l")
mtext("% Removal of marked fish",1,line=1.8,at=0.5, cex=0.8)
text(x=0.02, y=max(Ret.nats)*0.97, label="(d)")

par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(x=BSh, y=pHOSh, ylim=c(0,1), ylab="pHOSeff or pNOB at equilibrium", xlab="", type="l")#, axes=FALSE)
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(2,at=(0:5)*0.2, labels=(0:5)*0.2,line=0)
lines(x=BSh, y=pNOBh, lty="dashed")
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.2,1.2), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.2,1.2), labels=NA, line=0)
text(x=30, y=0.97, label="(a)")
legend(x=0.07, y=0.6, legend=c("pHOSeff", "pNOB"), lty=c("solid", "dashed"), bty="n")

plot(x=BSh, y=PNIh, ylim=c(0,1), ylab="PNI at equilbrium", xlab="", type="l")#, axes=FALSE)
text(x=30, y=0.97, label="(b)")
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(2,at=(0:5)*0.2, labels=(0:5)*0.2,line=0)
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.2,1.2), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.2,1.2), labels=NA, line=0)

plot(x=BSh, y=RperSh, ylim=c(0,2.5), ylab="Recruits/spawner at equilbrium", xlab="", type="l")
text(x=50, y=2.5*0.97, label="(c)")
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(2,at=(0:5)*0.5, labels=(0:5)*0.5,line=0)
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.5,3), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.5,3), labels=NA, line=0)
title<-paste("Levers held constant: Mark rate=", h.per.mark*100, "%;  Removal of marked fish=",h.sel*100,"%", sep="")
mtext(title, side=1, line=2.8, at=max(BSh, na.rm=T), cex=0.8)

plot(x=BSh, y=Ret.nath, ylim=c(0,max(Ret.nath)), ylab="Recruits from HOS+NOS at equilbrium", xlab="", type="l")
text(x=20, y=max(Ret.nath)*0.97, label="(d)")
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(2,at=(0:6)*500,labels=(0:6)*500, line=0)
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=max(BSh, na.rm=T)/2, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.5,3), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.5,3), labels=NA, line=0)
dev.off()

#plot(x=(1:100)*0.003, y=pHOSh, ylim=c(0,1), ylab="pHOS or pNOB at equilibrium", xlab="", type="l", axes=FALSE, bty="o")
#lines(x=(1:100)*0.003, y=pNOBh, lty="dashed")
#axis(1,(0:10)*0.03,labels=0:10*0.03,line=0)
#mtext("Brood stock",1,line=2,at=0.15, cex=0.8)
#axis(2,at=(0:5)*0.2, labels=(0:5)*0.2,line=0)
#axis(1,(0:10)*0.03,labels=0:10*0.03,line=3)
#mtext("% Ave natural recruitment",1,line=5,at=0.15, cex=0.8)
#text(x=0.006, y=0.97, label="(a)")
#legend(x=0.07, y=0.5, legend=c("pHOS", "pNOB"), lty=c("solid", "dashed"), bty="n")
#title<-paste("Levers held constant: Mark rate=", h.per.mark*100, "%;  Removal of marked fish=",h.sel*100,"%", sep="")
#mtext(title, side=3, line=0.5, at=0.35, cex=0.8)
#plot(x=(1:100)*0.003, y=PNIh, ylim=c(0,1), ylab="PNI at equilbrium", xlab="Size of hatchery (% of ave natural returns)", type="l")
#text(x=0.006, y=0.97, label="(b)")
#plot(x=(1:100)*0.003, y=RperSh, ylim=c(0,max(RperSh)), ylab="Recruits/spawner at equilbrium", xlab="Size of hatchery (% of ave natural returns)", type="l")
#text(x=0.007, y=max(RperSh)*0.97, label="(c)")
#plot(x=(1:100)*0.003, y=Ret.nath, ylim=c(0,max(Ret.nath)), ylab="Returns from HOS+NOS at equilbrium", xlab="Size of hatchery (% of ave natural returns)", type="l")
#text(x=0.002, y=max(Ret.nath)*0.97, label="(d)")

