#Code to estimate effect of 3 hatchery levers: harchery size, marking, and selective harvest
#Last updated 13 Feb. 2017
#rm(list = ls())

source("R/FuncDefs.R")
source("R/SeqFormulationBH.R")
plot.UnivariateSA<-FALSE
plot.timeseries<-FALSE

run.lever.model<-function(per.mark,hatchery.size, sel, Theta.hatch, c, percent.hatch, HR, h, w, mar.surv, RS, mar.surv.hatch, sex.ratio){
#per.mark=0.75; hatchery.size=0.45; sel=0; Theta.hatch=80; c=40000; percent.hatch=0; HR=0.4; h=sqrt(0.25); w=sqrt(100); mar.surv=0.02; RS=0.8; mar.surv.hatch=0.005
#per.mark=1.0; hatchery.size=80*0.005; sel=90*0.01; Theta.hatch=80;mar.surv=0.02; RS=0.8; mar.surv.hatch=0.0024
#per.mark=1.0; hatchery.size=95*0.005; sel=2*0.01; Theta.hatch=80;mar.surv=0.02; RS=0.8; mar.surv.hatch=0.0024
  
#Levers
#per.mark<-1.0#Percentage of hatchery origin fish that are marked
#hatchery.size<-0.1#Hatchery size as a ppn of equilibrium capacity (Seq) of natural spawners
#sel<-0#Selective removals/harvest
  
#Biological Parameters
p<-175#Beverton-Holt productivity parameter
#c<-400000#Beverton-Holt capacity parameter (Rmax)
#mar.surv<-0.02# Marine survival of natural-origin fish
#RS<-0.8#Relative survival/fecundity of hatchery vs. natural-origin spawners in natural environment in BH model

#Hatchery parameters
bs.surv<-0.80#survival from capture to spawning for broodstock
fec<-4900#fecundity
release.surv<-0.88# combined survival of hatchery fish from egg-fry and fry-smolt
#mar.surv.hatch<-0.0024 #Marine survival of hatchery-origin fish
#percent.hatch<-0# 1-percent of hatchery fish that return to natural spawning grounds. Look at a range of values
rel.loss<-0.5#Ppn of fitness loss at smolt stage (compared with adult stage)
#HR<-0.4#Harvest rate

#Genetic parameters
sig<-sqrt(10)
#w<-sqrt(100)
#h<-sqrt(0.5)
Theta.nat<-100
#Theta.hatch<-80

#Set up vectors
P.nat<-NA#Mean phenotypic value in natural environment
P.hatch<-NA#Mean phenotypic value in the hatchery
HOS<-NA
HOR<-NA#Hatchery-origin returns
HOR.unmark<-NA
HOR.mark<-NA
HOR.rem<-NA
PNI<-NA
pHOS<-NA
pHOSEff<-NA
pNOB<-NA
NOS<-NA
HOB<-NA
NOB<-NA
ret.nat<-NA#Returns to the spawning grounds AFTER harvest
ret.nat.preharvest<-NA#Returns to the spawning grounds PRIOR TO harvest
ret.hatch<-NA#Returns to the hatchery AFTER harvest
ret.hatch.preharvest<-NA#Returns to the hatchery PRIOR TO harvest
BS<-NA#Brood stock
BS.mark<-NA
ppn.unmarked.spawners<-NA
Sp.nat<-NA#Total spawners in the natural environment
Sm.nat<-NA#Total smolts from the natural environment
Sm.hatch<-NA#Total smolts from the hatchery
fit.smolt<-NA#fitness mulitiplier at smolt stage
fit.adult<-NA#fitness mulitiplier at adult stage
RperS<-NA#Recruits from natural spawners  prior to harvest /(natural spawners=HOS+NOS)
catch<-NA
ext<-0

#Initial parameters
P.nat[1]<-Theta.nat
P.hatch[1]<-Theta.nat
Seq<-(c*mar.surv*(1-HR))*((p*mar.surv*(1-HR))-1)/(p*mar.surv*(1-HR)) #Assuming BH formulation #2 from Hilborn and Walters (1992), spawner abundances at equilibrium adult recruitment# See "SeqFormulation.r"
if(mar.surv.const==TRUE){Seq<-(c*0.02*(1-HR))*((p*0.02*(1-HR))-1)/(p*0.02*(1-HR)) #Assuming BH formulation #2 from Hilborn and Walters (1992), spawner abundances at equilibrium adult recruitment# See "SeqFormulation.r"
  }
Req<-Seq

#Dynamics
#Brood stock estimation in year prior to the first generation, required to estimate hatchery returns and HOS
#pred.hatch.recruits<-solver.HS(hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$fit
#NR.check<-ratio.func(pred.hatch.recruits, hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$Ret#just to CHECK what solver used for total natural returns
#HOS.check<-ratio.func(pred.hatch.recruits, hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$HOS[100]#just to CHECK what solver used for total natural returns
#NOS.check<-ratio.func(pred.hatch.recruits, hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$NOS[100]#just to CHECK what solver used for total natural returns
#BS.set<-pred.hatch.recruits/(fec*release.surv*mar.surv.hatch*(1HR))#The broodstock required to acheive predicted hatchery returns

BS.set<-hatchery.size*Req#Hatchery size as a function of longterm equilibrium recrutimetn in the absence of the hatchery
ret.hatch.minus1<-BS.set*bs.surv*sex.ratio*fec*release.surv*mar.surv.hatch*(1-HR)# Assuming no fitness effects here

HOR[1]<-ret.hatch.minus1*(1-percent.hatch)#Hatchery-origin returns to natural spawning grounds
HOR.mark[1]<-HOR[1]*per.mark
HOR.unmark[1]<-HOR[1]*(1-per.mark)#Number of HOR that are unmarked
HOR.rem[1]<-HOR[1]*per.mark*sel#Number of marked HOR that are selectively removed

if(BS.set>(0.33*(HOR[1]+Seq-HOR.rem[1]))){BS[1]<-0.33*(HOR[1]+Seq-HOR.rem[1])}
if(BS.set<=(0.33*(HOR[1]+Seq-HOR.rem[1]))){BS[1]<-BS.set}


ppn.unmarked.spawners[1]<-(HOR.unmark[1]+Seq)/(HOR.unmark[1]+HOR.mark[1]*(1-sel)+Seq)
if(ppn.unmarked.spawners[1]*BS[1]*3>=BS[1]){# If brood can be collected from unmarked fish with a limit of handling effort to sample size of 3XBS 

  pNOB[1]<-Seq/(Seq+HOR.unmark[1])#(Seq)/(Seq+HOR.unmark[1])#wild spawners/wild spawners + unmarked hatchery spawners
  NOS[1]<-Seq-BS[1]*(Seq/((HOR.unmark[1]+Seq)))#only a portion of BS since some come from unmarked hatchery fish
  if(NOS[1]<0)NOS[1]<-0
  HOS[1]<-HOR[1]-HOR.rem[1]-BS[1]*(HOR.unmark[1]/((HOR.unmark[1]+Seq)))#also subtract  brood take in proportion that HOR provide marked fish to returns
  
}#End of if(ppn.unmarked.spawners[1]*BS[1]*3>BS[1]){

if(ppn.unmarked.spawners[1]*BS[1]*3<BS[1]){#If brood cannot be collected from unmarked fish with a limit on handling effort to sample size of  3XBS
  BS.underhandling<-ppn.unmarked.spawners[1]*BS[1]*3#Number of unmarked BS taken within handling limits (<BS.set)
  ppnBS.overhandling<-1-BS.underhandling/BS[1]# Proprtion of BS taken after handling limit exceeded. BS taken in the ratio of HOR and NOR occuring
  ppnBS.underhandling<-BS.underhandling/BS[1]# Proportion of BS taken before handling limit exceeded. BS taken in the ration of unmarked HOR and NOR occuring
  
  BS.NO.underhandling<-BS[1]*(Seq/((HOR.unmark[1]+Seq)))*ppnBS.underhandling
  BS.NO.overhandling<-BS[1]*(Seq/(HOR[1]-HOR.rem[1]+Seq))*ppnBS.overhandling
  BS.HO.underhandling<-BS[1]*(HOR.unmark[1]/((HOR.unmark[1]+Seq)))*ppnBS.underhandling
  BS.HO.overhandling<-BS[1]*((HOR[1]-HOR.rem[1])/(HOR[1]-HOR.rem[1]+Seq))*ppnBS.overhandling
  
  NOS[1]<-Seq -  BS.NO.underhandling -  BS.NO.overhandling
  #if(NOS[1]<0)NOS[1]<-0

  pNOB[1]<-Seq/(Seq+HOR.unmark[1])*ppnBS.underhandling + Seq/(Seq+HOR[1]-HOR.rem[1])*ppnBS.overhandling #combing pNOB for ppn fish prior to and after limit exceeded
  HOS[1]<-HOR[1]-HOR.rem[1]- BS.HO.underhandling - BS.HO.overhandling# subtract  brood take
  
}#End of if(ppn.unmarked.spawners[1]*BS[1]*3<BS[1]){


NOB[1]<-BS[1]*pNOB[1]
HOB[1]<-BS[1]*(1-pNOB[1])
pHOS[1]<-HOS[1]/(HOS[1]+NOS[1])
pHOSEff[1]<-HOS[1]*RS/(HOS[1]*RS+NOS[1])
PNI[1]<-pNOB[1]/(pNOB[1]+pHOS[1])

Sp.nat[1]<-NOS[1]+HOS[1]

fit.smolt[1]<-fit.lifestage(P.nat[1], Theta.nat, w, sig, rel.loss)
Sm.nat[1]<-BH(HOS[1], NOS[1], RS, p, c)*fit.smolt[1]
Sm.hatch[1]<-Hatch.sm(BS[1]*bs.surv, fec, sex.ratio, release.surv)

fit.adult[1]<-fit.lifestage(P.nat[1], Theta.nat, w, sig, 1-rel.loss)
ret.nat[1]<-Sm.nat[1]*mar.surv*fit.adult[1]*(1-HR)#Returns to spawning grounds AFTER harvest
ret.hatch.preharvest[1]<-Sm.hatch[1]*mar.surv.hatch#Returns to hatchery AFTER harvest
ret.hatch[1]<-Sm.hatch[1]*mar.surv.hatch*(1-HR)#Returns to hatchery AFTER harvest
ret.nat.preharvest[1]<-Sm.nat[1]*mar.surv*fit.adult[1]
RperS[1]<-ret.nat.preharvest[1]/(HOS[1]+NOS[1])
catch[1]<-(ret.hatch.preharvest[1]+ret.nat.preharvest[1])*HR + ret.hatch[1]*(1-percent.hatch)*(per.mark)*sel#Includes harvest + fish selectively removed from river after harvest

for (i in 2:100){#for i generations)
  if(ext==0){#If population NOT extirpated (natural population extirpated because fitness impacts and 100% hatchery fish removed prior to spawning)
  HOR[i]<-ret.hatch[i-1]*(1-percent.hatch)#Hatchery-origin returns to natural spawning grounds
  HOR.unmark[i]<-HOR[i]*(1-per.mark)#Number of HOR that are unmarked
  HOR.mark[i]<-HOR[i]*(per.mark)#Number of HOR that are unmarked
  HOR.rem[i]<-HOR.mark[i]*sel#Number of marked HOR that are selectively removed

  if(ret.nat[i-1]>0){
    
    if(BS.set>(0.33*(ret.nat[i-1]+HOR[i]-HOR.rem[i]))){BS[i]<-0.33*(ret.nat[i-1]+HOR[i]-HOR.rem[i])}
    if(BS.set<=(0.33*(ret.nat[i-1]+HOR[i]-HOR.rem[i]))){BS[i]<-BS.set}
    if(i==100){if(BS[i]==BS.set){BS.mark<-0}; if(BS[i]!=BS.set){BS.mark<-1}}
    ppn.unmarked.spawners[i]<-(HOR.unmark[i]+ret.nat[i-1])/(HOR.unmark[i]+HOR.mark[i]*(1-sel)+ret.nat[i-1])
    if(ppn.unmarked.spawners[i]*BS[i]*3>=BS[i]){# If brood can be collected from unmarked fish with a limit of handling effort to sample size of 3XBS 
      NOS[i]<-ret.nat[i-1]-BS[i]*(ret.nat[i-1]/(HOR.unmark[i]+ret.nat[i-1]))#substract the portion of BS that is natural-origin
      if(NOS[i]<0)cat(NOS[i])#NOS[i]<-0
      HOS[i]<-HOR.unmark[i]+HOR.mark[i]*(1-sel)-BS[i]*(HOR.unmark[i]/(HOR.unmark[i]+ret.nat[i-1]))#substract the portion of BS that is hatchery-origin
      pNOB[i]<-(ret.nat[i-1])/(ret.nat[i-1]+HOR.unmark[i])#HOR.unmark[i])#wild spawners/wild spawners + unmarked hatchery spawners. BS taken entirely from spawning grounds
    }#End of if(ppn.unmarked.spawners[i]*BS[i]*3>BS[i]){
  
    if(ppn.unmarked.spawners[i]*BS[i]*3<BS[i]){#If brood cannot be collected from unmarked fish with a limit on handling effort to sample size of  3XBS
      #cat("handling limit reached in generation, with hatchery size, %mark, sel", i, hatchery.size, per.mark, sel )
      BS.underhandling<-ppn.unmarked.spawners[i]*BS[i]*3#Number of unmarked BS taken within handling limits (<BS[i])
      ppnBS.overhandling<-1-BS.underhandling/BS[i]# Proprtion of BS taken after handling limit exceeded. BS taken in the ratio of HOR and NOR occuring
      ppnBS.underhandling<-BS.underhandling/BS[i]# Proportion of BS taken before handling limit exceeded. BS taken in the ration of unmarked HOR and NOR occuring
    
      BS.NO.underhandling<-BS[i]*(ret.nat[i-1]/((HOR.unmark[i]+ret.nat[i-1])))*ppnBS.underhandling
      BS.NO.overhandling<-BS[i]*(ret.nat[i-1]/(HOR[i]-HOR.rem[i]+ret.nat[i-1]))*ppnBS.overhandling
      BS.HO.underhandling<-BS[i]*(HOR.unmark[i]/((HOR.unmark[i]+ret.nat[i-1])))*ppnBS.underhandling
      BS.HO.overhandling<-BS[i]*((HOR[i]-HOR.rem[i])/(HOR[i]-HOR.rem[i]+ret.nat[i-1]))*ppnBS.overhandling
    
      NOS[i]<-ret.nat[i-1]- BS.NO.underhandling -  BS.NO.overhandling
      #if(NOS[i]<0)NOS[i]<-0
      pNOB[i]<-(ret.nat[i-1])/(ret.nat[i-1]+HOR.unmark[i])*ppnBS.underhandling + (ret.nat[i-1])/(ret.nat[i-1]+HOR[i]-HOR.rem[i])*ppnBS.overhandling#combing pNOB for unmarked component and marked component (which as pNOB=0)  
      HOS[i]<-HOR[i]-HOR.rem[i]- BS.HO.underhandling - BS.HO.overhandling#subtract  brood take  
    }#End of if(ppn.unmarked.spawners[i]*BS[i]*3<BS[i]){

  }#End of if(ret.nat[i-1]>0)
  
  if(ret.nat[i-1]<=0){ppn.unmarked.spawners[i]<-0; NOS[i]<-0; HOS[i]<-HOR[i]-HOR.rem[i]; pNOB[i]<-0; ext<-1}#And BS[i]=BS[i] (not 0.33 of ret.nat+HOR)
  
  NOB[i]<-BS[i]*pNOB[i]
  HOB[i]<-BS[i]*(1-pNOB[i])
  pHOS[i]<-HOS[i]/(HOS[i]+NOS[i])
  pHOSEff[i]<-HOS[i]*RS/(HOS[i]*RS+NOS[i])
  PNI[i]<-pNOB[i]/(pNOB[i]+pHOSEff[i])#pNOB[i]/(pNOB[i]+pHOSEff[i])
  Sp.nat[i]<-NOS[i]+HOS[i]
  
  P.nat[i]<-Pnat(pHOSEff[i-1], P.nat[i-1], w, sig, Theta.nat, h, P.hatch[i-1])
  P.hatch[i]<-Phatch(pNOB[i-1], P.hatch[i-1], w, sig, Theta.hatch, h, P.nat[i-1])
  fit.smolt[i]<-fit.lifestage(P.nat[i], Theta.nat, w, sig, rel.loss)
  Sm.nat[i]<-BH(HOS[i], NOS[i], RS, p, c)*fit.smolt[i]
  Sm.hatch[i]<-Hatch.sm(BS[i]*bs.surv, fec, sex.ratio, release.surv)
  
  fit.adult[i]<-fit.lifestage(P.nat[i], Theta.nat, w, sig, 1-rel.loss)
  ret.nat[i]<-Sm.nat[i]*mar.surv*fit.adult[i]*(1-HR)
  ret.hatch[i]<-Sm.hatch[i]*mar.surv.hatch*(1-HR)
  ret.hatch.preharvest[i]<-Sm.hatch[i]*mar.surv.hatch
  ret.nat.preharvest[i]<-Sm.nat[i]*mar.surv*fit.adult[i]
  RperS[i]<-ret.nat.preharvest[i]/(HOS[i]+NOS[i])
  catch[i]<-(ret.hatch.preharvest[i]+ret.nat.preharvest[i])*HR + ret.hatch[i]*(1-percent.hatch)*(per.mark)*sel#Includes harvest + fish selectively removed from river after harvest
  }#End of if(ext==0) If population not extirpated

if(ext==1){
  NOS[i]<-0; HOS[i]<-HOR[i]-HOR.rem[i]; pNOB[i]<-0; HOB[i]<-0; pHOS[i]<-1; pHOSEff[i]<-1; PNI[i]<-0; 
  Sm.hatch[i]<-0;#Hatch.sm(BS.set*bs.surv, fec, release.surv)
  ret.hatch[i]<-0#Sm.hatch[i]*mar.surv.hatch*(1-HR)
  ret.hatch.preharvest[i]<-0#Sm.hatch[i]*mar.surv.hatch
  ret.nat.preharvest[i]<-0
  ret.nat[i]<-0
  RperS[i]<-0
  catch[i]<-0#(ret.hatch.preharvest[i]+ret.nat.preharvest[i])*HR + ret.hatch[i]*(1-percent.hatch)*(per.mark)*sel
}

}#End if i generations

  return(list(fit.smolt=fit.smolt, P.nat=P.nat, P.hatch=P.hatch, Theta.hatch=Theta.hatch, Theta.nat=Theta.nat, Sp.nat=Sp.nat, ret.nat=ret.nat, ret.nat.preharvest=ret.nat.preharvest, ret.hatch.preharvest=ret.hatch.preharvest, catch=catch, BS=rep(BS.set,100), Seq=Seq, NOS=NOS, NOB=NOB, HOS=HOS, HOB=HOB, pNOB=pNOB, pHOS=pHOS, PNI=PNI, per.mark=per.mark, hatchery.size=hatchery.size, RperS=RperS, sel=sel, mar.surv=mar.surv, c=c, pHOSeff=pHOSEff, fit.adult=fit.adult, BS.mark=BS.mark))
}#End of run.lever.model

panel.plots.lever<-function(){

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
legend(x=5,y=y.upper.lim*0.83, legend=c("Natural spawners", "NOS", "HOS", "Brood stock removals", "NOB", "HOB"), col=c(rep("red",3),rep("blue",3)), cex=0.7, lty=rep(c("solid", "dashed","dotted"),2), bty="n")
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

if(plot.timeseries==TRUE){
pdf("TimeseriesLevers19Jan2017.pdf")
res<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)#per.mark,hatchery.size, sel, Theta.hatch, c, percent.hatch, HR, h, w
res.nogenetics<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0, Theta.hatch=100, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever()

res<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0.5, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
res.nogenetics<-run.lever.model(per.mark=1.0,hatchery.size=0.3, sel=0.5, Theta.hatch=100, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever()

res<-run.lever.model(per.mark=0.5,hatchery.size=0.3, sel=0.5, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
res.nogenetics<-run.lever.model(per.mark=0.5,hatchery.size=0.3, sel=0.5, Theta.hatch=100, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
panel.plots.lever()
dev.off()
}

#per.mark=1;hatchery.size=0.1; sel=1; Theta.hatch=80
if(plot.UnivariateSA==TRUE){

PNIm<-NA; pHOSm<-NA; pNOBm<-NA; RperSm<-NA; Ret.natm<-NA; PNIs<-NA; pHOSs<-NA; pNOBs<-NA; RperSs<-NA; Ret.nats<-NA; PNIh<-NA; pHOSh<-NA; pNOBh<-NA; BSh<-NA; RperSh<-NA; Ret.nath<-NA;  
for(i in 1:100){
  m<-run.lever.model(per.mark=i*0.01, hatchery.size=0.3, sel=0, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
  m.BS<-m$BS
  m.hatchery.size<-m$hatchery.size
  m.sel<-m$sel
  PNIm[i]<-m$PNI[100]
  pHOSm[i]<-m$pHOSeff[100]
  pNOBm[i]<-m$pNOB[100]
  RperSm[i]<-m$RperS[100]
  Ret.natm[i]<-m$ret.nat.preharvest[100]
  s<-run.lever.model(per.mark=0.1, hatchery.size=0.3, sel=i*0.01, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
  s.BS<-s$BS
  s.hatchery.size<-s$hatchery.size
  s.per.mark<-s$per.mark
  PNIs[i]<-s$PNI[100]
  pHOSs[i]<-s$pHOS[100]
  pNOBs[i]<-s$pNOB[100]
  RperSs[i]<-s$RperS[100]  
  Ret.nats[i]<-s$ret.nat.preharvest[100]  
  h<-run.lever.model(per.mark=0.1, hatchery.size=i*0.003, sel=0, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024)
  h.per.mark<-h$per.mark
  h.sel<-h$sel
  BSh[i]<-h$BS[100]
  PNIh[i]<-h$PNI[100]
  pHOSh[i]<-h$pHOS[100]
  pNOBh[i]<-h$pNOB[100]
  RperSh[i]<-h$RperS[100]
  Ret.nath[i]<-h$ret.nat.preharvest[100]
}
 
pdf("EquilibriumLevers22Dec2016.pdf")
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
plot(x=(1:100)*0.01, y=pHOSs, ylim=c(0,1), ylab="pHOS or pNOB at equilibrium", xlab="", type="l")
mtext("% Removal of marked fish",1,line=1.8,at=0.5, cex=0.8)
lines(x=(1:100)*0.01, y=pNOBs, lty="dashed")
text(x=0.02, y=0.97, label="(a)")
legend(x=0.1, y=0.5, legend=c("pHOS", "pNOB"), lty=c("solid", "dashed"), bty="n")
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
plot(x=BSh, y=pHOSh, ylim=c(0,1), ylab="pHOS or pNOB at equilibrium", xlab="", type="l")#, axes=FALSE)
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=175, cex=0.8)
#axis(2,at=(0:5)*0.2, labels=(0:5)*0.2,line=0)
lines(x=BSh, y=pNOBh, lty="dashed")
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=175, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.2,1.2), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.2,1.2), labels=NA, line=0)
text(x=30, y=0.97, label="(a)")
legend(x=0.07, y=0.6, legend=c("pHOS", "pNOB"), lty=c("solid", "dashed"), bty="n")

plot(x=BSh, y=PNIh, ylim=c(0,1), ylab="PNI at equilbrium", xlab="", type="l")#, axes=FALSE)
text(x=30, y=0.97, label="(b)")
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=175, cex=0.8)
#axis(2,at=(0:5)*0.2, labels=(0:5)*0.2,line=0)
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=175, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.2,1.2), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.2,1.2), labels=NA, line=0)

plot(x=BSh, y=RperSh, ylim=c(0,2.5), ylab="Recruits/spawner at equilbrium", xlab="", type="l")
text(x=50, y=2.5*0.97, label="(c)")
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=175, cex=0.8)
#axis(2,at=(0:5)*0.5, labels=(0:5)*0.5,line=0)
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=175, cex=0.8)
#axis(1, c(-100,360), labels=NA, line=0); axis(2, c(-0.5,3), labels=NA, line=0); axis(3, c(-100,360), labels=NA, line=0); axis(4, c(-0.5,3), labels=NA, line=0)
title<-paste("Levers held constant: Mark rate=", h.per.mark*100, "%;  Removal of marked fish=",h.sel*100,"%", sep="")
mtext(title, side=1, line=2.8, at=400, cex=0.8)

plot(x=BSh, y=Ret.nath, ylim=c(0,max(Ret.nath)), ylab="Recruits from HOS+NOS at equilbrium", xlab="", type="l")
text(x=20, y=max(Ret.nath)*0.97, label="(d)")
#axis(1,(0:7)*50,labels=(0:7)*50,line=0)
mtext("Brood stock",1,line=1.8,at=175, cex=0.8)
#axis(2,at=(0:6)*500,labels=(0:6)*500, line=0)
axis(3,at=c(0,BSh[33], BSh[67], BSh[100]),labels=0:3*0.1,line=0, cex.axis=0.8, padj=1)
mtext("% Ave natural recruitment",3,line=1.2,at=175, cex=0.8)
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

}#End of if plot.UnivariateSA==TRUE