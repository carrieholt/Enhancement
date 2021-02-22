# Population model to evalaute genetic impact of hatcheries
# Date last revised: 6 Dec. 2016
# Created by: Carrie Holt
# Requires FuncDefs.r

source("FuncDefs.r")

Run.model<-function(Theta.hatch, RS, pNOB, pHOS){

#Input parameters: CHANGE TO MIKE's
#pNOB<-0.67
#pHOS<-0.33
sig<-sqrt(10)
w<-sqrt(100)
h<-sqrt(0.5)
Theta.nat<-100
#Theta.hatch<-80
p<-2.5#BH parameters
c<-500
#RS<-0.8#Relative survival of hatchery vs. natural spawners in BH model
#Hatchery parameters
fec<-100#fecundity
release.surv<-0.1#survival of hatchery fish from egg to smolt
mar.surv<-0.9# Marine survival of natural-origin fish
mar.surv.hatch<-mar.surv*0.9 #Marine survival of hatchery-origin fish
percent.hatch<-0.8# percent of hatchery fish that return to natural spawning grounds
rel.loss<-0.5#Ppn of fitness loss at smolt stage (compared with adult stage)

#Set up vectors
P.nat<-NA#Mean phenotypic value in natural environment
P.hatch<-NA#Mean phenotypic value in the hatchery
HOS<-NA
NOS<-NA
HOB<-NA
NOB<-NA
ret.nat<-NA#Returns to the spawning grounds
ret.hatch<-NA#Returns to the hatchery
BS<-NA#Brood stock
Sp.nat<-NA#Total spawners in the natural environment
Sm.nat<-NA#Total smolts from the natural environment
Sm.hatch<-NA#Total smolts from the hatchery
fit.smolt<-NA#fitness mulitiplier at smolt stage
fit.adult<-NA#fitness mulitiplier at adult stage

#Initial parameters
P.nat[1]<-Theta.nat
P.hatch[1]<-Theta.nat
Seq<-c*(p-1)/p #Assuming BH formulation #2 from Hilborn and Walters (1992)

#Brood stock estimation in year prior to the first generation, required to estimate hatchery returns and HOS
#BS is optimized to acheive pHOS
BS.minus1<-as.numeric(solver.BS (pNOB, pHOS, fec, release.surv, mar.surv.hatch, Seq, percent.hatch))#Broodstock in generation prior to 1st generation, assuming ret.nat[gen-1]=Seq
ret.hatch.minus1<-BS.minus1*fec*release.surv*mar.surv.hatch# Assuming no fitness effects here
HOS[1]<-ret.hatch.minus1*(1-percent.hatch)

ret.nat[1]<-Seq
BS[1]<-as.numeric(solver.BS(pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat[1], percent.hatch))# Assuming BS[1] = BS[gen -1] 
  #BS that can acheive pHOS given percent.hatch and the natural returns from the previous generation
NOB[1]<-BS[1]*pNOB
if(NOB[1]>ret.nat[1]*0.5){NOB[1]<-0.5*ret.nat; print("NOB>0.5*ret.nat in year 1")}
HOB[1]<-BS[1]-NOB[1]
NOS[1]<-ret.nat[1]-NOB[1]
Sp.nat[1]<-NOS[1]+HOS[1]

fit.smolt[1]<-fit.lifestage(P.nat[1], Theta.nat, w, sig, rel.loss)
Sm.nat[1]<-BH(HOS[1], NOS[1], RS, p, c)*fit.smolt[1]
Sm.hatch[1]<-Hatch.sm(BS[1], fec, release.surv)

fit.adult[1]<-fit.lifestage(P.nat[1], Theta.nat, w, sig, 1-rel.loss)
ret.nat[1]<-Sm.nat[1]*mar.surv*fit.adult[1]
ret.hatch[1]<-Sm.hatch[1]*mar.surv.hatch

for (i in 2:100){#for i generations
  BS[i]<-as.numeric(solver.BS(pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat[i-1], percent.hatch))
  NOB[i]<-BS[i]*pNOB
  if(NOB[i]>ret.nat[i-1]*0.5){NOB[i]<-0.5*ret.nat; cat("NOB>0.5*ret.nat in year",i)}
  HOB[i]<-BS[i]-NOB[i]
  NOS[i]<-ret.nat[i-1]-NOB[i]
  HOS[i]<-ret.hatch[i-1]*(1-percent.hatch)
  Sp.nat[i]<-NOS[i]+HOS[i]
  
  P.nat[i]<-Pnat(pHOS, P.nat[i-1], w, sig, Theta.nat, h, P.hatch[i-1])
  P.hatch[i]<-Phatch(pNOB, P.hatch[i-1], w, sig, Theta.hatch, h, P.nat[i-1])
  fit.smolt[i]<-fit.lifestage(P.nat[i], Theta.nat, w, sig, rel.loss)
  Sm.nat[i]<-BH(HOS[i], NOS[i], RS, p, c)*fit.smolt[i]
  Sm.hatch[i]<-Hatch.sm(BS[i], fec, release.surv)
  
  fit.adult[i]<-fit.lifestage(P.nat[1], Theta.nat, w, sig, 1-rel.loss)
  ret.nat[i]<-Sm.nat[i]*mar.surv*fit.adult[i]
  ret.hatch[i]<-Sm.hatch[i]*mar.surv.hatch
  }#End if i generations
  return(list(fit.smolt=fit.smolt, P.nat=P.nat, P.hatch=P.hatch, Theta.hatch=Theta.hatch, Theta.nat=Theta.nat, Sp.nat=Sp.nat, ret.nat=ret.nat, BS=BS, Seq=Seq, NOS=NOS, NOB=NOB, HOS=HOS, HOB=HOB, pNOB=pNOB, pHOS=pHOS))
}#End of run.model

panel.plots<-function(){
  par(mfrow=c(2,2), mar=c(3,4,2,1))
  plot(1:100,res$fit.smolt^2, type="l", ylim=c(0,1), ylab="Fitness", xlab="Generation")
  title<-paste("pHOS=",res$pHOS," pNOB=",res$pNOB)
  mtext(title, side=3, line=0.5, at=110, cex=0.8)
  plot(1:100,res$P.nat, type="l", col="blue", ylim=c(0,100), ylab="Phenotype", xlab="Generation")
  points(x=1:100,y=res$P.hatch, col="red", type="l")
  abline(h=res$Theta.nat, col="blue", lty="dashed")
  abline(h=res$Theta.hatch, col="red", lty="dashed")
  legend(x=5,y=40, legend=c("Natural origin", "Hatchery origin"), col=c("red", "blue"), cex=0.5, lty=c("solid", "solid"), bty="n")
  
  y.upper.lim<-res$Seq*1.2
  if(max(res$Sp.nat)>res$Seq*1.2){y.upper.lim<-max(res$Sp.nat)}
     
  plot(1:100, res$Sp.nat, type="l", col="red", ylim=c(0,y.upper.lim), ylab="Spawners (HOS+NOS)", xlab="Generation")
  points(x=1:100, y=res$NOS, col="red", type="l", lty="dashed")
  points(x=1:100, y=res$HOS, col="red", type="l",lty="dotted")
  points(x=1:100, y=res$BS, col="blue", type="l")
  points(x=1:100, y=res$NOB, col="blue", type="l",lty="dashed")
  points(x=1:100, y=res$HOB, col="blue", type="l",lty="dotted")
  legend(x=5,y=250, legend=c("Natural spawners", "NOS", "HOS", "Brood stock", "NOB", "HOB"), col=c(rep("red",3),rep("blue",3)), cex=0.5, lty=rep(c("solid", "dashed","dotted"),2), bty="n")
  
  y.upper.lim<-res$Seq*1.2
  if(max(res$ret.nat)>res$Seq*1.2){y.upper.lim<-max(res$ret.nat)*1.1}
  plot(1:100, res$ret.nat, type="l", col="red", ylim=c(0,y.upper.lim), ylab="Returns from natural spawners", xlab="Generation")
  points(x=1:100, res.nogenetics$ret.nat, type="l", col="dark green")
  #points(x=1:100, res.nogenetics.nodemog$ret.nat, type="l", col="purple")
  #legend(x=5,y=200, legend=c("Without genetic impacts \nassuming equal survival of HOS and NOS","Without genetic impacts & reduced survival of HOS",  "With genetic impacts & reduced survival of HOS"), cex=0.5, col=c("purple", "dark green", "red"), lty=c("solid", "solid"), bty="n")
  legend(x=5,y=200, legend=c("Without genetic impacts",  "With genetic impacts"), cex=0.5, col=c("dark green", "red"), lty=c("solid"), bty="n")
  
}

#pdf("Genetic Impacts Plots 1Dec2016.pdf")

#Model is run here: Wild scenario
res<-Run.model(Theta.hatch=80, RS=0.8, pNOB=0.9, pHOS=0.1)
res.nogenetics<-Run.model(Theta.hatch=100, RS=0.8, pNOB=0.9, pHOS=0.1)
res.nogenetics.nodemog<-Run.model(Theta.hatch=100, RS=1, pNOB=0.9, pHOS=0.1)
#fix the nodemog run- no correct that equil values >> c when pHOS high I don't think???
panel.plots()


#Model is run here: Natural scenario
res<-Run.model(Theta.hatch=80, RS=0.8, pNOB=0.67, pHOS=0.33)
res.nogenetics<-Run.model(Theta.hatch=100, RS=0.8, pNOB=0.67, pHOS=0.33)
res.nogenetics.nodemog<-Run.model(Theta.hatch=100, RS=1, pNOB=0.67, pHOS=0.33)
panel.plots()

#Model is run here: Production scenario
res<-Run.model(Theta.hatch=80, RS=0.8, pNOB=0.5, pHOS=0.5)
res.nogenetics<-Run.model(Theta.hatch=100, RS=0.8, pNOB=0.5, pHOS=0.5)
res.nogenetics.nodemog<-Run.model(Theta.hatch=100, RS=1, pNOB=0.5, pHOS=0.5)
panel.plots()

#Model is run here: Hatchery scenario
res<-Run.model(Theta.hatch=80, RS=0.8, pNOB=0.2, pHOS=0.8)
res.nogenetics<-Run.model(Theta.hatch=100, RS=0.8, pNOB=0.2, pHOS=0.8)
res.nogenetics.nodemog<-Run.model(Theta.hatch=100, RS=1, pNOB=0.2, pHOS=0.8)
panel.plots()

#dev.off()