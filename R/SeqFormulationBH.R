#Code to  calculate Seq when it combines an egg-smolt BH and dens-indep marine survival and harvest rates
#AND code to estimate EStimated hatchery recruitment abundances that would meet hatchery size specified (as a ratio  hatchery recruitment/total recruitment)
#Also estimates long-term average recruitment from spawning grounds (Ret[100]) given a specified input of hatchery fish

#Formulation of BH from Hilborn 1986, as used in AHA model

#Date last revised 19 Dec. 2016

#rm(list = ls())
source("R/FuncDefs.r")

#use this commented out code to run the function
#HR<-0.4
#RS<-1
#ratio<-0.3
#p<-175#Beverton-Holt productivity parameter
#c<-400000#Beverton-Holt capacity parameter (Rmax)
#mar.surv<-0.02# Marine survival of natural-origin fish
#per.hatch<-0#All hatchery fish return to natural spawning grounds
#HOR<-0#1674.189#Hatchery-origin returns #Change this to match Ratio to hatchery size
#fec<-4900#fecundity
#release.surv<-0.88# combined survival of hatchery fish from egg-fry and fry-smolt
#mar.surv.hatch<-0.0024 #Marine survival of hatchery-origin fish
#bs.surv<-0.80#survival from capture to spawning for broodstock


ratio.func<-function(HOR, ratio, HR, RS, p, c, mar.surv,per.hatch, fec, release.surv, mar.surv.hatch, bs.surv){

HOS<-NA; NOS<-NA; Ret<-NA; Adults<-NA
#HR<-0.4      #use this commented out code to work through the steps within this function
#RS<-1
#ratio<-0.5#~hatchery size as a  ratio  of predicted hatchery recruitment/predicted total recruitment
#p<-175#Beverton-Holt productivity parameter
#c<-400000#Beverton-Holt capacity parameter (Rmax)
#mar.surv<-0.01#0.02# Marine survival of natural-origin fish
#per.hatch<-0#All hatchery fish return to natural spawning grounds
#HOR<-0#112#Hatchery-origin returns #Change this to match Ratio to hatchery size
#fec<-4900#fecundity
#release.surv<-0.88# combined survival of hatchery fish from egg-fry and fry-smolt
#mar.surv.hatch<-0.0024 #Marine survival of hatchery-origin fish
#bs.surv<-0.80#survival from capture to spawning for broodstock

#Can add HR here to estimate what is sustainable under a given HR
Seq<-((c*mar.surv*(1-HR))*((p*mar.surv*(1-HR))-1)/(p*mar.surv*(1-HR)))#*mar.surv
NOS[1]<-Seq
HOS[1]<-HOR*(1-per.hatch)
Ret[1]<-mar.surv*BH(HOS[1], NOS[1], RS, p, c)
Adults[1]<-mar.surv*BH(HOS[1], NOS[1], RS, p, c)

for (i in 2:100){
  NOS[i]<-Ret[i-1]-HOR/(fec*release.surv*mar.surv.hatch*(1-HR)*bs.surv)# NOS= Returns minus BS, assuming BS comes from NOS only
  HOS[i]<-HOR*(1-per.hatch)
  Adults[i]<-mar.surv*BH(HOS[i], NOS[i], RS, p, c)#Returns before harvest
  Ret[i]<-mar.surv*BH(HOS[i], NOS[i], RS, p, c)*(1-HR)#Returns after harvest
}

ratio.est<-HOR/(Ret[100]+HOR) #estimated ratio  of predicted hatchery recruitment/predicted total recruitment
epsilon<-ratio.est-ratio#Deviations between estimated and specified ratio (or hatchery size)
nloglike=sum(dnorm(epsilon,0,sd=1, log=T)) 
return(list(HOR=HOR, epsilon=epsilon, Ret=Ret[100], nloglike=nloglike, HOS=HOS, NOS=NOS))  
}

fn.ratio <- function(HOR, ratio, HR, RS, p, c, mar.surv,per.hatch, fec, release.surv, mar.surv.hatch, bs.surv) -1.0*ratio.func(HOR, ratio, HR, RS, p, c, mar.surv,per.hatch, fec, release.surv, mar.surv.hatch, bs.surv)$nloglike  #gives the minLL

solver.HS <- function(ratio, HR, RS, p, c, mar.surv,per.hatch, fec, release.surv, mar.surv.hatch, bs.surv) {
  fit=optimize(f=fn.ratio,interval=c(0,c*mar.surv), ratio=ratio, HR=HR, RS=RS, p=p,c=c, mar.surv=mar.surv, per.hatch=per.hatch, fec=fec, release.surv=release.surv, mar.surv.hatch=mar.surv.hatch, bs.surv=bs.surv)  
  out<-fit$minimum
  return(list(fit=out))# EStimated hatchery recruitment abundances that would meet hatchery size specified (as a ratio  hatchery recruitment/total recruitment)
}

#Check that solver provides HOR that actually does achieve the target ratio (hatchery size)
#HOR<-solver.HS(ratio, HR, RS, p, c, mar.surv,per.hatch, fec, release.surv, mar.surv.hatch, bs.surv)$fit
#NR<-ratio.func(HOR, ratio, HR, RS, p, c, mar.surv,per.hatch, fec, release.surv, mar.surv.hatch, bs.surv)$Ret #Total returns from spawners (HOS and NOS) NOT accounting for genetic effects
#HOR/(HOR+NR)
