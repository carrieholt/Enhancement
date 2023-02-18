# Functions for model to evalaute genetic impact of hatcheries
# Date last revised: 14 Nov. 2016
# Created by: Carrie Holt
# Equations from HSRG (2009) Appendix C: AHA Model

BH<-function(HOS, NOS, RS, p, c){#Beverton-Holt spawner-recruitment function. Eqns. 1 and 2 HSRG(2009)
  #INPUT: HOS=hatchery-origin spanwers; NOS=natural-origin spawners; RS= relative surival of hatchery fish; p=density-independent survival (slope at origin of R/S); c=capacity (Rmax)
  #OUTPUT: Sm.nat=smolts
  #See Hilborn and Walters (1992) for a description of this BH formulation
  Sm.nat<-p*(NOS+RS*HOS)/(1+(p*(NOS+RS*HOS))/c)
  return(Sm.nat)
}

Fit.total<-function(P.nat, Theta.nat, w, sig){#Fitness over entire generation Eqn. 34 of HSRG(2009)
  #INPUT: P.nat=mean phenotypic value of the natural population in current generation; Theta-nat=Phenotypic optimum for natural population; w=SD of distribution of fitness along gradient in phenotypic values; sig=SD of of phenotype
  #OUTPUT: Fitness= population fitness over an entire generation
  Fitness<-exp(-0.5*((P.nat-Theta.nat)^2/(w^2+sig^2)))
  return(Fitness)
}

fit.lifestage<-function(P.nat, Theta.nat, w, sig, rel.loss){#Fitness by life-stage Eqn. 3 of HSRG (2009)
  #INPUT: Fit.total= fitness over the entire generation; rel.loss=ppn of total fitness attributable to a given life stage
  #OUTPUT: fit= proportion of total fitness attributable to a given life-stage
  fit<-Fit.total(P.nat, Theta.nat, w, sig)^rel.loss
  return(fit)
  
}

Pnat<-function(pHOS, P.nat.minus1, w, sig, Theta.nat, h, P.hatch.minus1){#Mean phenotypic value Eqn. 35 of HSRG (2009)
  #INPUT: pHOS=proportion of hachery-origin spawners, Pnat.mnius1=mean phenotypic value of natural spawners in previous generation, w=SD of distribution of fitness along gradient in phenotypic values; Theta.nat=Phenotypic optimum for natural population; 
  #     sig=SD of of phenotype; h=sqrt of heritability of trait; P.hatch.minus1=mean phenotypic value of hatchery fish in previous generation
  #OUTPUT: P.nat=mean phenotypic value of the natural population in current generation
  P.nat<-(1-pHOS)*(P.nat.minus1 +(((P.nat.minus1*w^2 + Theta.nat*sig^2)/(w^2 + sig^2)) - P.nat.minus1) *h^2) + pHOS*(P.hatch.minus1 + (((P.hatch.minus1*w^2 + Theta.nat*sig^2)/(w^2+sig^2))-P.hatch.minus1)*h^2)
  return(P.nat) 
}

Phatch<-function(pNOB, P.hatch.minus1, w, sig, Theta.hatch, h, P.nat.minus1){#Mean phenotypic value Eqn. 36 of HSRG (2009)
  #INPUT: pNOB=proportion of natural-origin broodstock, Phatch.mnius1=mean phenotypic value of hatchery fish in previous generation, w=SD of distribution of fitness along gradient in phenotypic values; Theta.hatch=Phenotypic optimum for hatchery population; 
  #     sig=SD of of phenotype; h=sqrt of heritability of trait; P.nat.minus1=mean phenotypic value of natural spawnersin previous generation
  #OUTPUT: P.hatch=mean phenotypic value of the hatchery fish in current generation
  P.hatch<-(1-pNOB)*(P.hatch.minus1 +(((P.hatch.minus1*w^2 + Theta.hatch*sig^2)/(w^2 + sig^2)) - P.hatch.minus1) *h^2) + pNOB*(P.nat.minus1 + (((P.nat.minus1*w^2 + Theta.hatch*sig^2)/(w^2+sig^2))-P.nat.minus1)*h^2)
  return(P.hatch) 
}

Hatch.sm<-function(BS, fec, sex.ratio, release.surv){#Number of hathchery smolts
  #INPUT: BS=broodstock, fec=fecundity, release.surv=survival from eggs to smolt 
  #OUTPUT: Sm.hatch= number of hatchery smolts
  Sm.hatch<-BS*fec*release.surv*sex.ratio
  return(Sm.hatch)
}



BS.func<-function(BS, pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat, percent.hatch){#Function to optimize total BS to acheive pHOS given pNOB
  #INPUT: BS=brood stock, pNOB, pHOS, BS=brood stock, fec=fecundity of hatchery spawners, release.surv=survival of hatchery eggs to smolt, mar.surv.hatch=survival of hatchery smolts to return, ret.nat=number of natural-origin returns, percent.hatch=% of hatchery-origin fish recovered or die at point of release
  #OUTPUT: neg log likelihood of deviations between estimated pHOS and specified pHOS
  # Used in functions: fn.BS and solver.BS to optimize broodstock to achieve pHOS, given other input parameters
  NOB<-BS*pNOB 
  if(NOB >0.5*ret.nat){NOB<-0.5*ret.nat}#pNOB not achieved
  HOB<-BS-NOB
  ret.hatch<-BS*fec*release.surv*mar.surv.hatch #number of hatchery returns, assuming equilibrium between previous and current generation
  HOS<-ret.hatch*(1-percent.hatch)# hatchery origin spawners if total number returns to the hatchery x ppn that stray to natural spawning grounds
  Sp.nat<-ret.nat-NOB+HOS#Total spawners in the natural environment
  pHOS.est<-HOS/Sp.nat#Estiamted pHOS given HOS and Sp.nat calculated
  epsilon<-pHOS.est-pHOS#Deviations between estimated and specified pHOS
  nloglike=sum(dnorm(epsilon,0,sd=1, log=T)) 
  return(list(BS=BS, epsilon=epsilon, nloglike=nloglike))  
}

fn.BS <- function(BS, pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat, percent.hatch) -1.0*BS.func(BS, pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat, percent.hatch)$nloglike  #gives the minLL

solver.BS <- function(pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat, percent.hatch) {
    fit=optimize(f=fn.BS,interval=c(0,ret.nat*10), pNOB=pNOB, pHOS=pHOS, fec=fec, release.surv=release.surv, mar.surv.hatch=mar.surv.hatch, ret.nat=ret.nat, percent.hatch=percent.hatch)	
    out<-fit$minimum
    if(out*pNOB>0.5*ret.nat){print("NOB too high for natural returns")}
    return(list(fit=out))
}





#Trials for solver.BS functions, compare to xls "Recruits Accounting for BH model 16Nov2016.xlsx"
#pNOB<-0.67
#pHOS<-0.33
#fec<-100
#release.surv<-0.1
#mar.surv.hatch<-0.9
#ret.nat<-300
#percent.hatch<-0.8
#solver.BS(pNOB, pHOS, fec, release.surv, mar.surv.hatch, ret.nat, percent.hatch)
  
