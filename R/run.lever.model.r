# Function to estimate effect of 3 hatchery levers: hatchery size, marking, and 
# selective harvest, based on HSRG model of genetic impacts of hatcheries on 
# fitness
# Code documented in Withler et al. 2018 CSAS Res. Doc.
# Last updated 23 Feb. 2023

# Source base functions
source("R/FuncDefs.R")
source("R/SeqFormulationBH.R")

plot.UnivariateSA<-FALSE
plot.timeseries<-FALSE

run.lever.model<-function( per.mark, #% of hatchery origin fish that are marked
                           hatchery.size, #Hatchery size as a ppn of Seq
                           sel, #Selective removals/harvest
                           Theta.hatch, 
                           c, #Beverton-Holt capacity parameter (Rmax)
                           percent.hatch, 
                           HR, #Harvest rate
                           h, #heritability
                           w, #selection
                           mar.surv, # Marine survival of natural-origin fish
                           RS, #Relative survival/fecundity of hatchery vs. NOS
                           mar.surv.hatch, #Marine survival of HOS 
                           sex.ratio, # sex ratio of brood stock
                           ppn.RR = NA, # brood stock limit as ppn of returns 
                           BS.ppnRR = FALSE, # is Brood stock determined as a ppn 
                           # of returns to river? (if FALSE, then set based on 
                           # hatchery.size)
                           BS.invRR = FALSE, # is Brood stock determined as an
                           # inverse of ppn of marked returns to river? 
                           p = 175, #Beverton-Holt productivity parameter
                           bs.surv = 0.80,  #survival from capture to spawning for broodstock
                           fec = 4900,  #fecundity
                           release.surv = 0.88 #combined survival of hatchery fish from egg-fry and 
                           ){

n_gen <- 100#11
# # Levers
# per.mark <- 0 #Percentage of hatchery origin fish that are marked
# hatchery.size <- 0.5 #Hatchery size as a ppn of equilibrium capacity (Seq) of
# #  natural spawners
# sel<-0#Selective removals/harvest
# ppn.RR <- 0.1 #allowable take of total returns to river per year for BS
# 
# Biological Parameters
# p <- 175#Beverton-Holt productivity parameter
# c <- 400000#Beverton-Holt capacity parameter (Rmax)
# mar.surv <- 0.02# Marine survival of natural-origin fish
# RS <- 0.8#Relative survival/fecundity of hatchery vs. natural-origin spawners
# # in natural environment in BH model

# Hatchery Parameters
# bs.surv <- 0.80 #survival from capture to spawning for broodstock
# fec <- 4900 #fecundity
# release.surv <- 0.88 #combined survival of hatchery fish from egg-fry and 
# fry-smolt
# mar.surv.hatch <- 0.0024 #Marine survival of hatchery-origin fish
# percent.hatch <- 0 # 1-percent of hatchery fish that return to natural
# # spawning grounds. Look at a range of values
rel.loss <- 0.5 #Ppn of fitness loss at smolt stage (compared with adult stage)
#HR <- 0.4 #Harvest rate

#Genetic parameters
sig <- sqrt(10)
# w <- sqrt(100) #Selection
# h <- sqrt(0.25) #Heritability
Theta.nat <- 100
# Theta.hatch <- 80

#Initialize vectors
P.nat <- NA #Mean phenotypic value in natural environment
P.hatch <- NA #Mean phenotypic value in the hatchery
HOS <- NA
HOR <- NA #Hatchery-origin returns
HOR.unmark <- NA
HOR.mark <- NA
HOR.rem <- NA
PNI <- NA
pHOS <- NA
pHOSEff <- NA
pNOB <- NA
NOS <- NA
HOB <- NA
NOB <- NA
ret.nat <- NA #Returns to the spawning grounds AFTER harvest
ret.nat.preharvest <- NA #Returns to the spawning grounds PRIOR TO harvest
ret.hatch <- NA #Returns to the hatchery AFTER harvest
ret.hatch.preharvest <- NA #Returns to the hatchery PRIOR TO harvest
BS <- NA #Brood stock collected after considering limit on ppn.RR
BS.mark <- NA # Is BS limit exceeded?
BS.limit <- 0
handling.limit.mark <- NA # Are returns < handling size limit?
handling.limit.year <- NA # Years in which handling limts are exceeded
ppn.unmarked.spawners <- NA
Sp.nat <- NA #Total spawners in the natural environment
Sm.nat <- NA #Total smolts from the natural environment
Sm.hatch <- NA #Total smolts from the hatchery
fit.smolt <- NA #fitness mulitiplier at smolt stage
fit.adult <- NA #fitness mulitiplier at adult stage
RperS <- NA #Recruits from natural spawners  prior to harvest /
# (natural spawners=HOS+NOS)
catch <- NA
ext <- 0

#Initial parameters
P.nat[1] <- Theta.nat
P.hatch[1] <- Theta.nat
# Seq assuming BH formulation #2 from Hilborn and Walters (1992), spawner 
# abundances at equilibrium adult recruitment# See "SeqFormulation.r"
Seq <- (c*mar.surv*(1-HR)) * ((p*mar.surv * (1 - HR)) - 1) / 
  (p * mar.surv * (1 - HR)) 
if (mar.surv.const == TRUE){
  Seq <- (c*0.02 * (1 - HR)) * ((p * 0.02 * (1 - HR))-1) / (p * 0.02*(1 - HR)) 
  }
Req <- Seq

# ----------------------------------------------------------------------------
# Population Dynamics
# ----------------------------------------------------------------------------

# NOT USED:
# Brood stock estimation in year prior to the first generation, required to 
# estimate hatchery returns and HOS
# pred.hatch.recruits<-solver.HS(hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$fit
# NR.check<-ratio.func(pred.hatch.recruits, hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$Ret#just to CHECK what solver used for total natural returns
# HOS.check<-ratio.func(pred.hatch.recruits, hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$HOS[100]#just to CHECK what solver used for total natural returns
# NOS.check<-ratio.func(pred.hatch.recruits, hatchery.size, HR, RS, p, c, mar.surv,percent.hatch,  fec, release.surv, mar.surv.hatch, bs.surv)$NOS[100]#just to CHECK what solver used for total natural returns
# BS.set<-pred.hatch.recruits/(fec*release.surv*mar.surv.hatch*(1HR))#The broodstock required to acheive predicted hatchery returns

# Hatchery size as a function of long-term equilibrium recruitment in the 
# absence of the hatchery
BS.set <- hatchery.size * Req
ret.hatch.minus1 <- BS.set * bs.surv * sex.ratio * fec * release.surv * 
  mar.surv.hatch *(1 - HR) # Assuming no fitness in the first year 

# Hatchery-origin returns to natural spawning grounds
HOR[1] <- ret.hatch.minus1 * (1 - percent.hatch)
# Number of HOR that are marked and unmarked
HOR.mark[1] <- HOR[1] * per.mark
HOR.unmark[1] <- HOR[1] * (1 - per.mark) 
#Number of marked HOR that are selectively removed
HOR.rem[1] <- HOR[1] * per.mark * sel 

# BS is equal to BS.set unless BS.set is larger than 0.33 * returns to river
# changed to 0.33 from ppn.RR, as ppn.RR is used for new analyes below, exclusively
if (BS.set > (0.33 * (HOR[1] + Seq - HOR.rem[1]) ) ){
  BS[1] <- 0.33 * (HOR[1] + Seq - HOR.rem[1])
  }
if (BS.set <= (0.33 * (HOR[1] + Seq - HOR.rem[1]) ) ){
  BS[1]<-BS.set
  }

# If BS is set based on ppn of returns to river regardless of Req (BS.Set)
# overwrite the BS[1] above
if(BS.ppnRR) {
  BS[1] <- ppn.RR * (HOR[1] + Seq - HOR.rem[1])
  # if(BS[1] > BS.set) {
  #   BS[1] <- BS.set
  #   BS.limit <- 1
  # }
}

# What is the ppn of spawners that are unmarked?
ppn.unmarked.spawners [1] <- (HOR.unmark[1] + Seq) / 
  (HOR.unmark[1] + HOR.mark[1] * (1 - sel) + Seq)


# If brood can be collected from unmarked fish with a limit of handling effort 
# to sample size of 3 X BS, only if BS is not taken in inv ppns of marking in spawner
if(!BS.invRR) {
  if (ppn.unmarked.spawners[1] * BS[1] * 3 >= BS[1]){
    # pNOB is #wild spawners/wild spawners + unmarked hatchery spawners
    pNOB[1] <- Seq / (Seq + HOR.unmark[1])
    # NOS includes only a portion of BS since some come from unmarked hatchery fish
    NOS[1] <- Seq - BS[1] * (Seq / ((HOR.unmark[1] + Seq)))
    if(NOS[1] < 0) NOS[1]<-0
    # HOS must subtract brood take in proportion that HOR provide marked fish to 
    # returns
    HOS[1] <- HOR[1] - HOR.rem[1] - BS[1] * (HOR.unmark[1] / 
                                               (HOR.unmark[1] + Seq))
    
  }#End of if(ppn.unmarked.spawners[1]*BS[1]*3>BS[1]){
  
  # If brood cannot be collected from unmarked fish with a limit on handling 
  # effort to sample size of  3XBS
  if(ppn.unmarked.spawners[1] * BS[1] * 3 < BS[1]) {
    #Number of unmarked BS taken within handling limits (<BS.set)
    BS.underhandling <- ppn.unmarked.spawners[1] * BS[1] * 3
    # Proportion of BS taken after handling limit exceeded. BS taken in the ratio 
    # of HOR and NOR occurring 
    ppnBS.overhandling <- 1 - BS.underhandling / BS[1]
    # Proportion of BS taken before handling limit exceeded. BS taken in the 
    # ration of unmarked HOR and NOR occurring
    ppnBS.underhandling <- BS.underhandling / BS[1]
    
    BS.NO.underhandling <- BS[1] * 
      (Seq / (HOR.unmark[1] + Seq) ) * ppnBS.underhandling
    BS.NO.overhandling <- BS[1] * 
      (Seq / (HOR[1] - HOR.rem[1] + Seq) ) * ppnBS.overhandling
    BS.HO.underhandling <- BS[1] * 
      (HOR.unmark[1] / (HOR.unmark[1] + Seq) ) * ppnBS.underhandling
    BS.HO.overhandling <- BS[1] * 
      ((HOR[1] - HOR.rem[1]) / (HOR[1] - HOR.rem[1] + Seq) ) * ppnBS.overhandling
    
    NOS [1] <- Seq - BS.NO.underhandling - BS.NO.overhandling
    #if(NOS[1]<0)NOS[1]<-0
    
    # Combing pNOB for ppn fish prior to and after limit exceeded
    pNOB[1] <- Seq / (Seq + HOR.unmark[1] ) * ppnBS.underhandling + 
      Seq / (Seq + HOR[1] - HOR.rem[1]) * ppnBS.overhandling 
    # HOS is HOR minus selectively harvested fish and brood take
    HOS[1] <- HOR[1] - HOR.rem[1] - BS.HO.underhandling - BS.HO.overhandling
    
  }#End of if(ppn.unmarked.spawners[1]*BS[1]*3<BS[1]){
}
if(BS.invRR){
  BS.HO <- BS[1] * ppn.unmarked.spawners[1]
  BS.NO <- BS[1] - BS.HO
  if(BS.HO >= (HOR.mark[1] * (1 - sel))){ 
    BS.HO <-   (HOR.mark[1] * (1 - sel))
  }
  if(BS.NO >= (HOR.unmark[1] + Seq)){ 
    BS.NO <-   (HOR.unmark[1]  + Seq)
    BS.HO <- BS[1] - BS.NO
  }
  
  NOS[1] <- Seq - BS.NO
  pNOB[1] <- Seq / (Seq + HOR[1] - HOR.rem[1]) 
  HOS[1] <- HOR[1] - HOR.rem[1] - BS.HO
}

NOB[1] <- BS[1] * pNOB[1]
HOB[1] <- BS[1] * (1 - pNOB[1])
pHOS[1] <- HOS[1] / (HOS[1] + NOS[1])
pHOSEff[1] <- HOS[1] * RS / (HOS[1] * RS + NOS[1])
PNI[1] <- pNOB [1] / (pNOB[1] + pHOS[1])

Sp.nat[1] <- NOS[1] + HOS[1]

fit.smolt[1] <- fit.lifestage( P.nat[1], Theta.nat, w, sig, rel.loss)
Sm.nat[1] <- BH(HOS[1], NOS[1], RS, p, c) * fit.smolt[1]
Sm.hatch[1] <- Hatch.sm( BS[1] * bs.surv, fec, sex.ratio, release.surv) #*  fit.smolt[1]

fit.adult[1] <- fit.lifestage( P.nat[1], Theta.nat, w, sig, 1 - rel.loss)
# Returns to spawning grounds AFTER harvest
ret.nat[1] <- Sm.nat[1] * mar.surv * fit.adult[1] * (1-HR)
#Returns to hatchery BEFORE harvest
ret.hatch.preharvest[1] <- Sm.hatch[1] * mar.surv.hatch
#Returns to hatchery AFTER harvest
ret.hatch[1] <- Sm.hatch[1] * mar.surv.hatch * (1 - HR) * fit.adult[1]
#Returns to hatchery BEFORE harvest
ret.nat.preharvest[1] <- Sm.nat[1] * mar.surv * fit.adult[1]
# Catch includes harvest + fish selectively removed from river after harvest
catch[1] <- (ret.hatch.preharvest[1] + ret.nat.preharvest[1]) * HR + 
  ret.hatch[1] * (1 - percent.hatch) * (per.mark) * sel
RperS[1] <- (ret.nat.preharvest[1] + catch[1]) / (HOS[1] + NOS[1])

for (i in 2:n_gen){#for i generations)
  # If population NOT extirpated (natural population extirpated because fitness 
  # impacts and 100% hatchery fish removed prior to spawning)
  if (ext == 0) {
    # Hatchery-origin returns to natural spawning grounds
    HOR[i] <- ret.hatch[i-1] * (1 - percent.hatch) 
    HOR.unmark[i] <- HOR[i] * (1 - per.mark) # Number of HOR that are unmarked
    HOR.mark[i] <- HOR[i] * (per.mark) # Number of HOR that are unmarked
    #Number of marked HOR that are selectively removed
    HOR.rem[i] <- HOR.mark[i] * sel 

    if(ret.nat[i-1]>0){
      # We limited the total broodstock to less than a third of the total 
      # returns to the river in any given year to avoid conservation concerns 
      # from removing too many fish for brood.
      # Changed from ppn.RR to 0.33 as ppn.RR used exclusively for new analyses below
      if (BS.set > (0.33 * (ret.nat[i-1] + HOR[i] - HOR.rem[i]) ) ) {
        BS[i] <- 0.33 * (ret.nat[i-1] + HOR[i] - HOR.rem[i])
        }
      if (BS.set <= (0.33 * (ret.nat[i-1] + HOR[i] - HOR.rem[i]) ) ) {
        BS[i] <- BS.set
        }
      
      # However, if BS is set based on ppn of returns to river regardless of Req
      # (BS.set) overwrite the BS[i] above
      if(BS.ppnRR) {
        BS[i] <- ppn.RR * (ret.nat[i-1] + HOR[i] - HOR.rem[i])
        # if(BS[i] > BS.set) {
        #   BS[i] <- BS.set
        #   BS.limit[i] <- 1
        # }
        
      }
      
      # Check if BS (and handing limit?) are actually available on the spawning 
      # grounds. Sometimes BS.set is a allowed by ppn.RR rule, but can't be 
      # achieved because there are no/too few unmarked spawners
      
      # BS.mark indicates if the target bs identified in 'hatchery.size' 
      # management lever is achieved in the last year/equilibrium state. In 
      # some cases (BS.mark =1), the limitation of BS<1/3 of returns to river is 
      # used exclusively over 100 generations
      if (i == 100) {
        if (BS[i] == BS.set) {BS.mark <- 0}
        if (BS[i] != BS.set) {BS.mark <- 1}
        }
      
      # The ppn of unmarked spawner
      ppn.unmarked.spawners[i] <- (HOR.unmark[i] + ret.nat[i-1]) / 
        (HOR.unmark[i] + HOR.mark[i] * (1 - sel) + ret.nat[i-1])
      

    
      # Can brood  be collected from unmarked fish with a limit of handling 
      # effort to sample size of 3XBS? Only use when BS is ppn of Seq or RR
      # not when BS if ppn of invRR
      if(!BS.invRR){
        
        handling.limit <- BS[i] * 3
        
        # If handling limit is available on the spawning grounds, sample for BS
        if ( handling.limit <  
             (HOR.unmark[i] + HOR.mark[i] * (1 - sel) + ret.nat[i-1]) ){
          if(i == 100){handling.limit.mark <- 0}
          
          # If handling limits are not exceeded to find unmarked brood for the BS 
          # target given ppn of unmarked spawners, then collect unmarked BS
          if ( (ppn.unmarked.spawners[i] * handling.limit) >= BS[i]) {
            # NOS is natural returns minus the portion of BS that is natural-origin
            NOS[i] <- ret.nat[i-1] - 
              BS[i] * (ret.nat[i-1] / (HOR.unmark[i] + ret.nat[i-1]))
            # Can the brood be collected from unmarked fish within a limit of the 
            # total number unmarked spawners on the spawning grounds?
            # If not, NOS is negative
            if (NOS[i] < 0) {browser(); cat(NOS[i])} #NOS[i]<-0
            
            # HOS is HORs minus the portion of BS that is hatchery-origin
            HOS[i] <- HOR.unmark[i] + HOR.mark[i] * (1 - sel) - 
              BS[i] * (HOR.unmark[i] / (HOR.unmark[i] + ret.nat[i-1]) )
            
            # Brood are taken from unmarked fish only, on the spawning grounds
            # pNOB is natrual return/natural return + unmarked hatchery spawners
            pNOB[i]<-(ret.nat[i-1])/(ret.nat[i-1]+HOR.unmark[i])
          }#End of if(ppn.unmarked.spawners[i]* handling.limit > BS[i]){
          
          #If brood cannot be collected from unmarked fish with a limit on 
          # handling effort to sample size of  3XBS
          if(ppn.unmarked.spawners[i] * BS[i] * 3 < BS[i]) {
            #Number of unmarked BS taken within handling limits (<BS[i])
            BS.underhandling <- ppn.unmarked.spawners[i] * BS[i] * 3
            # Proportion of BS taken after handling limit exceeded. BS taken in 
            # the ratio of HOR and NOR occurring
            ppnBS.overhandling <- 1 - BS.underhandling / BS[i]
            # Proportion of BS taken before handling limit exceeded. BS taken in 
            # the ration of unmarked HOR and NOR occurring
            ppnBS.underhandling <- BS.underhandling/BS[i]
            
            BS.NO.underhandling <- BS[i] * 
              (ret.nat[i-1] / ( (HOR.unmark[i] + ret.nat[i-1]) ) ) * 
              ppnBS.underhandling
            
            BS.NO.overhandling <- BS[i] * 
              (ret.nat[i-1] / (HOR[i] - HOR.rem[i] + ret.nat[i-1]) ) * 
              ppnBS.overhandling
            
            BS.HO.underhandling <- BS[i] * 
              (HOR.unmark[i] / ( (HOR.unmark[i] + ret.nat[i-1]) ) ) * 
              ppnBS.underhandling
            
            BS.HO.overhandling <- BS[i] * 
              ((HOR[i] - HOR.rem[i]) / (HOR[i] - HOR.rem[i] + ret.nat[i-1])) * 
              ppnBS.overhandling
            
            NOS[i] <- ret.nat[i-1] - BS.NO.underhandling -  BS.NO.overhandling
            #if(NOS[i]<0)NOS[i]<-0
            
            # Combing pNOB for unmarked component and marked component 
            pNOB[i] <- (ret.nat[i-1]) / 
              (ret.nat[i-1] + HOR.unmark[i]) * ppnBS.underhandling + 
              (ret.nat[i-1]) / 
              (ret.nat[i-1] + HOR[i] - HOR.rem[i]) * ppnBS.overhandling
            HOS[i] <- HOR[i] - HOR.rem[i] - BS.HO.underhandling - 
              BS.HO.overhandling
            #subtract  brood take  
          }#End of if(ppn.unmarked.spawners[i]*BS[i]*3<BS[i]){
          
        } # End of if (handling.limit > 
        
        # If handling limit is NOT available on the spawning grounds, take bs in 
        # proportions on the spawning grounds
        if ( handling.limit >=  
             (HOR.unmark[i] + HOR.mark[i] * (1 - sel) + ret.nat[i-1]) ){
          handling.limit.year[i] <- 1
          if(i == 100) {handling.limit.mark <- 1}
          BS.NO <- BS[i] * 
            (ret.nat[i-1] / (HOR[i] - HOR.rem[i] + ret.nat[i-1]) ) 
          BS.HO <- BS[i] * 
            ((HOR[i] - HOR.rem[i]) / (HOR[i] - HOR.rem[i] + ret.nat[i-1])) 
          
          NOS[i] <- ret.nat[i-1] - BS.NO
          #if(NOS[i]<0)NOS[i]<-0
          
          pNOB[i] <- (ret.nat[i-1]) / (ret.nat[i-1]+HOR[i]-HOR.rem[i]) 
          HOS[i] <- HOR[i] - HOR.rem[i] - BS.HO
          #subtract  brood take  
          
        } # End of if ( handling.limit >= 
      }
      
    # If broodstrock is invRR, then Brood stock marked proportions = ppn of 
      # unmarked spawners, and total BS is a proportoin of Seq
    if(BS.invRR){
      BS.HO <- BS[i] * ppn.unmarked.spawners[i]
      BS.NO <- BS[i] - BS.HO
      # Make sure there are enough marked/unmarked fish to take BS.HO and BS.NO
      if(BS.HO >= (HOR.mark[i] * (1 - sel))){ 
        BS.HO <-   (HOR.mark[i] * (1 - sel))
      }
      if(BS.NO >= (HOR.unmark[i] + ret.nat[i-1])){ 
        BS.NO <-   (HOR.unmark[i]  + ret.nat[i-1])
        BS.HO <- BS[i] - BS.NO
      }
      
      NOS[i] <- ret.nat[i-1] - BS.NO
      pNOB[i] <- (ret.nat[i-1]) / (ret.nat[i-1]+HOR[i]-HOR.rem[i]) 
      HOS[i] <- HOR[i] - HOR.rem[i] - BS.HO
    }
  }#End of if(ret.nat[i-1]>0)
  
  if(ret.nat[i-1]<=0){
    ppn.unmarked.spawners[i]<-0; 
    NOS[i]<-0; 
    HOS[i]<-HOR[i]-HOR.rem[i]; 
    pNOB[i]<-0; 
    ext<-1
    }#And BS[i]=BS[i] (not 0.33 of ret.nat+HOR)
  
  NOB[i]<-BS[i]*pNOB[i]
  HOB[i]<-BS[i]*(1-pNOB[i])
  pHOS[i]<-HOS[i]/(HOS[i]+NOS[i])
  pHOSEff[i]<-HOS[i]*RS/(HOS[i]*RS+NOS[i])
  PNI[i]<-pNOB[i]/(pNOB[i]+pHOSEff[i])#pNOB[i]/(pNOB[i]+pHOSEff[i])
  Sp.nat[i]<-NOS[i]+HOS[i]
  
  # P.nat[i]<-Pnat(pHOS[i-1], P.nat[i-1], w, sig, Theta.nat, h, P.hatch[i-1])
  P.nat[i]<-Pnat(pHOSEff[i-1], P.nat[i-1], w, sig, Theta.nat, h, P.hatch[i-1])
  P.hatch[i]<-Phatch(pNOB[i-1], P.hatch[i-1], w, sig, Theta.hatch, h, P.nat[i-1])
  fit.smolt[i]<-fit.lifestage(P.nat[i], Theta.nat, w, sig, rel.loss)
  Sm.nat[i]<-BH(HOS[i], NOS[i], RS, p, c)*fit.smolt[i]
  Sm.hatch[i]<-Hatch.sm(BS[i]*bs.surv, fec, sex.ratio, release.surv)#*fit.smolt[i]
  
  fit.adult[i]<-fit.lifestage(P.nat[i], Theta.nat, w, sig, 1-rel.loss)
  ret.nat[i]<-Sm.nat[i]*mar.surv*fit.adult[i]*(1-HR)
  ret.hatch[i]<-Sm.hatch[i]*mar.surv.hatch*(1-HR)*fit.adult[i]
  ret.hatch.preharvest[i]<-Sm.hatch[i]*mar.surv.hatch
  ret.nat.preharvest[i]<-Sm.nat[i]*mar.surv*fit.adult[i]
  catch[i]<-(ret.hatch.preharvest[i]+ret.nat.preharvest[i])*HR + ret.hatch[i]*(1-percent.hatch)*(per.mark)*sel#Includes harvest + fish selectively removed from river after harvest
  RperS[i]<-(ret.nat.preharvest[i] + catch[i])/(HOS[i]+NOS[i])
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

  return( list (fit.smolt = fit.smolt, P.nat = P.nat, P.hatch = P.hatch, 
                Theta.hatch = Theta.hatch, Theta.nat = Theta.nat, 
                Sp.nat = Sp.nat, ret.nat = ret.nat, 
                ret.nat.preharvest = ret.nat.preharvest, 
                ret.hatch.preharvest = ret.hatch.preharvest, catch = catch, 
                BS = BS,#rep (BS.set, 100), 
                Seq = Seq, NOS = NOS, NOB = NOB, 
                HOS = HOS, HOB = HOB, pNOB = pNOB, pHOS = pHOS, PNI = PNI, 
                per.mark = per.mark, hatchery.size = hatchery.size, 
                RperS = RperS, sel = sel, mar.surv = mar.surv, c = c, 
                pHOSeff = pHOSEff, fit.adult = fit.adult, BS.mark = BS.mark,
                handling.limit.mark = handling.limit.mark, 
                handling.limit.year = handling.limit.year,
                Sm.nat = Sm.nat, BS.limit = BS.limit) )

}#End of run.lever.model

panel.plots.lever<-function(res, res.nogenetics){

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
#pdf("TimeseriesLevers19Jan2017.pdf")
# pdf("TimeseriesLevers.pdf")
dir <- here::here("Results")
scenario <- "pHOSeff"#"pHOSeff"# "pHOS"
png(paste(dir, "/TimeseriesLevers", scenario, ".png", sep=""), width=6, 
    height=6, units="in", res=1000)
# pdf(paste(dir, "/TimeseriesLevers", scenario, ".pdf", sep=""))
res <- run.lever.model (per.mark = 1, hatchery.size = 0.5, sel = 0, 
                        Theta.hatch = 80, c = 400000, percent.hatch = 0, 
                        HR = 0.4, h = sqrt(0.25), w = sqrt(100), 
                        mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, 
                        sex.ratio = 0.5,  ppn.RR = 0.1, BS.ppnRR = FALSE, BS.invRR = TRUE)
res.nogenetics <- run.lever.model (per.mark = 1, hatchery.size = 0.5, sel = 0, 
                                   Theta.hatch = 100, c = 400000, 
                                   percent.hatch = 0, HR = 0.4, h = sqrt(0.25), 
                                   w = sqrt(100), mar.surv = 0.02, RS = 0.8, 
                                   mar.surv.hatch = 0.0024, sex.ratio = 0.5,  
                                   ppn.RR = 0.1, BS.ppnRR = FALSE, BS.invRR = TRUE)
panel.plots.lever(res=res, res.nogenetics=res.nogenetics)

dev.off()
}

#per.mark=1;hatchery.size=0.1; sel=1; Theta.hatch=80
if(plot.UnivariateSA==TRUE){

PNIm<-NA; pHOSm<-NA; pNOBm<-NA; RperSm<-NA; Ret.natm<-NA; PNIs<-NA; pHOSs<-NA; pNOBs<-NA; RperSs<-NA; Ret.nats<-NA; PNIh<-NA; pHOSh<-NA; pNOBh<-NA; BSh<-NA; RperSh<-NA; Ret.nath<-NA;  
for(i in 1:100){
  m<-run.lever.model(per.mark=i*0.01, hatchery.size=0.3, sel=0, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024, sex.ratio=0.5)
  m.BS<-m$BS
  m.hatchery.size<-m$hatchery.size
  m.sel<-m$sel
  PNIm[i]<-m$PNI[100]
  pHOSm[i]<-m$pHOSeff[100]
  pNOBm[i]<-m$pNOB[100]
  RperSm[i]<-m$RperS[100]
  Ret.natm[i]<-m$ret.nat.preharvest[100]
  s<-run.lever.model(per.mark=0.1, hatchery.size=0.3, sel=i*0.01, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024, sex.ratio=0.5)
  s.BS<-s$BS
  s.hatchery.size<-s$hatchery.size
  s.per.mark<-s$per.mark
  PNIs[i]<-s$PNI[100]
  pHOSs[i]<-s$pHOS[100]
  pNOBs[i]<-s$pNOB[100]
  RperSs[i]<-s$RperS[100]  
  Ret.nats[i]<-s$ret.nat.preharvest[100]  
  h<-run.lever.model(per.mark=0.1, hatchery.size=i*0.003, sel=0, Theta.hatch=80, c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), mar.surv=0.02, RS=0.8, mar.surv.hatch=0.0024, sex.ratio=0.5)
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