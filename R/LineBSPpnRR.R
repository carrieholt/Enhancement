#***************************************************************************
# Code to run line plots over various targets on brood stock as specified 
# ppns of returns to river that determine hatchery size, assuming 0% marking
# Date last revised 7 July 2023
#***************************************************************************

library(ggplot2)
library(RColorBrewer)

# Source population dynamics models
source(here::here("R", "run.lever.model.r"))
source(here::here("R", "FuncDefs.R"))

mar.surv.const<-TRUE#Is the marine suvival assumed to be 0.02 when estimating 
# Seq, or derived from inputed marine survival?

BS.ppnRR <- TRUE
BS.invRR <- TRUE

lineSA <- function(scenario, 
                   h, 
                   w, 
                   mar.surv, 
                   mar.surv.hatch, 
                   RS, 
                   Theta.hatch, 
                   p, 
                   c, 
                   bs.surv, 
                   fec, 
                   release.surv,
                   BS.ppnRR = TRUE,
                   BS.invRR = FALSE
                   ){
  # First, set up empty vectors
  PNI <- NA
  pHOS <- NA
  pNOB <- NA   
  RperS <- NA
  Ret.nat <- NA
  Ret.hat <- NA 
  Catch <- NA
  fit <- NA
  BSmark <- NA
  HLmark <- NA
  hatchery.size <- NA
  BS <- NA
  
  sel <- 0#0.5
  per.mark <- 1#0#1
  percent.hatch<-0# 1-percent of hatchery fish that return to natural spawning 
  HR<-0.4#Harvest rate
  mar.surv.const<-TRUE#Is the marine suvival assumed to be 0.02 when estimating 
  # Seq, or derived from inputed marine survival?
  sex.ratio <- 0.5
  
  for(i in 1:100){
    # if(BS.invRR == TRUE) {c <- i * 10000}
    # if(BS.invRR == TRUE) {per.mark <- i * 0.01}
    mark0 <- run.lever.model (per.mark = per.mark, 
                              hatchery.size = 0.15, # Not used if BS.ppnRR=TRUE
                              sel = sel, 
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
                              ppn.RR = i*0.005, #Not used if BS.ppnRR = FALSE
                              BS.ppnRR = BS.ppnRR,
                              BS.invRR = BS.invRR,
                              p = p, #Beverton-Holt productivity parameter
                              bs.surv = bs.surv,  #survival from capture to spawning for broodstock
                              fec = fec,  #fecundity
                              release.surv = release.surv #combined survival of hatchery fish from egg-fry and 
    )
    
    hatchery.size[i] <- mark0$hatchery.size
    BS[i] <- mark0$BS[100]
    
    PNI[i] <- mark0$PNI[100]
    pHOS[i] <- mark0$pHOSeff[100]
    pNOB[i] <- mark0$pNOB[100]
    RperS[i] <- mark0$RperS[100]
    Ret.nat[i] <- mark0$ret.nat.preharvest[100]
    Ret.hat[i] <- mark0$ret.hatch.preharvest[100]
    Catch[i] <- mark0$catch[100]
    fit[i] <- mark0$fit.adult[100]
    BSmark[i] <- mark0$BS.mark
    HLmark[i] <- mark0$handling.limit.mark
  }
  return(list( scenario = scenario,
               hatchery.size = hatchery.size, 
               BS = BS, 
               PNI = PNI, 
               pHOS = pHOS, 
               pNOB = pNOB, 
               RperS = RperS, 
               Ret.nat = Ret.nat, 
               Ret.hat = Ret.hat, 
               Catch = Catch, 
               fit = fit, 
               BSmark = BSmark, 
               HLmark = HLmark))
}


run.scenarios <- function(BS.ppnRR = TRUE, BS.invRR = FALSE){
  
  baseCase <- lineSA (scenario = "Base case",
                      h = sqrt(0.25), 
                      w = sqrt(100), 
                      mar.surv = 0.02, 
                      mar.surv.hatch = 0.0024, 
                      RS = 0.8, 
                      Theta.hatch = 80, 
                      p = 175, 
                      c = 400000,
                      bs.surv = 0.8, 
                      fec = 4900, 
                      release.surv = 0.88, 
                      BS.ppnRR = BS.ppnRR,
                      BS.invRR = BS.invRR)
  
  
  
  # Selectoin pressure w
  # Withler et al. 2018 used sqrt(100) with sqrt(1000) and sqrt(40) as sensitivity analyses
  # Bradford et al. 2023 used sqrt(50) with sqrt(30) as sensitivity 
  
  w30 <- lineSA (scenario = "Strong selection",
                 h = sqrt(0.25), 
                 w = sqrt(30), 
                 mar.surv = 0.02, 
                 mar.surv.hatch = 0.0024, 
                 RS = 0.8, 
                 Theta.hatch = 80, 
                 p = 175, 
                 c = 400000,
                 bs.surv = 0.8, 
                 fec = 4900, 
                 release.surv = 0.88, 
                 BS.ppnRR = BS.ppnRR,
                 BS.invRR = BS.invRR)
  
  
  w1000 <- lineSA (scenario = "Weak selection",
                   h = sqrt(0.25), 
                   w = sqrt(1000), 
                   mar.surv = 0.02, 
                   mar.surv.hatch = 0.0024, 
                   RS = 0.8, 
                   Theta.hatch = 80, 
                   p = 175, 
                   c = 400000,
                   bs.surv = 0.8, 
                   fec = 4900, 
                   release.surv = 0.88, 
                   BS.ppnRR = BS.ppnRR,
                   BS.invRR = BS.invRR)
  # Heritability
  # Withler et al. 2018 used sqrt(0.25) with sqrt(0.05) and sqrt(0.5) as sensitivity analyses
  # Bradford et al. 2023 used sqrt(0.25)
  h0.05 <- lineSA (scenario = "Low heritability",
                   h = sqrt(0.05), 
                   w = sqrt(100), 
                   mar.surv = 0.02, 
                   mar.surv.hatch = 0.0024, 
                   RS = 0.8, 
                   Theta.hatch = 80, 
                   p = 175, 
                   c = 400000,
                   bs.surv = 0.8, 
                   fec = 4900, 
                   release.surv = 0.88, 
                   BS.ppnRR = BS.ppnRR,
                   BS.invRR = BS.invRR)
  
  h0.5 <- lineSA (scenario = "High heritability",
                  h = sqrt(0.5), 
                  w = sqrt(100), 
                  mar.surv = 0.02, 
                  mar.surv.hatch = 0.0024, 
                  RS = 0.8, 
                  Theta.hatch = 80, 
                  p = 175, 
                  c = 400000,
                  bs.surv = 0.8, 
                  fec = 4900, 
                  release.surv = 0.88, 
                  BS.ppnRR = BS.ppnRR,
                  BS.invRR = BS.invRR)
  
  
  # Optimal trait value for hatchery (natural = 100)
  # Withler et al. used 80
  # Bradford et al. used 60
  thetaHatch60 <- lineSA (scenario = "Large difference in optimal trait",
                          h = sqrt(0.25), 
                          w = sqrt(100), 
                          mar.surv = 0.02, 
                          mar.surv.hatch = 0.0024, 
                          RS = 0.8, 
                          Theta.hatch = 60, 
                          p = 175, 
                          c = 400000,
                          bs.surv = 0.8, 
                          fec = 4900, 
                          release.surv = 0.88, 
                          BS.ppnRR = BS.ppnRR,
                          BS.invRR = BS.invRR)
  
  # Relative reproductive success of hatcher-origin fish
  # Withler et al. use 0.8 with 0.2 in sensitivity analysis
  # Bradford et al. 2023 used 0.8
  RS0.2 <- lineSA (scenario = "Low RRS",
                   h = sqrt(0.25), 
                   w = sqrt(100), 
                   mar.surv = 0.02, 
                   mar.surv.hatch = 0.0024, 
                   RS = 0.2, 
                   Theta.hatch = 80, 
                   p = 175, 
                   c = 400000,
                   bs.surv = 0.8, 
                   fec = 4900, 
                   release.surv = 0.88, 
                   BS.ppnRR = BS.ppnRR,
                   BS.invRR = BS.invRR)
  
  #Beverton-Holt productivity parameter
  # Withler et al. 2018 used 175
  # Bradford et al. 2023 use 138 smolts/spawner
  
  p138 <- lineSA (scenario = "Low productivity",
                  h = sqrt(0.25), 
                  w = sqrt(100), 
                  mar.surv = 0.02, 
                  mar.surv.hatch = 0.0024, 
                  RS = 0.8, 
                  Theta.hatch = 80, 
                  p = 138, 
                  c = 400000,
                  bs.surv = 0.8, 
                  fec = 4900, 
                  release.surv = 0.88, 
                  BS.ppnRR = BS.ppnRR,
                  BS.invRR = BS.invRR)
  
  # mar.surv 
  # Withler et al. 2018 used 0.02 (0.01 and 0.05 in sensitivity analyses)
  # Bradford et al. 2023 used 0.03 
  marSurv0.01 <- lineSA (scenario = "Low survival",
                         h = sqrt(0.25), 
                         w = sqrt(100), 
                         mar.surv = 0.01, 
                         mar.surv.hatch = 0.0024, 
                         RS = 0.8, 
                         Theta.hatch = 80, 
                         p = 175, 
                         c = 400000,
                         bs.surv = 0.8, 
                         fec = 4900, 
                         release.surv = 0.88, 
                         BS.ppnRR = BS.ppnRR,
                         BS.invRR = BS.invRR)
  marSurv0.05 <- lineSA (scenario = "High survival",
                         h = sqrt(0.25), 
                         w = sqrt(100), 
                         mar.surv = 0.05, 
                         mar.surv.hatch = 0.0024, 
                         RS = 0.8, 
                         Theta.hatch = 80, 
                         p = 175, 
                         c = 400000,
                         bs.surv = 0.8, 
                         fec = 4900, 
                         release.surv = 0.88, 
                         BS.ppnRR = BS.ppnRR,
                         BS.invRR = BS.invRR)
  
  # mar.surv.hatch
  # Withler et al. 2018 used 0.0024 (0.001 and 0.005 in sensitivity analyses)
  # Bradford et al. 2023 used 0.005
  marSurvHatch0.001 <- lineSA (scenario = "Low hatchery survival",
                               h = sqrt(0.25), 
                               w = sqrt(100), 
                               mar.surv = 0.02, 
                               mar.surv.hatch = 0.001, 
                               RS = 0.8, 
                               Theta.hatch = 80, 
                               p = 175, 
                               c = 400000,
                               bs.surv = 0.8, 
                               fec = 4900, 
                               release.surv = 0.88, 
                               BS.ppnRR = BS.ppnRR,
                               BS.invRR = BS.invRR)
  marSurvHatch0.005 <- lineSA (scenario = "High hatchery survival",
                               h = sqrt(0.25), 
                               w = sqrt(100), 
                               mar.surv = 0.02, 
                               mar.surv.hatch = 0.005, 
                               RS = 0.8, 
                               Theta.hatch = 80, 
                               p = 175, 
                               c = 400000,
                               bs.surv = 0.8, 
                               fec = 4900, 
                               release.surv = 0.88, 
                               BS.ppnRR = BS.ppnRR,
                               BS.invRR = BS.invRR)
  # Fecundity (hard coded in run.lever.model). Withler used 4900
  # Bradford et al. 2023 use 4000
  fec4000 <- lineSA (scenario = "Low fecundity",
                     h = sqrt(0.25), 
                     w = sqrt(100), 
                     mar.surv = 0.02, 
                     mar.surv.hatch = 0.0024, 
                     RS = 0.8, 
                     Theta.hatch = 80, 
                     p = 175, 
                     c = 400000,
                     bs.surv = 0.8, 
                     fec = 4000, 
                     release.surv = 0.88, 
                     BS.ppnRR = BS.ppnRR,
                     BS.invRR = BS.invRR)
  # Survival from egg to release in hatchery (hard coded in run.lever.model)
  # Withler et al. 2018 used 0.88
  # Bradford et al . 2023 used 0.8
  releaseSurv0.8 <- lineSA (scenario = "Low release survival",
                            h = sqrt(0.25), 
                            w = sqrt(100), 
                            mar.surv = 0.02, 
                            mar.surv.hatch = 0.0024, 
                            RS = 0.8, 
                            Theta.hatch = 80, 
                            p = 175, 
                            c = 400000,
                            bs.surv = 0.8, 
                            fec = 4900, 
                            release.surv = 0.8, 
                            BS.ppnRR = BS.ppnRR,
                            BS.invRR = BS.invRR)
  
  #capacity
  # Witherl et al. 2018 used 400000
  # Bradford et al. 2023 use 100000
  capacity100000 <- lineSA (scenario = "Low capacity",
                            h = sqrt(0.25), 
                            w = sqrt(100), 
                            mar.surv = 0.02, 
                            mar.surv.hatch = 0.0024, 
                            RS = 0.8, 
                            Theta.hatch = 80, 
                            p = 175, 
                            c = 100000,
                            bs.surv = 0.8, 
                            fec = 4900, 
                            release.surv = 0.88, 
                            BS.ppnRR = BS.ppnRR,
                            BS.invRR = BS.invRR)
  
  # Pessimistic scenario
  pessimistic <- lineSA (scenario = "Pessimistic",
                         h = sqrt(0.25), 
                         w = sqrt(30), 
                         mar.surv = 0.02, 
                         mar.surv.hatch = 0.005, 
                         RS = 0.8, 
                         Theta.hatch = 60, 
                         p = 175, 
                         c = 400000,
                         bs.surv = 0.8, 
                         fec = 4900, 
                         release.surv = 0.88, 
                         BS.ppnRR = BS.ppnRR,
                         BS.invRR = BS.invRR)
  
  optimistic <- lineSA (scenario = "Optimistic",
                        h = sqrt(0.25), 
                        w = sqrt(100), 
                        mar.surv = 0.02, 
                        mar.surv.hatch = 0.001, 
                        RS = 0.2, 
                        Theta.hatch = 80, 
                        p = 175, 
                        c = 400000,
                        bs.surv = 0.8, 
                        fec = 4000, 
                        release.surv = 0.80, 
                        BS.ppnRR = BS.ppnRR,
                        BS.invRR = BS.invRR)
  out <- list(w30 = w30, 
              RS0.2 = RS0.2, 
              thetaHatch60 = thetaHatch60, 
              marSurvHatch0.001 = marSurvHatch0.001, 
              marSurvHatch0.005 = marSurvHatch0.005, 
              fec4000 = fec4000,
              releaseSurv0.8 = releaseSurv0.8, 
              pessimistic = pessimistic, 
              optimistic = optimistic, 
              baseCase = baseCase)
  return(out)
}


out <- run.scenarios(BS.ppnRR = BS.ppnRR, BS.invRR = BS.invRR)


# If BS.ppnRR = TRUE
if(BS.ppnRR == TRUE){
  df <- data.frame (data.frame(scenario = out$w30$scenario, PNI = out$w30$PNI, BSppnR = 1:100*(0.005))) |>
    # rbind (data.frame(scenario = w1000$scenario, PNI = w1000$PN, BSppnR = 1:100*(0.005))) |> # NO IMPACT REMOVE
    # rbind (data.frame(scenario = h0.05$scenario, PNI = h0.05$PNI, BSppnR = 1:100*(0.005))) |> # NO IMPACT REMOVE
    # rbind (data.frame(scenario = h0.5$scenario, PNI = h0.5$PNI, BSppnR = 1:100*(0.005))) |># NO IMPACT REMOVE
    rbind (data.frame(scenario = out$RS0.2$scenario, PNI = out$RS0.2$PNI, BSppnR = 1:100*(0.005))) |>
    rbind (data.frame(scenario = out$thetaHatch60$scenario, PNI = out$thetaHatch60$PNI, BSppnR = 1:100*(0.005))) |>
    # rbind (data.frame(scenario = p138$scenario, PNI = p138$PNI, BSppnR = 1:100*(0.005))) |>
    # rbind (data.frame(scenario = marSurv0.01$scenario, PNI = marSurv0.01$PNI, BSppnR = 1:100*(0.005))) |># NO IMPACT REMOVE
    # rbind (data.frame(scenario = marSurv0.05$scenario, PNI = marSurv0.05$PNI, BSppnR = 1:100*(0.005))) |># NO IMPACT REMOVE
    rbind (data.frame(scenario = out$marSurvHatch0.001$scenario, PNI = out$marSurvHatch0.001$PNI, BSppnR = 1:100*(0.005))) |>
    rbind (data.frame(scenario = out$marSurvHatch0.005$scenario, PNI = out$marSurvHatch0.005$PNI, BSppnR = 1:100*(0.005))) |>
    rbind (data.frame(scenario = out$fec4000$scenario, PNI = out$fec4000$PNI, BSppnR = 1:100*(0.005))) |>
    rbind (data.frame(scenario = out$releaseSurv0.8$scenario, PNI = out$releaseSurv0.8$PNI, BSppnR = 1:100*(0.005))) |>
    # rbind (data.frame(scenario = capacity100000$scenario, PNI = capacity100000$PNI, BSppnR = 1:100*(0.005)))
    rbind (data.frame(scenario = out$pessimistic$scenario, PNI = out$pessimistic$PNI, BSppnR = 1:100*(0.005))) |>
    rbind (data.frame(scenario = out$optimistic$scenario, PNI = out$optimistic$PNI, BSppnR = 1:100*(0.005))) |>
    rbind (data.frame(scenario = out$baseCase$scenario, PNI = out$baseCase$PNI, BSppnR = 1:100*(0.005))) 
  
}
  
# if (BS.invRR == TRUE){
#   df.permark <- data.frame (data.frame(scenario = out$w30$scenario, PNI = out$w30$PNI, perMark = 1:100*(0.01))) |>
#     # rbind (data.frame(scenario = w1000$scenario, PNI = w1000$PN, BSppnR = 1:100*(0.005))) |> # NO IMPACT REMOVE
#     # rbind (data.frame(scenario = h0.05$scenario, PNI = h0.05$PNI, BSppnR = 1:100*(0.005))) |> # NO IMPACT REMOVE
#     # rbind (data.frame(scenario = h0.5$scenario, PNI = h0.5$PNI, BSppnR = 1:100*(0.005))) |># NO IMPACT REMOVE
#     rbind (data.frame(scenario = out$RS0.2$scenario, PNI = out$RS0.2$PNI, perMark = 1:100*(0.01))) |>
#     rbind (data.frame(scenario = out$thetaHatch60$scenario, PNI = out$thetaHatch60$PNI, perMark = 1:100*(0.01))) |>
#     # rbind (data.frame(scenario = p138$scenario, PNI = p138$PNI, BSppnR = 1:100*(0.005))) |>
#     # rbind (data.frame(scenario = marSurv0.01$scenario, PNI = marSurv0.01$PNI, BSppnR = 1:100*(0.005))) |># NO IMPACT REMOVE
#     # rbind (data.frame(scenario = marSurv0.05$scenario, PNI = marSurv0.05$PNI, BSppnR = 1:100*(0.005))) |># NO IMPACT REMOVE
#     rbind (data.frame(scenario = out$marSurvHatch0.001$scenario, PNI = out$marSurvHatch0.001$PNI, perMark = 1:100*(0.01))) |>
#     rbind (data.frame(scenario = out$marSurvHatch0.005$scenario, PNI = out$marSurvHatch0.005$PNI, perMark = 1:100*(0.01))) |>
#     rbind (data.frame(scenario = out$fec4000$scenario, PNI = out$fec4000$PNI, perMark = 1:100*(0.01))) |>
#     rbind (data.frame(scenario = out$releaseSurv0.8$scenario, PNI = out$releaseSurv0.8$PNI, perMark = 1:100*(0.01))) |>
#     # rbind (data.frame(scenario = capacity100000$scenario, PNI = capacity100000$PNI, BSppnR = 1:100*(0.005)))
#     rbind (data.frame(scenario = out$pessimistic$scenario, PNI = out$pessimistic$PNI, perMark = 1:100*(0.01))) |>
#     rbind (data.frame(scenario = out$optimistic$scenario, PNI = out$optimistic$PNI, perMark = 1:100*(0.01))) |>
#     rbind (data.frame(scenario = out$baseCase$scenario, PNI = out$baseCase$PNI, perMark = 1:100*(0.01))) 
# }

nlevels <- 3
cp<-colorRampPalette(brewer.pal(9,'Blues'))
cp4<-rev(cp(nlevels+2)[-c(1,2)])
# cp4[c(1,nlevels/3, nlevels/1.5)]
#http://www.javascripter.net/faq/hextorgb.htm

if(BS.ppnRR){
  colLines <- ggplot(df, aes(x=BSppnR, y=PNI, colour=scenario)) +
    geom_line() + 
    theme_bw() +
    annotate('rect', xmin=0, xmax=0.55, ymin=0, ymax=0.49, alpha=.3, fill=cp4[1]) + 
    annotate('rect', xmin=0, xmax=0.55, ymin=0.5, ymax=0.71, alpha=.2, fill=cp4[2]) + 
    annotate('rect', xmin=0, xmax=0.55, ymin=0.72, ymax=1, alpha=.1, fill=cp4[3]) +
    xlab("Brood stock as a proportion of returns to river") +
    ylab("PNI") +
    scale_x_continuous(breaks = seq(0,0.5,0.1))
}

# if(BS.invRR){
#   colLines <- ggplot(df.permark, aes(x=perMark, y=PNI, colour=scenario)) +
#     geom_line() + 
#     theme_bw() +
#     annotate('rect', xmin=0, xmax=1.05, ymin=0, ymax=0.49, alpha=.3, fill=cp4[1]) + 
#     annotate('rect', xmin=0, xmax=1.05, ymin=0.5, ymax=0.71, alpha=.2, fill=cp4[2]) + 
#     annotate('rect', xmin=0, xmax=1.05, ymin=0.72, ymax=1, alpha=.1, fill=cp4[3]) +
#     xlab("Pecent of hatchery fish marked") + 
#     ylab("PNI") +
#     scale_x_continuous(breaks = seq(0,1,0.2))
#   
# }

# Change order of plotting scenarios
df$scenario <- factor(df$scenario,  c("High hatchery survival",
                                      "Low hatchery survival",
                                      "Large difference in optimal trait",
                                      "Low fecundity",
                                      "Low release survival",
                                      "Low RRS",
                                      "Strong selection",
                                      "Pessimistic",
                                      "Optimistic",
                                      "Base case") )
scenarios <- c(
  "High hatchery survival",
  "Low hatchery survival",
  "Large difference in optimal trait",
  "Low fecundity",
  "Low release survival",
  "Low RRS",
  "Strong selection",
  "Pessimistic",
  "Optimistic",
  "Base case")
pal <- setNames(c(rep("grey",9), "red"), scenarios)

if(BS.ppnRR){
  greyLines <- ggplot(df, aes(x=BSppnR, y=PNI, colour=scenario)) +
    geom_line() +
    theme_bw() +
    annotate('rect', xmin=0, xmax=0.55, ymin=0, ymax=0.49, alpha=.3, fill=cp4[1]) + 
    annotate('rect', xmin=0, xmax=0.55, ymin=0.5, ymax=0.71, alpha=.2, fill=cp4[2]) + 
    annotate('rect', xmin=0, xmax=0.55, ymin=0.72, ymax=1, alpha=.1, fill=cp4[3]) +
    xlab("Brood stock as a proportion of returns to river") + 
    ylab("PNI") +
    scale_colour_manual(values = pal, breaks= c( "High hatchery survival", 
                                                 "Low hatchery survival",
                                                 "Large difference in optimal trait",
                                                 "Low fecundity",
                                                 "Low release survival",
                                                 "Low RRS",
                                                 "Strong selection",
                                                 "Pessimistic",
                                                 "Optimistic",
                                                 "Base case") ) + 
    theme(legend.position="none")  +
    scale_x_continuous(breaks = seq(0,0.5,0.1))
}
# if(BS.invRR){
#   greyLines <- ggplot(df.permark, aes(x=perMark, y=PNI, colour=scenario)) +
#     geom_line() +
#     theme_bw() +
#     annotate('rect', xmin=0, xmax=1.05, ymin=0, ymax=0.49, alpha=.3, fill=cp4[1]) + 
#     annotate('rect', xmin=0, xmax=1.05, ymin=0.5, ymax=0.71, alpha=.2, fill=cp4[2]) + 
#     annotate('rect', xmin=0, xmax=1.05, ymin=0.72, ymax=1, alpha=.1, fill=cp4[3]) +
#     # xlab("Brood stock as a proportion of returns to river") + 
#     xlab("Percent of hatchery fish marked") + 
#     ylab("PNI") +
#     scale_colour_manual(values = pal, breaks= c( "High hatchery survival", 
#                                                  "Low hatchery survival",
#                                                  "Large difference in optimal trait",
#                                                  "Low fecundity",
#                                                  "Low release survival",
#                                                  "Low RRS",
#                                                  "Strong selection",
#                                                  "Pessimistic",
#                                                  "Optimistic",
#                                                  "Base case") ) + 
#     theme(legend.position="none")  +
#     scale_x_continuous(breaks = seq(0,1,0.2))
# }

# ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "colLines.png"), colLines, width=6, height=6, units="in")
  # ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "greyLines.png"), greyLines, width=5, height=6, units="in")
  # ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "colLinesPerMark1.png"), colLines, width=6, height=6, units="in")
  # ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "greyLinesPerMark1.png"), greyLines, width=5, height=6, units="in")
  # ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "colLinesPerMark1Sel0.5.png"), colLines, width=6, height=6, units="in")
  # ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "greyLinesPerMark1Sel0.5.png"), greyLines, width=5, height=6, units="in")

# ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "colLines_invRR.png"), colLines, width=6, height=6, units="in")
# ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "greyLines_invRR.png"), greyLines, width=5, height=6, units="in")
# ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "colLinesel0.5_invRR.png"), colLines, width=6, height=6, units="in")
# ggsave(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses" , "greyLinesSel0.5_invRR.png"), greyLines, width=5, height=6, units="in")
