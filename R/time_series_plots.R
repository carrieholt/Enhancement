# Time-series plots to evaluate plausible fitness parameters, i.e., 
# those that give long or short-term trajectoriss that agree with previous 
# literature
# Adapted from "HatcheryLevers.R", which gives similar trajectories limited to 
# 100 yrs

# Also used to investigate PNI and fitness outcomes of different brood take 
# limits as a proportion of returns to river
# See "plot.fig.SSR" for plots for CSAS SSR (2024)

# Date created: 5 Oct 2023
# Date last revised: 18 Jan 2024


library(ggplot2)
library(dplyr)
library(ggpubr)
# library(RColorBrewer)
library(viridis)
source("R/FuncDefs.r")
source("R/SeqFormulationBH.r")
source("R/run.lever.model.r")

c<-400000; # to match Lian's base case for May 2022
percent.hatch<-0 
HR<-0.2#0.4; # to match Lian's base case
h=sqrt(0.25)
w<-sqrt(100)
mar.surv.const<-TRUE

run_alone <- TRUE
if(run_alone){
  p = 175 # to match Lian's base case stage 0+
  per.mark= 0 #1.0,# to match Lian's base case
  hatchery.size=0.1 # not used
  BS.ppnRR = FALSE
  ppn.RR = 0.33
  BS.invRR = FALSE
  sel=0
  Theta.hatch=60 
  c=c
  sex.ratio=0.5 
  percent.hatch=percent.hatch
  HR=HR 
  mar.surv=0.02
  RS=0.8
  mar.surv.hatch=0.0024 
  bs.surv = 0.80  #survival from capture to spawning for broodstock
  fec = 4900 #fecundity
  release.surv = 0.88 #combined survival of hatchery fish from egg-fry and
  
}

res<-run.lever.model(p = 175, # to match Lian's base case stage 0+
                     per.mark= 0, #1.0,# to match Lian's base case
                     hatchery.size=0.1, # not used
                     BS.ppnRR = FALSE,
                     ppn.RR = 0.33,
                     sel=0, 
                     Theta.hatch=60, 
                     c=c, 
                     sex.ratio=0.5, 
                     percent.hatch=percent.hatch, 
                     HR=HR, 
                     h=h, 
                     w=w, 
                     mar.surv=0.02, 
                     RS=0.8, 
                     mar.surv.hatch=0.0024)# to match Lian's base case stage 0+

res.nogenetics <- run.lever.model(p = 175, # to match Lian's base case
                                  per.mark=0, #1.0,# to match Lian's base case
                                  hatchery.size=0.1, 
                                  BS.ppnRR = FALSE,
                                  ppn.RR = 0.33,
                                  sel=0, 
                                  Theta.hatch=100, 
                                  c=c, 
                                  sex.ratio=0.5, 
                                  percent.hatch=percent.hatch, 
                                  HR=HR, 
                                  h=h, 
                                  w=w, 
                                  mar.surv=0.02, 
                                  RS=0.8, 
                                  mar.surv.hatch=0.0024)#0.0024)# to match Lian's base case

time_series_plots<-function(res=res, res.nogenetics=res.nogenetics, ngen=100,  BS.ppnRR = FALSE, ppn.RR = ppn.RR){
  if(ngen == 11) {time_frame <- 0:10}
  if(ngen != 11) {time_frame <- 1:ngen}
  par(mfrow=c(3,2), mar=c(4,4,2,1))
  plot(time_frame,res$fit.smolt[1:ngen]^2, type="l", ylim=c(0,1), ylab="Fitness", xlab="Generation")
  #points(x=1:100,y=res$P.nat/100, col="red",type="l")
  #points(x=1:100,y=res$P.hatch/100, col="blue", type="l")
  #abline(h=res$Theta.nat/100, col="red", lty="dotted")
  #abline(h=res$Theta.hatch/100, col="blue", lty="dotted")
  #legend(x=5,y=0.40, legend=c("Population fitness", "Phenotype: natural origin fish", "Phenotype: hatchery origin fish"), col=c("black", "red", "blue"), cex=0.7, lty=c("solid", "solid"), bty="n")
  text(x=1, y=0.93, labels="(a)")
  # title<-paste("Mark rate=",res$per.mark*100,"%;  Brood stock=",round(res$BS[1],0), " (",res$hatchery.size*100, "% ave natural recruitment);  Removal of marked fish=",res$sel*100,"%", sep="")
  if (!BS.ppnRR) title<-paste("Theta.hatch=",res$Theta.hatch,";  selection pressure, w^2=", w^2, "; Brood stock=", round(res$BS[1],0), " (",res$hatchery.size*100, "% ave natural recruitment)", sep="")
  if(BS.ppnRR) title<-paste("Theta.hatch=",res$Theta.hatch,";  selection pressure, w^2=", w^2, "; Brood stock= ", ppn.RR, " of returns to river", sep="")
  title<-paste("Rule on max ppn of returns=",res$BS.ppnRR.cap,"; Brood stock cap=", round(res$BS[100],0), " (",res$hatchery.size*100, "% ave natural recruitment)", sep="")
  mtext(title, side=3, line=0.5, at=ngen, cex=0.75)
  
  y.upper.lim<-max(res$Sp.nat, na.rm=T)#res$Seq*res$mar.surv*40
  if(max(res$Sp.nat)>res$Seq*1.2){y.upper.lim<-max(res$Sp.nat)}
  
  plot(time_frame, res$Sp.nat[1:ngen], type="l", col="red", ylim=c(0,y.upper.lim), ylab="Spawners (HOS+NOS)", xlab="Generation")
  points(x=time_frame, y=res$NOS[1:ngen], col="red", type="l", lty="dashed")
  points(x=time_frame, y=res$HOS[1:ngen], col="red", type="l",lty="dotted")
  points(x=time_frame, y=res$BS[1:ngen], col="blue", type="l")
  points(x=time_frame, y=res$NOB[1:ngen], col="blue", type="l",lty="dashed")
  points(x=time_frame, y=res$HOB[1:ngen], col="blue", type="l",lty="dotted")
  legend(x=5,y=y.upper.lim*0.83, legend=c("Natural spawners", "NOS", "HOS", "Brood stock removals", "NOB (accounting for mortality prior to spawning)", "HOB"), col=c(rep("red",3),rep("blue",3)), cex=0.7, lty=rep(c("solid", "dashed","dotted"),2), bty="n")
  text(x=1, y=y.upper.lim*0.95, labels="(b)")
  
  # y.max <- max(c(res$NOS, res$NOB), na.rm=T)
  # plot(time_frame, res$NOS[1:ngen], type="l", col="red",lty="dashed", ylim=c(0,y.max), ylab="HOS, NOS(scaled from (b))", xlab="Generation")
  # points(x=time_frame, y=res$NOB[1:ngen], col="blue", type="l",lty="dashed")
  # legend(x=5,y=y.max*0.83, legend=c("NOS","NOB (accounting for mortality prior to spawning)"), col=c(rep("red",1),rep("blue",1)), cex=0.7, lty=rep(c("dashed"),2), bty="n")
  # text(x=1, y=y.max*0.95, labels="(c)")
  # 
  plot(time_frame,res$HOS+res$HOB, type="l", col="red", ylab="Returns", xlab="Generation", ylim = c(0,max(c(res$HOS+res$HOB, res$NOS+res$NOB), na.rm=T) ) )
  points(x=1:ngen, y=(res$NOB + res$NOS), col="blue",type="l")
  legend(x=ngen/2,y=max(res$HOS+res$HOB)*0.85, legend=c("Hatchery returns", "Natural returns"), col=c(rep("red",1),rep("blue",1)), cex=0.7, lty=rep(c("solid"),2), bty="n")
  
  plot(time_frame, res$pHOSeff[1:ngen], type="l", col="red", ylim=c(0,1), ylab="Proportion (pHOSeff or pNOB)", xlab="Generation")
  lines(time_frame, res$pNOB[1:ngen], col="blue")
  lines(time_frame, res$PNI[1:ngen], col="black")
  legend(x=5,y=0.5, legend=c("pHOSeff","pNOB", "PNI"), cex=0.7, col=c("red", "blue", "black"), lty=c("solid", "solid", "solid"), bty="n")
  text(x=1, y=0.93, labels="(d)")
  
  # plot(time_frame, res$PNI[1:ngen], type="l", col="black", ylim=c(0,1), ylab="PNI",  xlab="Generation")
  # text(x=1, y=0.95, labels="(d)")
  
  plot(time_frame, res$RperS[1:ngen], ylim=c(min(res$RperS)*0.8, max(res$RperS.hatch)*1.2), type="l", col="black",  ylab="Recruits/spawner on spawning grounds", xlab="Generation")
  lines(time_frame, res$RperS.hatch[1:ngen], col="blue")
  text(x=1, y=max(res$RperS)*1.18, labels="(e)")
  # plot(time_frame, res$Sm.nat[1:ngen],  type="l", col="black",  ylab="Natural smolt production", xlab="Generation")
  # text(x=1, y=max(res$Sm.nat)*1.18, labels="(e)")
  
  y.upper.lim<-max(c(res$ret.nat.preharvest, res.nogenetics$ret.nat.preharvest), na.rm=T)#res$c*res$mar.surv
  #if(max(res$ret.nat)>res$Seq*1.2){y.upper.lim<-max(res$ret.nat)*1.1}
  plot(time_frame, res$ret.nat.preharvest[1:ngen], type="l", col="red", ylim=c(0,y.upper.lim), ylab="Recruitment from HOS+NOS", xlab="Generation")
  points(x=time_frame, res.nogenetics$ret.nat.preharvest[1:ngen], type="l", col="blue")
  legend(x=5,y=y.upper.lim*0.8, legend=c("Without genetic impacts",  "With genetic impacts"), cex=0.7, col=c("blue", "red"), lty=c("solid", "solid"), bty="n")
  text(x=1, y=y.upper.lim*0.95, labels="(f)")
  
}

if (file.exists(here::here("Results")) == FALSE){
  dir.create(here::here("Results"))
}
if (file.exists(here::here("Results/FitnessPars")) == FALSE){
  dir.create(here::here("Results/FitnessPars"))
}

return_time_series_plots<-function(res=res, ngen=100,  BS.ppnRR = FALSE, ppn.RR = ppn.RR){
  if(ngen == 11) {time_frame <- 0:10}
  if(ngen != 11) {time_frame <- 1:ngen}
  par(mfrow=c(1,1))#, mar=c(4,4,2,1))
  
  plot(time_frame,res$HOS+res$HOB, type="l", col="red", ylab="Returns", xlab="Generation", ylim = c(0,max(res$HOS+res$HOB, res$NOS+res$NOB)))
  points(x=1:ngen, y=(res$NOB + res$NOS), col="blue",type="l")
  legend(x=ngen/2,y=max(res$HOS+res$HOB)*0.85, legend=c("Hatchery returns", "Natural returns"), col=c(rep("red",1),rep("blue",1)), cex=0.7, lty=rep(c("solid"),2), bty="n")
  #points(x=1:100,y=res$P.hatch/100, col="blue", type="l")
  #abline(h=res$Theta.nat/100, col="red", lty="dotted")
  #abline(h=res$Theta.hatch/100, col="blue", lty="dotted")
  #legend(x=5,y=0.40, legend=c("Population fitness", "Phenotype: natural origin fish", "Phenotype: hatchery origin fish"), col=c("black", "red", "blue"), cex=0.7, lty=c("solid", "solid"), bty="n")
  # title<-paste("Mark rate=",res$per.mark*100,"%;  Brood stock=",round(res$BS[1],0), " (",res$hatchery.size*100, "% ave natural recruitment);  Removal of marked fish=",res$sel*100,"%", sep="")
  if (!BS.ppnRR) title<-paste("Theta.hatch=",res$Theta.hatch,";  selection pressure, w^2=", w^2, "; Brood stock=", round(res$BS[1],0), " (",res$hatchery.size*100, "% ave natural recruitment)", sep="")
  if(BS.ppnRR) title<-paste("Theta.hatch=",res$Theta.hatch,";  selection pressure, w^2=", w^2, "; Brood stock= ", ppn.RR, " of returns to river", sep="")
  mtext(title, side=3, line=0.5, at=ngen/2, cex=0.75)
 
}


# pdf(here::here("Results/FitnessPars/Timeseries-w2_40_hatchsize_0.3_sexRatio1_PnatpHOS.pdf"))
# pdf(here::here("Results/FitnessPars/Timeseries-w2_40_hatchsize_0.3_sexRatio1.pdf"))
# ***Need to change ngen to 500 or 11 in run.lever.model.r and in time-series function above
# pdf(here::here("Results/FitnessPars/Timeseries-w2_100_hatchsize_0.1_500gen.pdf"))
# pdf(here::here("Results/FitnessPars/Timeseries-w2_100_hatchsize_0.1_100gen.pdf"))
hs <- 0.1
w <- sqrt(100)#sqrt(40)
#Theta.hatch = 80
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                       BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
                       Theta.hatc = 80, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = w,
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                                  BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = w, mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024)
time_series_plots(res, res.nogenetics)
#Theta.hatch = 70
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                       BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
                       Theta.hatc = 70, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = w,
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                                  BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = w, mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024)
time_series_plots(res, res.nogenetics)
#Theta.hatch = 60
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                       BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
                       Theta.hatc = 60, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = w,
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                                  BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = w, mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024)
time_series_plots(res, res.nogenetics)
#Theta.hatch = 50
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                       BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
                       Theta.hatc = 50, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = w,
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                                  BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = w, mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024)
time_series_plots(res, res.nogenetics)
# #Theta.hatch = 40
# res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
#                        BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
#                        Theta.hatc = 40, c = c, sex.ratio = 0.5,
#                        percent.hatch = percent.hatch, HR = HR, h = h, w = w,
#                        mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
# res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
#                                   BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
#                                   Theta.hatch = 100, c = c, sex.ratio = 0.5, 
#                                   percent.hatch = percent.hatch, HR = HR, h = h,
#                                   w = w, mar.surv = 0.02, RS = 0.8, 
#                                   mar.surv.hatch = 0.0024)
# time_series_plots(res, res.nogenetics)
# 
# #Theta.hatch = 20
# res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
#                        BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
#                        Theta.hatc = 20, c = c, sex.ratio = 0.5,
#                        percent.hatch = percent.hatch, HR = HR, h = h, w = w,
#                        mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
# res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
#                                   BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
#                                   Theta.hatch = 100, c = c, sex.ratio = 0.5, 
#                                   percent.hatch = percent.hatch, HR = HR, h = h,
#                                   w = w, mar.surv = 0.02, RS = 0.8, 
#                                   mar.surv.hatch = 0.0024)
# time_series_plots(res, res.nogenetics)
# #Theta.hatch = 0
# res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
#                        BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
#                        Theta.hatc = 0, c = c, sex.ratio = 0.5,
#                        percent.hatch = percent.hatch, HR = HR, h = h, w = w,
#                        mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024)
# res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
#                                   BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
#                                   Theta.hatch = 100, c = c, sex.ratio = 0.5, 
#                                   percent.hatch = percent.hatch, HR = HR, h = h,
#                                   w = w, mar.surv = 0.02, RS = 0.8, 
#                                   mar.surv.hatch = 0.0024)
# time_series_plots(res, res.nogenetics)

# dev.off()

# run for w2=100, theta.hatch=0, 20, 40, 60, 80.
# run for w2= 40, theta.hatch=0, 20, 40, 60, 80.
# run for w2= 40, theta.hatch=0, 20, 40, 60, 80, sex.ratio = 1
# run for w2= 40, theta.hatch=0, 20, 40, 60, 80,sex.ratio = 1, Pnat calculated from pHOS instead of pHOSeff


# PLots of time-series for various ppns of returns to river showing dyanmics prior to equilibrium
# pdf(here::here("Results/Timeseries-ppnRR-withRperS-noGenetics.pdf"))
ppn.RR <-  0.1
Theta.hatch <- 100#80
hatchery.size <- 0.1
BS.invRR <- FALSE
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                       BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0,
                       Theta.hatch = Theta.hatch, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = sqrt(100),
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                                  BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = sqrt(100), mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
time_series_plots(res, res.nogenetics, BS.ppnRR = TRUE, ppn.RR=ppn.RR)
# return_time_series_plots(res,  BS.ppnRR = TRUE, ppn.RR=ppn.RR)

ppn.RR <-  0.2
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                       BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0,
                       Theta.hatch = Theta.hatch, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = sqrt(100),
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                                  BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = sqrt(100), mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
time_series_plots(res, res.nogenetics, BS.ppnRR = TRUE, ppn.RR=ppn.RR)
# return_time_series_plots(res,  BS.ppnRR = TRUE, ppn.RR=ppn.RR)

ppn.RR <-  0.3
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                       BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0,
                       Theta.hatch = Theta.hatch, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = sqrt(100),
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                                  BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = sqrt(100), mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
time_series_plots(res, res.nogenetics, BS.ppnRR = TRUE, ppn.RR=ppn.RR)
# return_time_series_plots(res,  BS.ppnRR = TRUE, ppn.RR=ppn.RR)

ppn.RR <-  0.4
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                       BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0,
                       Theta.hatch = Theta.hatch, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = sqrt(100),
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                                  BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = sqrt(100), mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
time_series_plots(res, res.nogenetics, BS.ppnRR = TRUE, ppn.RR=ppn.RR)
# return_time_series_plots(res,  BS.ppnRR = TRUE, ppn.RR=ppn.RR)

ppn.RR <-  0.5
res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                       BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0,
                       Theta.hatch = Theta.hatch, c = c, sex.ratio = 0.5,
                       percent.hatch = percent.hatch, HR = HR, h = h, w = sqrt(100),
                       mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hatchery.size, 
                                  BS.ppnRR = TRUE, ppn.RR = ppn.RR, sel = 0, 
                                  Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                  percent.hatch = percent.hatch, HR = HR, h = h,
                                  w = sqrt(100), mar.surv = 0.02, RS = 0.8, 
                                  mar.surv.hatch = 0.0024, BS.invRR = BS.invRR)
time_series_plots(res, res.nogenetics, BS.ppnRR = TRUE, ppn.RR=ppn.RR)
# return_time_series_plots(res,  BS.ppnRR = TRUE, ppn.RR=ppn.RR)

# dev.off()

# Results on time-series were BS is determined as a ppn of eq pop size, but 
# there is a max cap on BS at beginning as ppn RR, BS.ppnRR.cap

plot.timeseries <- FALSE
if(plot.timeseries){
  pdf(here::here("Results/Timeseries-BScaps_hatchsize_0.1.pdf"))
  hs <- 0.1
  w <- sqrt(100)#sqrt(40)
  S.init.ppn <- 0.1
  # BS.ppnRR.cap = 0.1
  for (i in 1:10){
    BS.ppnRR.cap <- i*0.05
    print(BS.ppnRR.cap)
    res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                           BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
                           Theta.hatc = 80, c = c, sex.ratio = 0.5,
                           percent.hatch = percent.hatch, HR = HR, h = h, w = w,
                           mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, 
                           BS.ppnRR.cap = BS.ppnRR.cap, S.init.ppn = S.init.ppn)
    res.nogenetics <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                                      BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0, 
                                      Theta.hatch = 100, c = c, sex.ratio = 0.5, 
                                      percent.hatch = percent.hatch, HR = HR, h = h,
                                      w = w, mar.surv = 0.02, RS = 0.8, 
                                      mar.surv.hatch = 0.0024, BS.ppnRR.cap = BS.ppnRR.cap, 
                                      S.init.ppn = S.init.ppn)
    time_series_plots(res, res.nogenetics)
    
  }
  dev.off()
  
}

# Suggested Fig 1 from Tim

plot.fig.SSR <- TRUE # Figures for CSAS SSR (2024)
if(plot.SSR){
  
  
  
  # hs <- 0.1
  for (hs in seq(0.1,0.5,0.1)){
    w <- sqrt(100)#sqrt(40)
    S.init.ppn <- 0.1
    # BS.ppnRR.cap = 0.1
    PNI.dum <- matrix(NA, nrow=100, ncol=3)
    fitness.dum <- matrix(NA, nrow=100, ncol=3)
    bs.dum <- matrix(NA, nrow=100, ncol=3)
    PNI.df <- NA
    fitness.df <- NA
    bs.df <- NA
    
    for (i in 1:10){
      BS.ppnRR.cap <- i*0.05 # actually limit not cap
      print(BS.ppnRR.cap)
      res <- run.lever.model(p = 175, per.mark = 0, hatchery.size = hs, 
                             BS.ppnRR = FALSE, ppn.RR = 0.33, sel = 0,
                             Theta.hatc = 80, c = c, sex.ratio = 0.5,
                             percent.hatch = percent.hatch, HR = HR, h = h, w = w,
                             mar.surv = 0.02, RS = 0.8, mar.surv.hatch = 0.0024, 
                             BS.ppnRR.cap = BS.ppnRR.cap, S.init.ppn = S.init.ppn)
      
      # Extract PNI from model runs
      PNI.dum[,1] <- BS.ppnRR.cap  # actually limit not cap
      PNI.dum[,2] <- 1:100
      PNI.dum[,3] <- res$PNI
      if (i==1) PNI.df <- PNI.dum
      if (i>1) PNI.df <- rbind(PNI.df, PNI.dum)
      
      # Extract life-cycle fitness from model runs
      fitness.dum[,1] <- BS.ppnRR.cap  # actually limit not cap
      fitness.dum[,2] <- 1:100
      fitness.dum[,3] <- res$fit.smolt^2
      if (i==1) fitness.df <- fitness.dum
      if (i>1) fitness.df <- rbind(fitness.df, fitness.dum)
      
      # Extarct brood size from model runs
      bs.dum[,1] <- BS.ppnRR.cap
      bs.dum[,2] <- 1:100
      bs.dum[,3] <- res$BS
      if (i==1) bs.df <- bs.dum
      if (i>1) bs.df <- rbind(bs.df, bs.dum)
      
    }
    
    PNI.df <- as.data.frame(PNI.df)
    PNI.df <- PNI.df |> dplyr::transmute(Limit = V1, Year = V2, PNI=V3)
    fitness.df <- as.data.frame(fitness.df)
    fitness.df <- fitness.df |> 
      dplyr::transmute(Limit = V1, Year = V2, Fitness = V3)
    bs.df <- as.data.frame(bs.df)
    bs.df <- bs.df |> dplyr::transmute(Limit = V1, Year = V2, MaxHatcherySize=V3)
    
    
    # nb.cols <- 18
    # # mycolours <- colorRampPalette(brewer.pal(8, "GnBu"))(nb.cols)[9:18]
    # mycolours <- colorRampPalette(brewer.pal(8, "GnBu"))(nb.cols)[9:18]
    
    PNI.plot <- PNI.df |> ggplot(aes(Year, PNI, colour = as.factor(Limit))) + 
      geom_line(size=1.1) + 
      scale_color_viridis(discrete = TRUE, direction = -1, option="viridis") +
      # theme_classic() + 
      theme_minimal() +
      theme(legend.position = "none") +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14)) + 
      annotate("text", x = 1, y = 1.03, label = "(a)") +
      ylim(0,1.03) 
    
    fitness.plot <- fitness.df |> ggplot(aes(Year, Fitness, 
                                             colour = as.factor(Limit))) + 
      geom_line(size=1.1) + 
      scale_color_viridis(discrete = TRUE, direction = -1, option="viridis") + 
      # scale_color_manual(values = mycolours) +
      # theme_classic() + 
      theme_minimal() +
      theme(legend.position = "none") +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14)) + 
      annotate("text", x = 1, y = 1.03, label = "(b)") +
      ylim(0,1.03) 
    legend.plot <- PNI.df |> ggplot(aes(Year, PNI, colour = as.factor(Limit))) + 
      geom_line(size=1.1) + 
      scale_color_viridis(discrete = TRUE, direction = -1, option="viridis", 
                          name="Limit on brood take\nas a proportion of\nreturns") 
    leg <- get_legend(legend.plot)
    legend <- as_ggplot(leg) 
    
    multi.panel <- ggarrange(PNI.plot, fitness.plot, legend, nrow = 1, ncol = 3, 
                             widths = c(4,4,1.5))
    #--------------------------------------------------------------------------
    # Figure 1- printed 5 times for each cap on hatchery size
    file.name <- paste("Fig1-HatcheryCap-", hs, sep="")
    # ggsave( paste(here::here("Results"), "/", file.name, ".png", sep = ""), 
    #         multi.panel, width = 10, height =5, bg="white")
    #--------------------------------------------------------------------------
    
    # Data for Fig. 2
    fitness.df.last <- fitness.df |> filter(Year==100) |> select(!Year) |> 
      mutate (Cap = hs)
    if(hs == 0.1) df.1 <- fitness.df.last
    if(hs > 0.1) df.1 <- rbind(df.1, fitness.df.last)
    
    bs.df.last <- bs.df |> filter(Year==100) |> select(!Year) |> 
      mutate(Cap = hs)
    if(hs==0.1) df.2 <- bs.df.last
    if(hs > 0.1) df.2 <- rbind(df.2, bs.df.last)
   
  }

#--------------------------------------------------------------------------
# Figure 2
# df.1 |> ggplot(aes(Limit, Fitness, colour = as.factor(Cap))) +
#     geom_line(size = 1.1)
  
  Fig2.df <- left_join(df.1,df.2, join_by(Limit == Limit, Cap == Cap)) 
  
  Fig2a <- Fig2.df |> ggplot(aes(Limit, Fitness, colour = as.factor(Cap))) +
    geom_line(size = 1.1) + 
    geom_point(size = 2, shape = 21, fill="white") +
    scale_color_viridis(discrete = TRUE, direction = -1, option="viridis") + 
    theme_minimal() +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) + 
    annotate("text", x = 0.01, y = 1.03, label = "(a)") +
    labs(x="Limit on brood size as\na proportion of returns")
    
  Fig2b <- Fig2.df |> ggplot(aes(MaxHatcherySize, Fitness, colour = as.factor(Cap))) +
    geom_line(size = 1.1) + 
    geom_point(size = 2, shape = 21, fill="white") +
    scale_color_viridis(discrete = TRUE, direction = -1, option="viridis") + 
    theme_minimal() +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) + 
    annotate("text", x = 1, y = 1.03, label = "(b)") + 
    labs(x="Maximum\nhatchery size", y=element_blank())
    
  Fig2.legend <- Fig2.df |> ggplot(aes(MaxHatcherySize, Fitness, 
                                       colour = as.factor(Cap))) +
    geom_line(size = 1.1) + 
    geom_point(size = 2, shape = 21, fill="white") +
    theme_minimal() +
    scale_color_viridis(discrete = TRUE, direction = -1, option="viridis", 
                        name="Cap on brood\nsize as a\nproportion of\nequilibrium\nabundances") 

  Fig2.leg <- get_legend(Fig2.legend)
  legend <- as_ggplot(Fig2.leg) 
  
  Fig2.multi.panel <- ggarrange(Fig2a, Fig2b, legend, nrow = 1, ncol = 3, 
                           widths = c(4, 3.7, 1.5))
  
  Fig2.multi.panel

  ggsave( here::here("Results/Fig2.png"), 
           Fig2.multi.panel, width = 10, height =5, bg="white")
  
}
