# Code to run contour plots over original sensitivity analyses in Withler et al. 
# 2018, but limiting brood take to be less than specified ppns of returns to 
# river
# Date last revised 24 Feb 2023

library(akima)
library(RColorBrewer)
mar.surv.const<-TRUE#Is the marine suvival assumed to be 0.02 when estimating Seq, or derived from inputed marine survival?
source("R/run.lever.model.r")
source("R/ContourSA.r")

# c<-400000#Beverton-Holt capacity parameter (Rmax)
# percent.hatch<-0# 1-percent of hatchery fish that return to natural spawning grounds. Look at a range of values
# HR<-0.4#Harvest rate
# w<-sqrt(100)#sqrt(1000)
# h<-sqrt(0.25)#sqrt(0.05)#sqrt(0.25)#sqrt(0.5)
# RS<-0.8#0.2
# sex.ratio <- 0.5
# ppn.RR <- 0.6


if (file.exists(here::here("Results")) == FALSE){
  dir.create(here::here("Results"))
}

if (file.exists(here::here("Results", "ppnReturnsRiver")) == FALSE){
  dir.create(here::here("Results", "ppnReturnsRiver"))
}
if (file.exists(here::here("Results", "ppnReturnsRiver")) == TRUE){
  if (file.exists(here::here("Results", "ppnReturnsRiver", "MainResults")) == FALSE){
    dir.create(here::here("Results", "ppnReturnsRiver", "MainResults"))
  }
  if (file.exists(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses")) == FALSE){
    dir.create(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses"))
  }
  if (file.exists(here::here("Results", "ppnReturnsRiver", "Appendices")) == FALSE){
    dir.create(here::here("Results", "ppnReturnsRiver", "Appendices"))
  }
}




Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.1, HR=0.4,
             h=sqrt(0.25), w=sqrt(100), RS=0.8,
             scenario="10percent", title="\n10% of returns=limit to BS", dir="Results/ppnReturnsRiver/MainResults/")

Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.3, HR=0.4,
             h=sqrt(0.25), w=sqrt(100), RS=0.8,
             scenario="30percent", title="\n30% of returns=limit to BS", dir="Results/ppnReturnsRiver/MainResults/")

Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.6, HR=0.4,
             h=sqrt(0.25), w=sqrt(100), RS=0.8,
             scenario="60percent", title="\n60% of returns=limit to BS", dir="Results/ppnReturnsRiver/MainResults/")
