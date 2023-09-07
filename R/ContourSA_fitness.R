#Code to run contour plots over senstivity analyses
#Date last revised 23 Feb 2023

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
# ppn.RR <- 0.33
# pdf.label<-"Contourplots6Nov2017w1000.pdf"


if (file.exists(here::here("Results")) == FALSE){
  dir.create(here::here("Results"))
  }
if (file.exists(here::here("Results", "2018")) == FALSE){
  dir.create(here::here("Results", "2018"))
  }
if (file.exists(here::here("Results", "2018", "MainResults")) == FALSE){
  dir.create(here::here("Results", "2018", "MainResults"))
  }
if (file.exists(here::here("Results", "2018", "SensitivityAnalyses")) == FALSE){
  dir.create(here::here("Results", "2018", "SensitivityAnalyses"))
  }
if (file.exists(here::here("Results", "2018", "Appendices")) == FALSE){
  dir.create(here::here("Results", "2018", "Appendices"))
  }

  
if (file.exists(here::here("Results", "UpdatedSexRatio")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatio"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatio", "MainResults")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatio", "MainResults"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatio", "SensitivityAnalyses")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatio", "SensitivityAnalyses"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatio", "Appendices")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatio", "Appendices"))
  }
  
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeff")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeff"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeff", "MainResults")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeff", "MainResults"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeff", "SensitivityAnalyses")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeff", "SensitivityAnalyses"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeff", "Appendices")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeff", "Appendices"))
  }
  
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh", "MainResults")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh", "MainResults"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh", "SensitivityAnalyses")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh", "SensitivityAnalyses"))
  }
if (file.exists(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh", "Appendices")) == FALSE){
  dir.create(here::here("Results", "UpdatedSexRatiopHOSeffPNIThresh", "Appendices"))
  }



if (file.exists(here::here("Results", "ppnReturnsRiver")) == FALSE){
  dir.create(here::here("Results", "ppnReturnsRiver"))
  }
if (file.exists(here::here("Results", "ppnReturnsRiver", "MainResults")) == FALSE){
  dir.create(here::here("Results", "ppnReturnsRiver", "MainResults"))
  }
if (file.exists(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses")) == FALSE){
  dir.create(here::here("Results", "ppnReturnsRiver", "SensitivityAnalyses"))
  }
if (file.exists(here::here("Results", "ppnReturnsRiver", "Appendices")) == FALSE){
  dir.create(here::here("Results", "ppnReturnsRiver", "Appendices"))
  }



Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.33, HR=0.4, 
             h=sqrt(0.25), w=sqrt(100), RS=0.8,  
             scenario="", title="", dir="Results/UpdatedSexRatiopHOSeffPNIThresh/MainResults/")
Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.33, HR=0.4, 
             h=sqrt(0.05), w=sqrt(100), RS=0.8, 
             scenario="LowHeritability",  title="\nlow heritability", 
             dir="Results/UpdatedSexRatiopHOSeffPNIThresh/")
Run.contours(c=400000, percent.hatch=0,  sex.ratio = 0.5, ppn.RR = 0.33, HR=0.4, 
             h=sqrt(0.5), w=sqrt(100), RS=0.8, 
             scenario="HighHeritability", title="\nhigh heritability", 
             dir="Results/UpdatedSexRatiopHOSeffPNIThresh/")
# following runs also have base-case heritability. These plots match with paper
# ResDoc used high heritability for the following 3 sensitivity analyses, h=sqrt(5)
Run.contours(c=400000, percent.hatch=0,  sex.ratio = 0.5, ppn.RR = 0.33, HR=0.4, 
             h=sqrt(0.25), w=sqrt(1000), RS=0.8, 
             scenario="WeakSelection", title="\nweak selection", 
             dir="Results/UpdatedSexRatiopHOSeffPNIThresh/")
Run.contours(c=400000, percent.hatch=0,  sex.ratio = 0.5, ppn.RR = 0.33, HR=0.4, 
             h=sqrt(0.25), w=sqrt(40), RS=0.8, 
             scenario="StrongSelection", title="\nstrong selection", 
             dir="Results/UpdatedSexRatiopHOSeffPNIThresh/")
Run.contours(c=400000, percent.hatch=0,  sex.ratio = 0.5, ppn.RR = 0.33, HR=0.4, 
             h=sqrt(0.25), w=sqrt(100), RS=0.2,  scenario="rs0.2", 
             title=paste(": low", gamma), dir="Results/UpdatedSexRatiopHOSeffPNIThresh/")

# # 
# Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.3, HR=0.4,
#              h=sqrt(0.25), w=sqrt(100), RS=0.8,
#              scenario="30percent", title="\n30% of returns=limit to BS", dir="Results/ppnReturnsRiver/MainResults/")
# 
# Run.contours(c=400000, percent.hatch=0, sex.ratio = 0.5, ppn.RR = 0.6, HR=0.4,
#              h=sqrt(0.25), w=sqrt(100), RS=0.8,
#              scenario="60percent", title="\n60% of returns=limit to BS", dir="Results/ppnReturnsRiver/MainResults/")
