#Code to run contour plots over senstivity analyses
#Date last revised 5 Apr2017

library(akima)
library(RColorBrewer)
mar.surv.const<-TRUE#Is the marine suvival assumed to be 0.02 when estimating Seq, or derived from inputed marine survival?
source("R/run.lever.model.r")
source("R/ContourSA.r")

c<-400000#Beverton-Holt capacity parameter (Rmax)
percent.hatch<-0# 1-percent of hatchery fish that return to natural spawning grounds. Look at a range of values
HR<-0.4#Harvest rate
w<-sqrt(100)#sqrt(1000)
h<-sqrt(0.5)#sqrt(0.05)#sqrt(0.25)#sqrt(0.5)
RS<-0.2#0.8#
pdf.label<-"Contourplots6Nov2017w1000.pdf"


if (file.exists(here::here("Results")) == FALSE){
  dir.create(here::here("Results"))
  if (file.exists(here::here("Results", "2018")) == FALSE){
    dir.create(here::here("Results", "2018"))
    if (file.exists(here::here("Results", "2018", "MainResults")) == FALSE){
      dir.create(here::here("Results", "2018", "MainRestults"))
    }
    if (file.exists(here::here("Results", "2018", "SensitivityAnalyses")) == FALSE){
      dir.create(here::here("Results", "2018", "SensitivityAnalyses"))
    }
    if (file.exists(here::here("Results", "2018", "Appendices")) == FALSE){
      dir.create(here::here("Results", "2018", "Appendices"))
    }
  }
}

#Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), RS=0.8, pdf.label="Contourplots17Feb2017basecase.pdf")
# Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), RS=0.8, pdf.label="Contourplots21Feb2017_highheritability.pdf")
# Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(1000), RS=0.8, pdf.label="Contourplots21Feb2017_weakselection.pdf")
# Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), RS=0.5, pdf.label="Contourplots21Feb2017_lowhatchRS.pdf")

Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.25), w=sqrt(100), 
             RS=0.8, pdf.label="Contourplotsbasecase.pdf", scenario="", 
             title="", dir="Results/2018/MainResults")
Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.05), w=sqrt(100), 
             RS=0.8, pdf.label="Contourplots_lowheritability.pdf", 
             scenario="lowheritability", title="\nlow heritability", 
             dir="Results/2018/")
Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), 
             RS=0.8, pdf.label="Contourplots_highheritability.pdf", 
             scenario="highheritability", title="\nhigh heritability", 
             dir="Results/2018/")
# following runs ALSO have high heritabilty. These plots match with paper
Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(1000), 
             RS=0.8, pdf.label="Contourplots_weakselection.pdf", 
             scenario="weakselection", title="\nweak selection", 
             dir="Results/2018/")
Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(40), 
             RS=0.8, pdf.label="Contourplots_strongselection.pdf", 
             scenario="strongselection", title="\nstrong selection", 
             dir="Results/2018/")
Run.contours(c=400000, percent.hatch=0, HR=0.4, h=sqrt(0.5), w=sqrt(100), 
             RS=0.2, pdf.label="Contourplots_lowhatchRS.pdf", scenario="rs0.2", 
             title=paste(": low", gamma), dir="Results/2018/")
