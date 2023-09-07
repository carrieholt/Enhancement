#***************************************************************************
#Data for Contour plots for sensitivity analyses on management levers
#Date last revised 5 July 2023
#Creates data.frames of results for 3 combinations of management levers (3 rows), and 4 combinations of marine survival (high and low for natural/hatchery pop), for 3x4 plot
# Used by ContourSA_ms_plots.r

#First run.lever.model over all combinations, and output 4 lists (of lists)
#***************************************************************************


# Need to install older versin of akima for contour plots
# detach("package:akima", unload=TRUE)
# install.packages("C:/Users/HoltC/AppData/Local/Programs/R/R-4.2.1/library/akima_0.6-2.tar.gz", repos = NULL, type="source")

library(akima)
library(RColorBrewer)

mar.surv.const<-TRUE#Is the marine suvival assumed to be 0.02 when estimating Seq, or derived from inputed marine survival?
c=400000; percent.hatch=0; HR=0.4; h=sqrt(0.25); w=sqrt(100);RS=0.8; 
sex.ratio=0.5
p <- 175

mar.surv.vec=c(0.01,0.05,0.02,0.02); mar.surv.hatch.vec=c(0.0024,0.0024,0.001, 0.005)
for (k in 1:4){
  mar.surv<-mar.surv.vec[k]; mar.surv.hatch<-mar.surv.hatch.vec[k]
  PNImh<-matrix(NA,100,100); pHOSmh<-matrix(NA,100,100); pNOBmh<-matrix(NA,100,100);  RperSmh<-matrix(NA,100,100); Ret.natmh<-matrix(NA,100,100); Ret.hatmh<-matrix(NA,100,100); Catchmh<-matrix(NA,100,100); BSmarkmh<-matrix(NA,100,100);
  PNIsh<-matrix(NA,100,100); pHOSsh<-matrix(NA,100,100); pNOBsh<-matrix(NA,100,100);  RperSsh<-matrix(NA,100,100); Ret.natsh<-matrix(NA,100,100); Ret.hatsh<-matrix(NA,100,100); Catchsh<-matrix(NA,100,100); BSmarksh<-matrix(NA,100,100);
  PNIms<-matrix(NA,100,100); pHOSms<-matrix(NA,100,100); pNOBms<-matrix(NA,100,100);  RperSms<-matrix(NA,100,100); Ret.natms<-matrix(NA,100,100); Ret.hatms<-matrix(NA,100,100); Catchms<-matrix(NA,100,100); BSmarkms<-matrix(NA,100,100);
  BSmarkmh<-matrix(NA,100,100)
  HLmarkmh<-matrix(NA,100,100)
  BSmarksh<-matrix(NA,100,100)
  HLmarksh<-matrix(NA,100,100)
  BSmarkms<-matrix(NA,100,100)
  HLmarkms<-matrix(NA,100,100)
  
  
  mh.hatchery.size<-NA; mh.sel<-NA; mh.BS<-NA
  sh.hatchery.size<-NA; sh.per.mark<-NA; sh.BS<-NA
  ms.hatchery.size<-NA; ms.BS<-NA
  
  for(i in 1:100){
    print(i)
    for (j in 1:100){
      mh<-run.lever.model(per.mark=i*0.01, hatchery.size=0.15, ppn.RR=j*0.005, BS.ppnRR = TRUE, sel=0, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, RS=RS, mar.surv=mar.surv, mar.surv.hatch=mar.surv.hatch, sex.ratio=sex.ratio, p=p)
      mh.sel<-mh$sel
      if (i==1){mh.hatchery.size[j]<-mh$hatchery.size; mh.BS[j]<-mh$BS[100]; mh.Seq <- mh$Seq}
      PNImh[i,j]<-mh$PNI[100]
      pHOSmh[i,j]<-mh$pHOSeff[100]
      pNOBmh[i,j]<-mh$pNOB[100]
      RperSmh[i,j]<-mh$RperS[100]
      Ret.natmh[i,j]<-mh$ret.nat.preharvest[100]
      Ret.hatmh[i,j]<-mh$ret.hatch.preharvest[100]
      Catchmh[i,j]<-mh$catch[100]
      BSmarkmh[i,j]<-mh$BS.mark
      HLmarkmh[i,j]<-mh$handling.limit.mark
      
      sh<-run.lever.model(per.mark=0.5, hatchery.size=0.15, ppn.RR=j*0.005, BS.ppnRR = TRUE, sel=i*0.01, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, RS=RS, mar.surv=mar.surv, mar.surv.hatch=mar.surv.hatch, sex.ratio=sex.ratio, p=p)
      sh.per.mark<-sh$per.mark
      if (i==1){sh.hatchery.size[j]<-sh$hatchery.size;  sh.BS[j]<-sh$BS[100]; sh.Seq <- sh$Seq}
      PNIsh[i,j]<-sh$PNI[100]
      pHOSsh[i,j]<-sh$pHOSeff[100]
      pNOBsh[i,j]<-sh$pNOB[100]
      RperSsh[i,j]<-sh$RperS[100]  
      Ret.natsh[i,j]<-sh$ret.nat.preharvest[100]  
      Ret.hatsh[i,j]<-sh$ret.hatch.preharvest[100]  
      Catchsh[i,j]<-sh$catch[100]
      BSmarksh[i,j]<-sh$BS.mark
      HLmarksh[i,j]<-sh$handling.limit.mark
      
      ms<-run.lever.model(per.mark=i*0.01, hatchery.size=0.15, ppn.RR=0.3, BS.ppnRR = TRUE, sel=j*0.01, Theta.hatch=80, c=c, percent.hatch=percent.hatch, HR=HR, h=h, w=w, RS=RS, mar.surv=mar.surv, mar.surv.hatch=mar.surv.hatch, sex.ratio=sex.ratio, p=p)
      ms.hatchery.size<-ms$hatchery.size
      ms.BS<ms$BS[100]
      PNIms[i,j]<-ms$PNI[100]
      pHOSms[i,j]<-ms$pHOSeff[100]
      pNOBms[i,j]<-ms$pNOB[100]
      RperSms[i,j]<-ms$RperS[100]
      Ret.natms[i,j]<-ms$ret.nat.preharvest[100]
      Ret.hatms[i,j]<-ms$ret.hatch.preharvest[100]
      Catchms[i,j]<-ms$catch[100]
      BSmarkms[i,j]<-ms$BS.mark
      HLmarkms[i,j]<-ms$handling.limit.mark
      
    }
  }
  if(k==1){lowms<-list( mh.sel=mh.sel, mh.hatchery.size=mh.hatchery.size, mh.BS=mh.BS, mh.Seq=mh.Seq, sh.per.mark=sh.per.mark, sh.hatchery.size=sh.hatchery.size, sh.BS=sh.BS, sh.Seq=sh.Seq, ms.hatchery.size=ms.hatchery.size, ms.BS=ms.BS, PNImh=PNImh, PNIsh=PNIsh, PNIms=PNIms, Ret.natmh=Ret.natmh, Ret.natsh=Ret.natsh, Ret.natms=Ret.natms, Ret.hatmh=Ret.hatmh, Ret.hatsh=Ret.hatsh, Ret.hatms=Ret.hatms, pNOBmh=pNOBmh, pNOBsh=pNOBsh, pNOBms=pNOBms, pHOSmh=pHOSmh, pHOSsh=pHOSsh, pHOSms=pHOSms, BSmarkmh=BSmarkmh, BSmarksh=BSmarksh, BSmarkms=BSmarkms,  HLmarkmh= HLmarkmh,  HLmarksh=HLmarksh,  HLmarkms= HLmarkms) } 
  if(k==2){highms<-list( mh.sel=mh.sel, mh.hatchery.size=mh.hatchery.size, mh.BS=mh.BS, mh.Seq=mh.Seq, sh.per.mark=sh.per.mark, sh.hatchery.size=sh.hatchery.size, sh.BS=sh.BS, sh.Seq=sh.Seq, ms.hatchery.size=ms.hatchery.size, ms.BS=ms.BS, PNImh=PNImh, PNIsh=PNIsh, PNIms=PNIms, Ret.natmh=Ret.natmh, Ret.natsh=Ret.natsh, Ret.natms=Ret.natms, Ret.hatmh=Ret.hatmh, Ret.hatsh=Ret.hatsh, Ret.hatmsh=Ret.hatms, pNOBmh=pNOBmh, pNOBsh=pNOBsh, pNOBms=pNOBms, pHOSmh=pHOSmh, pHOSsh=pHOSsh, pHOSms=pHOSms, BSmarkmh=BSmarkmh, BSmarksh=BSmarksh, BSmarkms=BSmarkms,  HLmarkmh= HLmarkmh,  HLmarksh=HLmarksh,  HLmarkms= HLmarkms)}  
  if(k==3){lowmshatch<-list( mh.sel=mh.sel, mh.hatchery.size=mh.hatchery.size, mh.BS=mh.BS, mh.Seq=mh.Seq, sh.per.mark=sh.per.mark, sh.hatchery.size=sh.hatchery.size, sh.BS=sh.BS, sh.Seq=sh.Seq, ms.hatchery.size=ms.hatchery.size, ms.BS=ms.BS, PNImh=PNImh, PNIsh=PNIsh, PNIms=PNIms, Ret.natmh=Ret.natmh, Ret.natsh=Ret.natsh, Ret.natms=Ret.natms, Ret.hatmh=Ret.hatmh, Ret.hatsh=Ret.hatsh, Ret.hatms=Ret.hatms, pNOBmh=pNOBmh, pNOBsh=pNOBsh, pNOBms=pNOBms, pHOSmh=pHOSmh, pHOSsh=pHOSsh, pHOSms=pHOSms, BSmarkmh=BSmarkmh, BSmarksh=BSmarksh, BSmarkms=BSmarkms,  HLmarkmh= HLmarkmh,  HLmarksh=HLmarksh,  HLmarkms= HLmarkms)}  
  if(k==4){highmshatch<-list( mh.sel=mh.sel, mh.hatchery.size=mh.hatchery.size, mh.BS=mh.BS, mh.Seq=mh.Seq, sh.per.mark=sh.per.mark, sh.hatchery.size=sh.hatchery.size, sh.BS=sh.BS, sh.Seq=sh.Seq, ms.hatchery.size=ms.hatchery.size, ms.BS=ms.BS, PNImh=PNImh, PNIsh=PNIsh, PNIms=PNIms, Ret.natmh=Ret.natmh, Ret.natsh=Ret.natsh, Ret.natms=Ret.natms, Ret.hatmh=Ret.hatmh, Ret.hatsh=Ret.hatsh, Ret.hatms=Ret.hatms, pNOBmh=pNOBmh, pNOBsh=pNOBsh, pNOBms=pNOBms, pHOSmh=pHOSmh, pHOSsh=pHOSsh, pHOSms=pHOSms, BSmarkmh=BSmarkmh, BSmarksh=BSmarksh, BSmarkms=BSmarkms,  HLmarkmh= HLmarkmh,  HLmarksh=HLmarksh,  HLmarkms= HLmarkms)}  
}#of of for (k in 1:4)



