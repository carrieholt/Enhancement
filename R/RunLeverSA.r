#RunLeverSA
#Date last revised 6 Jan 2017

source(here::here("R", "FuncDefs.r"))
source(here::here("R", "run.lever.model.r"))

Data_DF <- data.frame(Per.mark=numeric(96), Hatchery.size=numeric(96), Sel=numeric(96), c=numeric(96), mar.surv=numeric(96), RS=numeric(96), mar.surv.hatch=numeric(96), HR=numeric(96), h=numeric(06), w=numeric(96), PNI=numeric(96),pHOS=numeric(96), pNOB=numeric(96), Rec=numeric(96), RperS=numeric(96), HP=numeric(96), C=numeric(96) )
Data_DF$Per.mark <- rep(rep(c(0.1,0.5), 4),12)
Data_DF$Hatchery.size <- rep(rep(rep(c(0.05, 0.3), each=2),2),12)
Data_DF$Sel <- rep(c(rep(0,4),rep(0.5,4)),12)
Data_DF$c <- c(rep(400000,8),rep(40000,8),rep(4000000,8),rep(400000,(96-3*8)))
Data_DF$mar.surv<-c(rep(0.02,3*8),rep(0.01,8), rep(0.1, 8), rep(0.02,(96-8*5)))
Data_DF$RS<-c(rep(0.8,5*8),rep(0.6,8), rep(0.8,(96-8*6)))
Data_DF$mar.surv.hatch<-c(rep(0.0024,6*8),rep(0.001,8), rep(0.005, 8), rep(0.0024,(96-8*8)))
Data_DF$HR<-c(rep(0.4,8*8),rep(0,8), rep(0.6,8), rep(0.4,(96-8*10)))
Data_DF$h<-c(rep(sqrt(0.25),10*8),rep(sqrt(0.5),8), rep(sqrt(0.25),(96-8*11)))
Data_DF$w<-c(rep(sqrt(100),11*8),rep(sqrt(1000),8))
                                  
for (i in 1:96){
  print(i)
  dum<-run.lever.model(per.mark=Data_DF$Per.mark[i],hatchery.size=Data_DF$Hatchery.size[i], sel=Data_DF$Sel[i], Theta.hatch=80, c=Data_DF$c[i], percent.hatch=0, HR=Data_DF$HR[i], h=Data_DF$h[i], w=Data_DF$w[i], mar.surv=Data_DF$mar.surv[i], RS=Data_DF$RS[i], mar.surv.hatch=Data_DF$mar.surv.hatch[i])#Add c, mar.surv, RS, mar.surv.hatch, 
  Data_DF$PNI[i]<-round(dum$PNI[100],2)
  Data_DF$pHOS[i]<-round(dum$pHOS[100],2)
  Data_DF$pNOB[i]<-round(dum$pNOB[100],2)                               
  Data_DF$Rec[i]<-round(dum$ret.nat.preharvest[100],0)
  Data_DF$RperS[i]<-round(dum$RperS[100],2)
  Data_DF$HP[i]<-round(dum$ret.hatch.preharvest[100],0)
  Data_DF$C[i]<-round(dum$catch[100],0)                                       
}

write.table(as.matrix(Data_DF), "SA6Jan2017.txt", sep="\t")

