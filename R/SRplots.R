#SR plots
#14FEb2017
source("FuncDefs.r")
par(mfrow=c(1,1))
HR<-0.4
RS<-1
p<-175#Beverton-Holt productivity parameter
c<-400000#Beverton-Holt capacity parameter (Rmax)
mar.surv<-0.01# Marine survival of natural-origin fish
per.hatch<-0#All hatchery fish return to natural spawning grounds

plot((1:1000)*100,mar.surv*BH(0, (1:1000)*100, RS, p, c)*(1-HR), xlim=c(0,10000), ylim=c(0,10000), xlab="Spawners", ylab="Adult recruits")

points((1:1000)*100,mar.surv*5*BH(0, (1:1000)*100, RS, p, c)*(1-HR), col="red")

lines(c(0,10000),c(0,10000))
