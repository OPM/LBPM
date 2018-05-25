library(ggplot2)
library(extrafont)
library(reshape2)
library(pspline)

setwd("~/Data/Multiphase/Piston/")

mainfont="Helvetica"

# additional definitions
dswdt=expression(L~frac(ds^{bar(bar(w))},dt))
dawndt=expression(frac(d~epsilon^{bar(bar(wn))},dt))
dansdt=expression(frac(V,2~pi~R)~frac(d~epsilon^{bar(bar(ws))},dt))
LHSA=expression(frac((p^{n}-p^{w})~R,gamma^{wn})-J[w]^{wn}~R)
RHSAa=expression(frac(eta~R,gamma^{wn})~frac(d~s^{bar(bar(w))},dt))
RESA=expression(frac((p^{n}-p^{w})~R,gamma^{wn})-J[w]^{wn}~R-c^{wn}~frac(eta~R,gamma^{wn})~frac(d~s^{bar(bar(w))},dt))
RHSAb=expression(frac(gamma^{wn}~(epsilon^{bar(bar(wn))}~R-epsilon[eq]^{bar(bar(wn))}~R),p^{wn}~R))
RHSB=expression(frac(1,2)~G[zz]^{wn}~(~v[z]^{bar(w)}~+~v[z]^{bar(n)}))

#pression(frac(d~s^{bar(bar(w))},dt^{*}))

# define the TCAT variables, length scale and IFT
source("~/Programs/LBPM-WIA/example/R/DefsTCAT.R")
gamma = 0.058
D=14.5
volume=4e6
PI=3.14159265
Jeq=1.711

tau = 0.7
visc = (tau-0.5)/3
C1 <- read.csv(file="~/Data/Multiphase/Piston/Case1/timelog.tcat",head=TRUE,sep=" ")
C1$Case <- "A"
C1$U <- predict(sm.spline(C1$time,C1$sw),C1$time,1)*400
C1$dswdt <- predict(sm.spline(C1$time,C1$sw),C1$time,1)
C1$Ca<-C1$U*visc/gamma
C1$dawndt <-predict(sm.spline(C1$time,C1$awn),C1$time,1)
C1$dansdt <-predict(sm.spline(C1$time,C1$ans),C1$time,1)*volume/(2*PI*D)
C1$pc<-(C1$pn-C1$pw)*D/gamma

tau = 0.7
visc = (tau-0.5)/3
C2<-read.csv(file="~/Data/Multiphase/Piston/Case2/timelog.tcat",head=TRUE,sep=" ")
C2$Case <- "B"
C2$U <- predict(sm.spline(C2$time,C2$sw),C2$time,1)*400
C2$dswdt <- predict(sm.spline(C2$time,C2$sw),C2$time,1)
C2$Ca<-C2$U*visc/gamma
C2$dawndt <-predict(sm.spline(C2$time,C2$awn),C2$time,1)
C2$dansdt <-predict(sm.spline(C2$time,C2$ans),C2$time,1)*volume/(2*PI*D)
C2$pc<-(C2$pn-C2$pw)*D/gamma

tau = 0.7
visc = (tau-0.5)/3
C3 <- read.csv(file="~/Data/Multiphase/Piston/Case3/timelog.tcat",head=TRUE,sep=" ")
C3$Case <- "C"
C3$U <- predict(sm.spline(C3$time,C3$sw),C3$time,1)*400
C3$dswdt <- predict(sm.spline(C3$time,C3$sw),C3$time,1)
C3$Ca<-C3$U*visc/gamma
C3$dawndt <-predict(sm.spline(C3$time,C3$awn),C3$time,1)
C3$dansdt <-predict(sm.spline(C3$time,C3$ans),C3$time,1)*volume/(2*PI*D)
C3$pc<-(C3$pn-C3$pw)*D/gamma

tau = 1.0
visc = (tau-0.5)/3
C4 <- read.csv(file="~/Data/Multiphase/Piston/Case4/timelog.tcat",head=TRUE,sep=" ")
C4$Case <- "D"
C4$U <- predict(sm.spline(C4$time,C4$sw),C4$time,1)*400
C4$dswdt <- predict(sm.spline(C4$time,C4$sw),C4$time,1)
C4$Ca<-C4$U*visc/gamma
C4$dawndt <-predict(sm.spline(C4$time,C4$awn),C4$time,1)
C4$dansdt <-predict(sm.spline(C4$time,C4$ans),C4$time,1)*volume/(2*PI*D)
C4$pc<-(C4$pn-C4$pw)*D/gamma

tau = 1.0
visc = (tau-0.5)/3
C5 <- read.csv(file="~/Data/Multiphase/Piston/Case5/timelog.tcat",head=TRUE,sep=" ")
C5$Case <- "E"
C5$U <- predict(sm.spline(C5$time,C5$sw),C5$time,1)*400
C5$dswdt <- predict(sm.spline(C5$time,C5$sw),C5$time,1)
C5$Ca<-C5$U*visc/gamma
C5$dawndt <-predict(sm.spline(C5$time,C5$awn),C5$time,1)
C5$dansdt <-predict(sm.spline(C5$time,C5$ans),C5$time,1)*volume/(2*PI*D)
C5$pc<-(C5$pn-C5$pw)*D/gamma

tau = 1.0
visc = (tau-0.5)/3
C6 <- read.csv(file="~/Data/Multiphase/Piston/Case6/timelog.tcat",head=TRUE,sep=" ")
C6$Case <- "F"
C6$U <- predict(sm.spline(C6$time,C6$sw),C6$time,1)*400
C6$dswdt <- predict(sm.spline(C6$time,C6$sw),C6$time,1)
C6$Ca<-C6$U*visc/gamma
C6$dawndt<-predict(sm.spline(C6$time,C6$awn),C6$time,1)
C6$dansdt<-predict(sm.spline(C6$time,C6$ans),C6$time,1)*volume/(2*PI*D)
C6$pc<-(C6$pn-C6$pw)*D/gamma

#Bind all of the data into one frame for ggplot (combine as Caseed rows)
Full<-rbind(C1,C2,C3,C4,C5,C6)

Full$Ca<-Full$Ca*(-1) 
Full$vawnz<-Full$vawnz*(-1)
Full$vawnsz<-Full$vawnsz*(-1)
#Full$awn<-Full$awn*volume/(PI*D*D*D)
Full$sgkvpmawns<-Full$sgkvpmawns*(-1)
Full$trJwn<-D*Full$trJwn
Full$U<-Full$U*(-1)
Full$awn<-Full$awn/D
Full$aws<-Full$aws/D
Full$ans<-Full$ans/D
Full$lwns<-Full$lwns/D/D
#Full$RHSAeq<-Full$pc*400/D-Jeq
Full$RHSA<-Full$pc*400/D-Full$trJwn*400/D

Full$LHSeq<-Full$pc-Full$trJwn
Full$LHS<-Full$pc-Full$trJwn
Full$NonDimA<-Full$Ca*D/400
Full$NonDimB<-(Full$awn*D-872.95/volume)/Full$pc
Full$ResA<-Full$LHS-2808.1693*Full$NonDimA+0.1038

myfitA<-lm(NonDimA~LHS+NonDimB,data=Full)
summary.lm(myfitA)
Err <- summary(myfitA)$coefficients[1,1]
cwn <- summary(myfitA)$coefficients[2,1]
kwn <- summary(myfitA)$coefficients[3,1]

Full$RHSB<-0.5*Full$Gwnzz*(Full$vawz+Full$vanz)

#postscript(family="ComputerModern",encoding="TeXtext.enc")
#loadfonts(device="postscript")
#font_install('fontcm')


p<-ggplot(Full,aes(LHS,NonDimA,colour=Case))+
	geom_point() + ylab(expression(frac(ds^{bar(bar(w))},dt~paste("*")))) + xlab(C^{wn}~Delta~P^{c}) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))+
	geom_abline(slope=cwn,intercept=Err,colour="gray")
						      
ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,NonDimA,colour=Case))+
	geom_line() + xlab("time") + ylab(expression(frac(ds^{bar(bar(w))},dt~paste("*")))) + 
	theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))
 
ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim-A.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,coefA1+coefA2*LHS,colour=Case))+
	geom_line() + xlab("time") + 
	ylab(expression(c^{wn~paste("*")}~(p^{wn~paste("*")}-p^{c~paste("*")})+epsilon^{~paste("*")})) + 
	theme_bw() + theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))
 
ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim-B.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,coefA3*NonDimB,colour=Case))+
	geom_line() + xlab("time") + 
	ylab(expression(Delta~epsilon^{bar(bar(wn))~paste("*")}/p^{c~paste("*")})) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))
  
ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim-C.pdf",p,width=4.5,heigh=3.5)


Steady<-subset(Full,time>37000 & time<39000,select=c(RHSA,Ca,vawnz,RHSB))
fitSteady<-lm(Ca~RHSA,data=Steady)

# Plot the saturation against time
p<-ggplot(Full,aes(time,sw,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(sw) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-sw.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,awn,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewn) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-ewn.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,ans,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ens) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-ens.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,aws,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ews) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-ews.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,lwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-ewns.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0075)) +
	theme(legend.title=element_text(size=rel(3)))
					      
ggsave("~/Data/Multiphase/Piston/piston-vaw.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vanz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vanz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
  
ggsave("~/Data/Multiphase/Piston/piston-van.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawnz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-vawn.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawnsz,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnsz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-vawns.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,sgkvpmawns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(cwns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-cwns.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,KNwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(KNwns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-KNwns.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,KGwns,colour=Case),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(KGwns) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-KGwns.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,U,colour=Case)) +
	geom_line() + xlab("time") + ylab(dswdt) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0075)) +
	theme(legend.title=element_text(size=rel(3)))

ggsave("~/Data/Multiphase/Piston/piston-dswdt.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,Gwnxx,colour=Case))+
	geom_line() + xlab("time") + ylab(Gwnxx) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-Gwnxx.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,Gwnzz,colour=Case))+
	geom_line() + xlab("time") + ylab(Gwnzz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-Gwnzz.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,Ca,colour=Case))+
	geom_line() + xlab("time") + ylab("Ca") + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-Ca.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,trJwn,colour=Case))+
	geom_line() + xlab("time") + ylab(JwnD) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))

						      
ggsave("~/Data/Multiphase/Piston/piston-Jwn.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,dansdt,colour=Case))+
	geom_line() + xlab("time") + ylab(dansdt) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0075)) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-dansdt.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,dawndt,colour=Case))+
	geom_line() + xlab("time") + ylab(dawndt) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(4))) +
	theme(legend.text=element_text(size=rel(4))) +
	theme(legend.title=element_text(size=rel(3)))
						      
ggsave("~/Data/Multiphase/Piston/piston-dawndt.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(NonDimA,LHS,colour=Case))+
	geom_point() + xlab(RHSAa) + ylab(LHSA) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3))) 
#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim-A.pdf",p,width=4.5,heigh=3.5)

#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim-A-eq.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(NonDimA,LHSeq,colour=Case))+
	geom_point() + xlab(expression(frac(ds^{bar(bar(w))},dt~paste("*")))) + ylab(expression(p^{paste(wn,"*")}-p^{paste(c,"*")})) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))+
	geom_abline(slope=2917.1011,intercept=-0.1527,colour="gray")
#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(NonDimB,LHS,colour=Case))+
	geom_point() + xlab(RHSAb) + ylab(LHSA) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3))) 
#+
#	coord_cartesian(ylim=c(0,2))

#ggsave("~/Data/Multiphase/Piston/piston-dswdt-nondim-B.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(NonDimB,ResA,colour=Case))+
	geom_point() + xlab(RHSAb) + ylab(RESA) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))

#ggsave("~/Data/Multiphase/Piston/piston-dswdt-residual-A.pdf",p,width=4.5,heigh=3.5)

#p<-ggplot(subset(Full,time>37000 & time<39000),aes(RHSB,vawnz,colour=Case))+
p<-ggplot(Full,aes(RHSB,vawnz,colour=Case))+
	geom_point() + xlab(RHSB) + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family=mainfont,size=rel(3))) +
	theme(legend.text=element_text(size=rel(4)))+
	theme(legend.title=element_text(size=rel(3)))+
	geom_abline(slope=1,intercept=0,colour="gray")
#	stat_smooth(aes(group=1),method=lm,fullrange=TRUE,se=FALSE,colour="gray")

ggsave("~/Data/Multiphase/Piston/piston-Gwn-vawn.pdf",p,width=4.5,heigh=3.5)

#p<-ggplot(SSA,aes(Ca,RHSA))+
#	geom_point() + xlab(LHSA) + ylab(RHSA) + theme_bw() +
#	theme(text=element_text(family=mainfont,size=rel(3))) +
#	theme(legend.text=element_text(size=rel(4))) +
#	coord_cartesian(xlim=c(-0,52)) +
#	theme(legend.title=element_text(size=rel(3)))
						      
