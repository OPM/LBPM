library(ggplot2)
library(extrafont)
library(reshape2)
library(pspline)

# define the TCAT variables, length scale and IFT
source("~/Programs/TCAT/R/DefsTCAT.R")
gamma = 0.058
D=15.5
volume=2e6
PI=3.14159265

tau = 1.0
visc = (tau-0.5)/3
C1 <- read.csv(file="~/Data/Multiphase/Piston/Case1/timelog.tcat",head=TRUE,sep=" ")
C1$label <- "A"
C1$U <- predict(sm.spline(C1$time,C1$sw),C1$time,1)*400
C1$Ca<-C1$U*visc/gamma

tau = 0.9
visc = (tau-0.5)/3
C2<-read.csv(file="~/Data/Multiphase/Piston/Case2/timelog.tcat",head=TRUE,sep=" ")
C2$label <- "B"
C2$U <- predict(sm.spline(C2$time,C2$sw),C2$time,1)*400
C2$Ca<-C2$U*visc/gamma

tau = 0.8
visc = (tau-0.5)/3
C3 <- read.csv(file="~/Data/Multiphase/Piston/Case3/timelog.tcat",head=TRUE,sep=" ")
C3$label <- "C"
C3$U <- predict(sm.spline(C3$time,C3$sw),C3$time,1)*400
C3$Ca<-C3$U*visc/gamma

tau = 1.0
visc = (tau-0.5)/3
C4 <- read.csv(file="~/Data/Multiphase/Piston/Case4/timelog.tcat",head=TRUE,sep=" ")
C4$label <- "D"
C4$U <- predict(sm.spline(C4$time,C4$sw),C4$time,1)*400
C4$Ca<-C4$U*visc/gamma

tau = 0.9
visc = (tau-0.5)/3
C5 <- read.csv(file="~/Data/Multiphase/Piston/Case5/timelog.tcat",head=TRUE,sep=" ")
C5$label <- "E"
C5$U <- predict(sm.spline(C5$time,C5$sw),C5$time,1)*400
C5$Ca<-C5$U*visc/gamma

tau = 0.8
visc = (tau-0.5)/3
C6 <- read.csv(file="~/Data/Multiphase/Piston/Case6/timelog.tcat",head=TRUE,sep=" ")
C6$label <- "F"
C6$U <- predict(sm.spline(C6$time,C6$sw),C6$time,1)*400
C6$Ca<-C6$U*visc/gamma

#Bind all of the data into one frame for ggplot (combine as labeled rows)
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

# Plot the saturation against time
p<-ggplot(Full,aes(time,sw,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(sw) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-sw.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,awn,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewn) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-ewn.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,ans,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ens) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-ens.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,aws,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ews) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-ews.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,lwns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewns) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-ewns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-vaw.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vanz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vanz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))
  
ggsave("~/Data/Multiphase/Piston/piston-van.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawnz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-vawn.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawnsz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnsz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-vawns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,sgkvpmawns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(cwns) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-cwns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,KNwns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(KNwns) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-KNwns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,KGwns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(KGwns) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))						      
ggsave("~/Data/Multiphase/Piston/piston-KGwns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,U,colour=label)) +
	geom_line() + xlab("time") + ylab("U") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))
ggsave("~/Data/Multiphase/Piston/piston-U.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,Gwnxx,colour=label))+
	geom_line() + xlab("time") + ylab(Gwnxx) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))
						      
ggsave("~/Data/Multiphase/Piston/piston-Gwnxx.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,Ca,colour=label))+
	geom_line() + xlab("time") + ylab("Ca") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))
						      
ggsave("~/Data/Multiphase/Piston/piston-Ca.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,trJwn,colour=label))+
	geom_line() + xlab("time") + ylab(JwnD) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	theme(legend.text=element_text(size=rel(4)))
						      
ggsave("~/Data/Multiphase/Piston/piston-Jwn.eps",p,width=4.5,heigh=3.5)
