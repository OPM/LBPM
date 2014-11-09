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
tau = 0.7
visc = (tau-0.5)/3

C1 <- read.csv(file="~/Data/Multiphase/Piston/Case1/timelog.tcat",head=TRUE,sep=" ")
C1$label <- "Case 1"
C1$U <- predict(sm.spline(C1$time,C1$sw),C1$time,1)*200

C2 <- read.csv(file="~/Data/Multiphase/Piston/Case2/timelog.tcat",head=TRUE,sep=" ")
C2$label <- "Case 2"
C2$U <- predict(sm.spline(C2$time,C2$sw),C2$time,1)*200

C3 <- read.csv(file="~/Data/Multiphase/Piston/Case3/timelog.tcat",head=TRUE,sep=" ")
C3$label <- "Case 3"
C3$U <- predict(sm.spline(C3$time,C3$sw),C3$time,1)*200

C4 <- read.csv(file="~/Data/Multiphase/Piston/Case4/timelog.tcat",head=TRUE,sep=" ")
C4$label <- "Case 4"
C4$U <- predict(sm.spline(C4$time,C4$sw),C4$time,1)*200

C5 <- read.csv(file="~/Data/Multiphase/Piston/Case5/timelog.tcat",head=TRUE,sep=" ")
C5$label <- "Case 5"
C5$U <- predict(sm.spline(C5$time,C5$sw),C5$time,1)*200

C6 <- read.csv(file="~/Data/Multiphase/Piston/Case6/timelog.tcat",head=TRUE,sep=" ")
C6$label <- "Case 6"
C6$U <- predict(sm.spline(C6$time,C6$sw),C6$time,1)*200

C7 <- read.csv(file="~/Data/Multiphase/Piston/Case1a/timelog.tcat",head=TRUE,sep=" ")
C7$label <- "Case 1a"
C7$U <- predict(sm.spline(C7$time,C7$sw),C7$time,1)*200

C8 <- read.csv(file="~/Data/Multiphase/Piston/Case2a/timelog.tcat",head=TRUE,sep=" ")
C8$label <- "Case 2a"
C8$U <- predict(sm.spline(C8$time,C8$sw),C8$time,1)*200

C9 <- read.csv(file="~/Data/Multiphase/Piston/Case3a/timelog.tcat",head=TRUE,sep=" ")
C9$label <- "Case 3a"
C9$U <- predict(sm.spline(C9$time,C9$sw),C9$time,1)*200

C10 <- read.csv(file="~/Data/Multiphase/Piston/Case4a/timelog.tcat",head=TRUE,sep=" ")
C10$label <- "Case 4a"
C10$U <- predict(sm.spline(C10$time,C10$sw),C10$time,1)*200

C11 <- read.csv(file="~/Data/Multiphase/Piston/Case5a/timelog.tcat",head=TRUE,sep=" ")
C11$label <- "Case 5a"
C11$U <- predict(sm.spline(C11$time,C11$sw),C11$time,1)*200

C12 <- read.csv(file="~/Data/Multiphase/Piston/Case6a/timelog.tcat",head=TRUE,sep=" ")
C12$label <- "Case 6a"
C12$U <- predict(sm.spline(C12$time,C12$sw),C12$time,1)*200

C13 <- read.csv(file="~/Data/Multiphase/Piston/Case7/timelog.tcat",head=TRUE,sep=" ")
C13$label <- "Case 7"
C13$U <- predict(sm.spline(C13$time,C13$sw),C13$time,1)*600

C14 <- read.csv(file="~/Data/Multiphase/Piston/Case8/timelog.tcat",head=TRUE,sep=" ")
C14$label <- "Case 8"
C14$U <- predict(sm.spline(C14$time,C14$sw),C14$time,1)*600

#Bind all of the data into one frame for ggplot (combine as labeled rows)
Full<-rbind(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12)
Big<-rbind(C13,C14)
 
Full$vawnz<-Full$vawnz*(-1)
Full$vawnsz<-Full$vawnsz*(-1)
Full$awn<-Full$awn*volume/(PI*D*D*D)
Full$sgkvpmawns<-Full$sgkvpmawns*(-1)
Full$trJwn<-D*Full$trJwn
Full$U<-Full$U*(-1)
Full$Ca<-Full$U*visc/gamma

Big$vawnz<-Big$vawnz*(-1)
Big$vawnsz<-Big$vawnsz*(-1)
Big$awn<-Big$awn*volume/(PI*D*D*D)
Big$sgkvpmawns<-Big$sgkvpmawns*(-1)
Big$trJwn<-D*Big$trJwn
Big$U<-Big$U*(-1)
Big$Ca<-Big$U*visc/gamma


# Plot the saturation against time
p<-ggplot(Full,aes(time,sw,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(sw) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/sw.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,awn,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewn) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/ewn.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,ans,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ens) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/ens.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,aws,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ews) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 

						      
ggsave("~/Data/Multiphase/Piston/ews.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,lwns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(ewns) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/ewns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/vaw.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vanz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vanz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/van.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawnz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/vawn.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,vawnsz,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(vawnsz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/vawns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,sgkvpmawns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab(cwns) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/cwns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Big,aes(time,KNwns,colour=label),coord_cartesian(xlim=c(0,25000)))+
	geom_line() + xlab("time") + ylab("KNwns") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 
						      
ggsave("~/Data/Multiphase/Piston/KNwns.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Full,aes(time,U,colour=label)) +
	geom_line() + xlab("time") + ylab("U") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(xlim=c(0,25000)) 

ggsave("~/Data/Multiphase/Piston/U.eps",p,width=4.5,heigh=3.5)

p<-ggplot(Big,aes(time,vawz,colour=label))+
	geom_line() + xlab("time") + ylab(vawz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(ylim=c(0,0.015)) 

ggsave("~/Data/Multiphase/Piston/big-vaw.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Big,aes(time,vanz,colour=label))+
	geom_line() + xlab("time") + ylab(vanz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(ylim=c(0,0.0185)) 

ggsave("~/Data/Multiphase/Piston/big-van.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Big,aes(time,vawnz,colour=label))+
	geom_line() + xlab("time") + ylab(vawnz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(ylim=c(0,0.015)) 


ggsave("~/Data/Multiphase/Piston/big-vawn.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Big,aes(time,vawnsz,colour=label))+
	geom_line() + xlab("time") + ylab(vawnsz) + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4)))  +
	coord_cartesian(ylim=c(0,0.015)) 


ggsave("~/Data/Multiphase/Piston/big-vawns.pdf",p,width=4.5,heigh=3.5)

p<-ggplot(Big,aes(time,U,colour=label))+
	geom_line() + xlab("time") + ylab("U") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=rel(4))) +
	coord_cartesian(ylim=c(0,0.015)) 

ggsave("~/Data/Multiphase/Piston/big-U.pdf",p,width=4.5,heigh=3.5)
