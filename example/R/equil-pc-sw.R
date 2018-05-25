source("~/Programs/TCAT/R/DefsTCAT.R")
library(ggplot2)
library(extrafont)
library(reshape2)

IFT=0.058
D=66.933792
Data <- read.csv(file="~/Data/Multiphase/PorousMedia/Sph1964/equil.tcat",head=TRUE,sep=" ")
Data$pc <- (Data$pn-Data$pw)*D/IFT

p<-ggplot(Data,aes(sw,trJwn*D,colour=awn))+
geom_point()+ xlab(sw) + ylab(JwnD) + theme_bw() +
scale_colour_gradient(name=ewnD,limits=c(0,0.8),low="red")+
theme(text=element_text(family="Helvetica",size=rel(4))) +
theme(legend.text=element_text(size=rel(3))) +
theme(legend.title=element_text(size=rel(4)))

ggsave("Jwn-sw-awn.pdf",p,width=4.5,height=3.5);

p<-ggplot(Data,aes(sw,pc,colour=awn))+
geom_point()+ xlab(sw) + ylab(JwnD) + theme_bw() +
scale_colour_gradient(name=ewnD,limits=c(0,0.8),low="red")+
theme(text=element_text(family="Helvetica",size=rel(4))) +
theme(legend.text=element_text(size=rel(3))) +
theme(legend.title=element_text(size=rel(4)))

ggsave("pc-sw-awn.pdf",p,width=4.5,height=3.5);

p<-ggplot(Data,aes(sw,lwns))+
geom_point()+ xlab(sw) + ylab(ewn) + theme_bw() +
theme(text=element_text(family="Helvetica",size=rel(4)))

ggsave("ewns.pdf",p,width=4.5,height=3.5);

p<-ggplot(Data,aes(sw,awn,colour=trJwn))+
geom_point()+ xlab(sw) + ylab(ewn) + theme_bw() +
scale_colour_gradient(name=ewnD,limits=c(7,15),low="red")+
theme(text=element_text(family="Helvetica",size=rel(4)))

ggsave("awn-sw-pc.pdf",p,width=4.5,height=3.5);
