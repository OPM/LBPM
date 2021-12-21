require("ggplot2")

Discs<-read.csv("FullMicromodel.discs",head=FALSE,sep=" ")

colnames(Discs)<-c("cx","cy","radius")

L=0.45

# Extract some subset from the interior of the discs
SubDiscs<-subset(Discs,Discs$cy>0.9-L)
SubDiscs<-subset(SubDiscs,SubDiscs$cy<0.9+L)
SubDiscs<-subset(SubDiscs,SubDiscs$cx>0.9-L)
SubDiscs<-subset(SubDiscs,SubDiscs$cx<0.9+L)

SubDiscs$cx<-SubDiscs$cx-0.9+L
SubDiscs$cy<-SubDiscs$cy-0.9+L

write.table(SubDiscs,file="DiscPack.in",quote=FALSE,row.names=FALSE,col.names=FALSE,sep=" ")

ShowPlot<-ggplot(SubDiscs)+geom_circle(aes(x0=cx,y0=cy,r=radius))

