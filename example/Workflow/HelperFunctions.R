require("ggplot2")

IngestRelperm<-function(PATH){

  FILE=paste0(PATH,"/subphase.csv")
  S<-read.csv(FILE,head=TRUE,sep=" ")
  S$Vw<-S$Vwc+S$Vwd
  S$Vn<-S$Vnc+S$Vnd
  S$Sw<-S$Vw/(S$Vn+S$Vw)
  S$Qwx<-S$Vw*(S$Pwc_x+S$Pwd_x)/(S$Mwc+S$Mwd)
  S$Qnx<-S$Vn*(S$Pnc_x+S$Pnd_x)/(S$Mnc+S$Mnd)
  S$Krn<-S$nun*S$Qnx/S$Fx
  S$Krw<-S$nuw*S$Qwx/S$Fx
  S$Case<-PATH
  return(S)
}

ReadRelperm<-function(PATH){
	FILE=paste0(PATH,"/relperm.csv")
	D<-read.csv(file=FILE,head=FALSE,sep=" ")
	colnames(D)<-c("time","nun","nuw","ift","Fx","Fy","Fz","Vn","Vw","unx","uny","unz","uwx","uwy","uwz")
	D$Sw<-D$Vw/(D$Vn+D$Vw)
	D$Krn<-D$Vn*D$nun*D$unx/D$Fx
	D$Krw<-D$Vw*D$nuw*D$uwx/D$Fx
	subset(D,D$time>100000)
	return(D)
}

ReadCase<-function(PATTERN){
    list<-list.files(pattern=PATTERN)
    Data=NULL
    for (k in 1:length(list)){
        print(list[k])
        tmp=ReadRelperm(list[k])
	tmp$Case<-list[k]
        Data<-rbind(Data,tmp)
    }
    return(Data)
}


D<-ReadCase("benth_w")
require("ggplot2")

B<-read.csv("bentheimer/drain.csv",head=TRUE,sep=" ")
B$Media<-"Bentheimer"
M1<-read.csv("mineral_model_1/drain.csv",head=TRUE,sep=" ")
M1$Media<-"Mineral Sandstone #1"

p<-ggplot()+
	geom_line(data=B,aes(Sw,2/R,colour=Media))+
	geom_line(data=M1,aes(Sw,2/R,colour=Media))+
	theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
	theme_bw()

ggsave("morph-drain.png",p,height=4.0,width=6.0)