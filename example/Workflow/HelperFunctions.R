require("ggplot2")

ReadSubphase<-function(PATH){
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

ReadTimelog<-function(PATH){
	FILE=paste0(PATH,"/timelog.csv")
	D<-read.csv(file=FILE,head=TRUE,sep=" ")
	D$time<-seq(1,nrow(D))
	return(D)
}

ReadRelperm<-function(PATH){
	FILE=paste0(PATH,"/relperm.csv")
	D<-read.csv(file=FILE,head=TRUE,sep=" ")
	D$Case<-PATH

	p<-ggplot(D)+
	  geom_line(aes(sat.water,eff.perm.oil,color="oil"))+
	  geom_line(aes(sat.water,eff.perm.water,color="water"))+
	  xlab("Water saturation")+ylab("Effective permeability")+
	  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
	  theme_bw()

	FILE=paste0(PATH,"-relperm.png")
 	ggsave(FILE,p,height=4.0,width=6.0)

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


ReadMorphDrain<-function(PATH){
	FILE=paste0(PATH,"/morphdrain.csv")
	B<-read.csv(FILE,head=TRUE,sep=" ")
	B$Media<-PATH

	p<-ggplot()+
	  geom_line(data=B,aes(Sw,2/R,colour=Media))+
	  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
	  theme_bw()

	FILE=paste0(PATH,"-morphdrain.png")
 	ggsave(FILE,p,height=4.0,width=6.0)
	return(B)
}

