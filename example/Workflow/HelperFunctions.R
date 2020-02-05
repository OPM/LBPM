require("ggplot2")

GG_THEME=theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ReadSubphase<-function(PATH){
  FILE=paste0(PATH,"/subphase.csv")
  S<-read.csv(FILE,head=TRUE,sep=" ")
  S$Vw<-S$Vwc+S$Vwd
  S$Vn<-S$Vnc+S$Vnd
  S$Aw<-S$Awc+S$Awd
  S$An<-S$Anc+S$And
  S$Hw<-S$Hwc+S$Hwd
  S$Hn<-S$Hnc+S$Hnd
  S$Xw<-S$Xwc+S$Xwd
  S$Xn<-S$Xnc+S$Xnd

  S$Sw<-S$Vw/(S$Vn+S$Vw)
  S$pw<-(S$pwc*S$Vwc+S$pwd*S$Vwd) / (S$Vwc+S$Vwd)
  S$pn<-(S$pnc*S$Vnc+S$pnd*S$Vnd) / (S$Vnc+S$Vnd)

  S$Qwx<-S$Vw*(S$Pwc_x+S$Pwd_x)/(S$Mwc+S$Mwd)
  S$Qnx<-S$Vn*(S$Pnc_x+S$Pnd_x)/(S$Mnc+S$Mnd)
  S$Qwy<-S$Vw*(S$Pwc_y+S$Pwd_y)/(S$Mwc+S$Mwd)
  S$Qny<-S$Vn*(S$Pnc_y+S$Pnd_y)/(S$Mnc+S$Mnd)
  S$Qwz<-S$Vw*(S$Pwc_z+S$Pwd_z)/(S$Mwc+S$Mwd)
  S$Qnz<-S$Vn*(S$Pnc_z+S$Pnd_z)/(S$Mnc+S$Mnd)

  S$Krn<-S$nun*S$Qnz/S$Fz
  S$Krw<-S$nuw*S$Qwz/S$Fz
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

