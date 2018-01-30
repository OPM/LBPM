require("ggplot2")

GGTHEME<-theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

PI=3.14159

CapVolume<-function(a,h){
	volume<-PI*h*(3*a^2+h^2)/6
	return(volume)
}

CapArea<-function(R,h){
	area<-2*PI*R*h
	return(area)
}

BulbVolume<-function(a1,a2,h){
	volume<-PI*h*(3*a1^2+3*a2^2+h^2)/6
	return(volume)
}

CylArea<-function(R,h){
	area<-2*PI*R*h
	return(area)
}

CylVolume<-function(R,h){
	volume<-PI*R^2*h
	return(volume)
}

# Geometry
L1 = 3.0
L2 = 3.0
r1=1.0
r2=0.8
#Radius of the bulb
Rb=1.2
# Length of the bulb
Hb=sqrt(Rb^2-r1^2)+sqrt(Rb^2-r2^2)

time=seq(0,100)
time_inlet<-time

M0_inlet<-CapVolume(r1,r1)+CylVolume(r1,time*0.01*L1)
M1_inlet<-CapArea(r1,r1)+CylArea(r1,time*0.01*L1)
M2_inlet<-2*CapArea(r1,r1)/r1+CylArea(r1,time*0.03)/r1
Jwn_inlet<-2/r1

# radius of the fluid interface while contact line is pinned
rt_pore<-r1+0.01*(Rb-r1)*time
# height of the spherical cap at the inlet side 
ht_pore<-(rt_pore-sqrt(rt_pore^2-r1^2))
# depth of penetration into the pore
xt_pore<-2*rt_pore-ht_pore

time_pore<-time+100
M0_pore<-CylVolume(r1,L1)+pi*xt_pore*(3*r1^2+xt_pore^2)/6
M1_pore<-CylArea(r1,L1)+2*pi*rt_pore*xt_pore
M2_pore<-CylArea(r1,L1)/r1 + 4*pi*xt_pore
Jwn_pore<-2/rt_pore

rt_pin<-Rb+(r2-Rb)*time*0.01
ht_pin<-rt_pin-sqrt(rt_pin^2-r2^2)
M0_pin<-CylVolume(r1,L1)+BulbVolume(r1,r2,Hb)+pi*ht_pin*(3*r2^2+ht_pin^2)/6
M1_pin<-CylArea(r1,L1) + 2*pi*Rb*Hb + 2*pi*rt_pin*ht_pin
M2_pin<-CylArea(r1,L1)/r1 + 4*pi*Hb + 4*pi*ht_pin
Jwn_pin<-2/rt_pin

time_outlet<-time_pore+100;
Jwn_outlet<-2/r2
M0_outlet<-CylVolume(r1,L1)+BulbVolume(r1,r2,Hb)+CapVolume(r2,r2)+CylVolume(r2,time*0.01*L2)
M1_outlet<-CylArea(r1,L1)+CapArea(Rb,Hb)+CapArea(r2,r2)+CylArea(r2,time*0.01*L2)
M2_outlet<-CylArea(r1,L1)/r1+2*CapArea(Rb,Hb)/Rb+2*CapArea(r2,r2)/r2+CylArea(r2,time*0.01*L2)/r2

M0<-as.vector(rbind(M0_inlet,M0_pore,M0_pin,M0_outlet))
M1<-as.vector(rbind(M1_inlet,M1_pore,M1_pin,M1_outlet))
M2<-as.vector(rbind(M2_inlet,M2_pore,M2_pin,M2_outlet))
Jwn<-as.vector(rbind(Jwn_inlet,Jwn_pore,Jwn_pin,Jwn_outlet))

Vol=max(M0)
sw<-1.0-M0/Vol

p1<-ggplot()+geom_line(aes(time,M0_inlet,color="inlet"))+
	geom_line(aes(time+100,M0_pore,color="pore"))+
	geom_line(aes(time+200,M0_pin,color="pin"))+
	geom_line(aes(time+300,M0_outlet,color="outlet"))

p2<-ggplot()+geom_line(aes(time,M1_inlet,color="inlet"))+
	geom_line(aes(time+100,M1_pore,color="pore"))+
	geom_line(aes(time+200,M1_pin,color="pin"))+
	geom_line(aes(time+300,M1_outlet,color="outlet"))

p3<-ggplot()+geom_line(aes(time,M2_inlet,color="inlet"))+
	geom_line(aes(time+100,M2_pore,color="pore"))+
	geom_line(aes(time+200,M2_pin,color="pin"))+
	geom_line(aes(time+300,M2_outlet,color="outlet"))


p4<-ggplot()+geom_line(aes(sw,Jwn))+GGTHEME+
	xlab(expression(s^w))+ylab(expression(p^c))

p5<-ggplot()+geom_line(aes(M0,M1))+GGTHEME+
	xlab(expression(M[0]))+ylab(expression(M[1]))

p6<-ggplot()+geom_line(aes(M0,M2))+GGTHEME+
	xlab(expression(M[0]))+ylab(expression(M[2]))

ggsave("InkBottle-Jwn-sw.png",p4,width=5.0,height=4.0)
ggsave("InkBottle-M0-M1.png",p5,width=5.0,height=4.0)
ggsave("InkBottle-M0-M2.png",p6,width=5.0,height=4.0)
