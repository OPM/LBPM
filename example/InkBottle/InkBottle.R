require(ggplot2)

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

time=seq(1,100)
time_inlet<-time

M0_inlet<-CapVolume(r1,r1)+CylVolume(r1,time*0.03)
M1_inlet<-CapArea(r1,r1)+CylArea(r1,time*0.03)
M2_inlet<-2*CapArea(r1,r1)/r1+CylArea(r1,time*0.03)/r1
Jwn_inlet<-2/r1

# radius of the fluid interface while contact line is pinned
rt_pore<-r1+0.01*(Rb-r1)*time
# height of the spherical cap at the inlet side 
ht_pore<-(-rt_pore+sqrt(rt_pore^2+r1^2))
# depth of penetration into the pore
xt_pore<-2*rt_pore-ht_pore

time_pore<-time+100;
M0_pore<-CylVolume(r1,L1)+CapVolume(r1,xt_pore)
M1_pore<-CylArea(r1,L1)+CapArea(rt,xt_pore)
M2_pore<-2*CapArea(r1,r1)/r1+CylArea(r1,L1)/r1+2*CapVolume(r1,xt_pore)/rt_pore
Jwn_pore<-2/rt_pore

time_outlet<-time_pore+100;
Jwn_outlet<-2/r2
M0_outlet<-CylVolume(r1,L1)+BubVolume(r1,r2,Hb)+CapVolume(r2,r2)+CylVolume(r2,time*0.03)
M1_outlet<-CylArea(r1,L1)+CapArea(Rb,Hb)+CapArea(r2,r2)+CylArea(r2,time*0.03)
M2_outlet<-2*CapArea(r2,r2)/r2+CylArea(r2,time*0.03)/r2

