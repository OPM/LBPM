import numpy as np
import math
import matplotlib.pyplot as plt

#physical constant
k_B_const = 1.380649e-23  #[J/K]
N_A_const = 6.02214076e23 #[1/mol]
e_const = 1.602176634e-19 #[C]
epsilon0_const = 8.85418782e-12 #[C/V/m]

#other material property parameters
epsilonr_water = 80.4
T=310.15 #[K]

#input ion concentration
C_Na_in  = 15e-3  #[mol/m^3]
C_Na_out = 145e-3 #[mol/m^3]
C_K_in   = 150e-3 #[mol/m^3]
C_K_out  = 4e-3   #[mol/m^3]
C_Cl_in  = 10e-3  #[mol/m^3]
C_Cl_out = 110e-3 #[mol/m^3] 

#calculating Debye length
#For the definition of Debye lenght in electrolyte solution, see:
#DOI:10.1016/j.cnsns.2014.03.005
#Eq(42) in Yoshida etal., Coupled LB method for simulator electrokinetic flows
prefactor= math.sqrt(epsilonr_water*epsilon0_const*k_B_const*T/2.0/N_A_const/e_const**2)
debye_length_in  = prefactor*np.sqrt(np.array([1.0/C_Na_in,1.0/C_K_in,1.0/C_Cl_in]))
debye_length_out = prefactor*np.sqrt(np.array([1.0/C_Na_out,1.0/C_K_out,1.0/C_Cl_out]))
print("Debye length inside membrane in [m]")
print(debye_length_in)
print("Debye length outside membrane in [m]")
print(debye_length_out)


#setup domain
cube_length_z  = 192
cube_length_xy = 64
#set LBPM domain resoluiton
h=0.01 #[um]
print("Image resolution = %.6g [um] (= %.6g [m])"%(h,h*1e-6))

domain=2*np.ones((cube_length_z,cube_length_xy,cube_length_xy),dtype=np.int8)
zgrid,ygrid,xgrid=np.meshgrid(np.arange(cube_length_z),np.arange(cube_length_xy),np.arange(cube_length_xy),indexing='ij')
domain_centre=cube_length_xy/2
make_bubble = np.logical_and(zgrid>=cube_length_z/4,zgrid<=cube_length_z*0.75)
domain[make_bubble]=1

##save domain
file_name= "Pseudo3D_double_plane_membrane_z192_xy64_InsideLabel1_OutsideLabel2.raw"
domain.tofile(file_name)
print("save file: "+file_name)

#debug plot
#plt.figure(1)
#plt.pcolormesh(domain[:,int(domain_centre),:])
#plt.colorbar()
#plt.axis("equal")
#plt.show()

##generate initial ion concentration - 3D
#domain_Na = C_Na_out*np.ones_like(domain,dtype=np.float64)
#domain_Na[make_bubble] = C_Na_in
#domain_K = C_K_out*np.ones_like(domain,dtype=np.float64)
#domain_K[make_bubble] = C_K_in
#domain_Cl = C_Cl_out*np.ones_like(domain,dtype=np.float64)
#domain_Cl[make_bubble] = C_Cl_in
#
#domain_Na.tofile("Pseudo3D_plane_membrane_concentration_Na_z192_xy64.raw")
#domain_K.tofile("Pseudo3D_plane_membrane_concentration_K_z192_xy64.raw")
#domain_Cl.tofile("Pseudo3D_plane_membrane_concentration_Cl_z192_xy64.raw")

##debug plot
#plt.figure(2)
#plt.subplot(1,3,1)
#plt.title("Na concentration")
#plt.pcolormesh(domain_Na[:,int(bubble_centre),:])
#plt.colorbar()
#plt.axis("equal")
#plt.subplot(1,3,2)
#plt.title("K concentration")
#plt.pcolormesh(domain_K[:,int(bubble_centre),:])
#plt.colorbar()
#plt.axis("equal")
#plt.subplot(1,3,3)
#plt.title("Cl concentration")
#plt.pcolormesh(domain_Cl[:,int(bubble_centre),:])
#plt.colorbar()
#plt.axis("equal")
#plt.show()
