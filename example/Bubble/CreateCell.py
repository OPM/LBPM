import numpy as np
import matplotlib.pylab as plt

D=np.ones((40,40,40),dtype="uint8")

cx = 20
cy = 20
cz = 20

for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            if (dist < 15.5 ) :
                D[i,j,k] = 2

D.tofile("cell_40x40x40.raw")

                
C1=np.zeros((40,40,40),dtype="double")
C2=np.zeros((40,40,40),dtype="double")
C3=np.zeros((40,40,40),dtype="double")
C4=np.zeros((40,40,40),dtype="double")
C5=np.zeros((40,40,40),dtype="double")
C6=np.zeros((40,40,40),dtype="double")

for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            #outside the cell
            C1[i,j,k] = 4.0e-6    # K 
            C2[i,j,k] = 150.0e-6  # Na
            C3[i,j,k] = 116.0e-6  # Cl
            C4[i,j,k] = 29.0e-6   # HC03
            #C5[i,j,k] = 2.4e-6    # Ca
            dist = np.sqrt((i-cx)*(i-cx) + (j-cx)*(j-cx) + (k-cz)*(k-cz))
            # inside the cell
            if (dist < 15.5 ) :
                C1[i,j,k] = 145.0e-6
                C2[i,j,k] = 12.0e-6
                C3[i,j,k] = 4.0e-6  
                C4[i,j,k] = 12.0e-6  # 12 mmol / L
                #C5[i,j,k] = 0.10e-6  # 100 nmol / L


# add up the total charge to make sure it is zero
TotalCharge = 0
for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            TotalCharge += C1[i,j,k] + C2[i,j,k] - C3[i,j,k] - C4[i,j,k] 

TotalCharge /= (40*40*40)

print("Total charge " + str(TotalCharge))


for i in range(0, 40):
    for j in range (0, 40):
        for k in range (0,40):
            if TotalCharge < 0 :
                # need more cation
                C5[i,j,k] = abs(TotalCharge)
                C6[i,j,k] = 0.0
            else :
                # need more anion
                C5[i,j,k] = 0.0
                C6[i,j,k] = abs(TotalCharge)
        

C1.tofile("cell_concentration_K_40x40x40.raw")
C2.tofile("cell_concentration_Na_40x40x40.raw")
C3.tofile("cell_concentration_Cl_40x40x40.raw")
C4.tofile("cell_concentration_HCO3_40x40x40.raw")
C5.tofile("cell_concentration_cation_40x40x40.raw")
C6.tofile("cell_concentration_anion_40x40x40.raw")

