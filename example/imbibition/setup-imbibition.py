import os
import sys
import csv
import glob

name=sys.argv[1]
#numpts=int(sys.argv[2])
process="drain"
cwd=os.getcwd()

print("Initializing imbibition sequence")
print("Drainage should initialize imbibition!")
print("Base name is "+name)

# Parameters for color LBM
tau=0.8
alpha=0.005
beta=0.95
phisol=-1.0
saturation=0.0
Fx=0.0
Fy=0.0
Fz=0.0
Restart=1
pBC=3
din=1.001
dout=0.999
maxtime=2000005
interval=50000
tolerance=1e-5

viscosity=(tau-0.5)/3
ift=5.796*alpha

radius=[]
#porsize.quartiles provides quartiles for the pore size
with open("poresize.quartiles","r") as f:
    for line in f:
        reader=csv.reader(f,delimiter=' ')
        for row in reader:
            radius.append(float(row[1]))
            #sw.append(float(row[1]))

#compute the pressure difference
# Use Q1 as the initial pressure
din=2*ift/radius[0]
# Use Q4 (maximum pore size) to set final pressure differnece
dout=2*ift/(radius[3]+1.0)

# iterator for the filenames matching drain in the current directory
#Source=glob.iglob('*drain_*')

#dirname=str(name)+"_"+str(process)+"_"+str(pt)
#print("Creating " + dirname)
#if not os.path.exists(dirname):
#    os.mkdir(dirname)
#os.chdir(dirname)

ParamFile = open("Color.in.imbibition","w")
ParamFile.write("%f\n" % tau)
ParamFile.write("%f " % alpha)
ParamFile.write("%f " % beta)
ParamFile.write("%f\n" % phisol)
ParamFile.write("%f\n" % saturation)
ParamFile.write("%f " % Fx)
ParamFile.write("%f " % Fy)
ParamFile.write("%f\n" % Fz)
ParamFile.write("%i " % Restart)
ParamFile.write("%i " % pBC)
ParamFile.write("%f " % din)
ParamFile.write("%f\n" % dout)
ParamFile.write("%i " % maxtime)
ParamFile.write("%i " % interval)
ParamFile.write("%f\n" % tolerance)
ParamFile.close()
 
#os.chdir(c)
