import os
import sys
import csv

name=sys.argv[1]
#numpts=int(sys.argv[2])
process="drain"
cwd=os.getcwd()

print("Initializing morphological drainage cases")
print("Base name is "+name)

# Parameters for color LBM
tau=0.7
alpha=0.005
beta=0.95
phisol=-1.0
saturation=0.0
Fx=0.0
Fy=0.0
Fz=0.0
Restart=0
pBC=1
din=1.001
dout=0.999
maxtime=100005
interval=50000
tolerance=1e-5

viscosity=(tau-0.5)/3
ift=5.796*alpha

radius=[]
sw=[]
with open("morphdrain.csv","r") as f:
    for line in f:
        reader=csv.reader(f,delimiter=' ')
        for row in reader:
            radius.append(float(row[0]))
            sw.append(float(row[1]))

numpts=len(radius)
print("Number of cases "+str(numpts))

for pt in range(0,numpts):
    #compute the pressure difference
    tag=pt+1
    dp=2*ift/radius[pt]
    din=1.0+0.5*dp
    dout=1.0-0.5*dp
    dirname=str(name)+"_"+str(process)+"_"+str(tag)
    print("Creating " + dirname)
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    ParamFile = open("Color.in","w")
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

    os.chdir(cwd)
