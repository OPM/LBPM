import os
import sys

name=sys.argv[1]
numpts=int(sys.argv[2])
process="drain"
cwd=os.getcwd()

print(name)
print(numpts)

tau=0.7
alpha=0.005
beta=0.95
for pt in range(0,numpts):
    dirname=str(name)+"_"+str(process)+"_"+str(pt)
   # os.path.join(cwd,name,process,"pt"),
    print(cwd,dirname)
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    ParamFile = open("Color.in","w")
    ParamFile.write("%f\n" % tau)
    ParamFile.write("%f " % alpha)
    ParamFile.write("%f\n" % beta)
    ParamFile.close()
    os.chdir(cwd)
