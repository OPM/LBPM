from skimage import io
import glob
import numpy as np

print("Concatenating tiff files into 3D volume")

IMAGES=sorted(glob.glob("*.tif"))
data=np.array([])

mydata=[]
for im in IMAGES:
    print(im)
    image=io.imread(im)
    mydata.append(image)
#    data=np.append(data,image,2)

data=np.array(mydata)
data.tofile("data.raw")


