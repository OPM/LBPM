#!/bin/env python
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################


import numpy as np

def segmentGrid32bit(inputFileName, outputFileName):
    inputFile=inputFileName
    value=np.fromfile(inputFile,dtype=np.float32)
    value_new = np.copy(value)

    value_new[value < 0.1] = 1
    value_new[value > -0.01 and value < 0.01] = 0
    value_new[value > 0.1] = 2


    value_new.tofile(outputFileName)

def segmentGrid8bit(inputFileName, outputFileName): 
    inputFile=inputFileName
    value=np.fromfile(inputFile,dtype=np.uint8)
    value_new = np.copy(value)
    
    value_new[value == 0] = 0
    value_new[value == 127] = 1
    value_new[value == 255] = 0
        
    value_new.tofile(outputFileName)


def segmentGrid8bitRedPhase(path): 
    import os
    for inputFile in os.listdir(path):
        if os.path.isfile(inputFile) and inputFile.endswith('.raw') and inputFile.startswith('id_t'):
            value=np.fromfile(inputFile,dtype=np.uint8)
            value_new = np.copy(value)
    
            value_new[value == 0] = 0
            value_new[value == 1] = 255
            value_new[value >= 2] = 0
        
            value_new.tofile('Red_'+inputFile)

def segmentGrid8bitRockPhase(inputFileName): 
    inputFile=inputFileName
    value=np.fromfile(inputFile,dtype=np.uint8)
    value_new = np.copy(value)
    
    value_new[value == 0] = 255
    value_new[value == 1] = 0
    value_new[value >= 2] = 0
    
    value_new.tofile('Rock_'+inputFile)



def unSegmentGrid8bit(inputFileName, outputFileName): 
    inputFile=inputFileName
    value=np.fromfile(inputFile,dtype=np.uint8)
    value_new = np.copy(value)
    
    value_new[value == 0] = 0
    value_new[value == 1] = 127
    value_new[value == 2] = 255
        
    value_new.tofile(outputFileName)



def segmentGrid32to8bit(inputFileName, outputFileName): 
## Something wrong with this...
    inputFile=inputFileName
    value=np.fromfile(inputFile,dtype=np.float32)
    value_new = np.copy(value)
    value_new.astype(np.uint8)
    
    
    value_new[value < -0.2] = 1
    value_new[value > -0.1, value < 0.1] = 0
    value_new[value > 0.1] = 2
        
    value_new.tofile(outputFileName)


def addVoidLayers(inputFileName, outputFileName, nx, ny, nz, layers):
    inputFile=open(inputFileName, 'r')
    value=np.fromfile(inputFile,dtype=np.uint8)
    inputFile.close()
    value_3D = np.reshape(value, (nz,ny,nx))
    value_3D = np.pad(value_3D, ((0,0), (0,0), (layers,layers)),
                          'constant', constant_values=((1,1)))
                          
    value_3D = np.pad(value_3D, ((1,1), (1,1), (0,0)),
                          'constant', constant_values=((0,0)))

    value_3D.tofile(outputFileName)

# from and including -- to and including
def crop(inputFileName, outputFileName, nx, ny, nz, x1, x2, y1, y2, z1, z2):
    inputFile=open(inputFileName, 'r')
    value=np.fromfile(inputFile,dtype=np.uint8)
    inputFile.close()
    value_3D = np.reshape(value, (nz,ny,nx))
    value_3D_crop = value_3D[z1-1:z2, y1-1:y2, x1-1:x2]
    value_3D_crop.tofile(outputFileName)


def circleTube(outputFileName,nx,ny,nz): 

    tube_3D_new = np.zeros((nz, ny, nx), dtype=np.uint8)
    centerZ = nz/2
    centerY = ny/2
    radius = ny/2
    
    
    # outer loop to handle number of rows 
    for i in range(0, nx): 
      
        # inner loop to handle number spaces 
        # values changing acc. to requirement 
        for j in range(0, ny):
            for k in range(0, nz): 
                if (((j-centerY)**2 + (k-centerZ)**2) < radius**2):
                    tube_3D_new[k,j,i] = tube_3D_new[k,j,i] + 1

    tube_3D_new.tofile(outputFileName)

def squareTube(outputFileName,nx,ny,nz): 

    tube_3D_new = np.zeros((nz, ny, nx), dtype=np.uint8) + 1
    tube_3D_new = np.pad(tube_3D_new, ((1,1), (1,1), (0,0)),
                          'constant', constant_values=((0,0)))

    
    tube_3D_new.tofile(outputFileName)

def makeVtk(headerVTK, path) :
    import os
    for f in os.listdir(path):
        if os.path.isfile(f) and f.endswith('.raw') and f.startswith('Red_id_t'):
            infile = open(f, 'rb')
            header = open(headerVTK, 'r')
            outFileName = f.split('.')[0]+'.vtk'
            outFile = open( outFileName, 'w')
            outFile.write(header.read())
            outFile.write(infile.read())
            outFile.close()
            infile.close()
            header.close()
                
            
                

    
  
# Driver Code
#segmentGrid8bitRedPhase('path')
#segmentGrid8bitRockPhase('id_t100008.raw')


#makeVtk('headerVTK.txt', 'path')

#circleTube('Tube.raw',100,101,101)

#squareTube('SquareTube.raw',100,100,100)

#crop('BentheimerSegmWW8bit2.raw', 'Outout2.raw', 1432, 461, 457, 100, 699, 5, 404, 5, 404) 
#addVoidLayers('Outout2.raw', 'BentheimerSegmWW8bit2CropCoated600x400x400.raw', 600, 400, 400, 0) 


