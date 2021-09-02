import sys

infile = sys.argv[1]
outfile = infile.replace(".vtu", "comp.vtu", 1)

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
myinput = XMLUnstructuredGridReader(FileName=[infile])
myinput.PointArrayStatus = ['circulation', 'radius', 'velocity']

# save data
SaveData(outfile, proxy=myinput, DataMode='Binary', CompressorType='ZLib', CompressionLevel='9')

