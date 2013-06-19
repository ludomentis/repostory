"""Heatmap generatro for Google Earth"""
import numpy as np
import math
import sys
import Image
import random

from pykml.factory import KML_ElementMaker as KML
from lxml import etree
from scipy.interpolate import Rbf
from upoints import utils

if len(sys.argv) != 6:
    print "Invalid arguments number! Usage:"
    print sys.argv[0]+" [filename] [width] [height] [pixel]"
    print "\t[filename]\t\tname of datafile"
    print "\t[width]\t\twidth of single square of aggregation (meters)"
    print "\t[height]\t\theight of single square of aggregation (meters)"
    print "\t[pixel] approximate image width (PIXELS)"
    sys.exit()

print "Reading " + sys.argv[1] + "...",
sys.stdout.flush()
DATA = np.genfromtxt(sys.argv[1])

print "done! " + str(len(DATA)) + " rows read!"

#caluclating size of aggregation in terms of lat and lon
RADIUS = 1000*utils.calc_radius(float(DATA[0][0])) #earth radius in meters
SS_W = float(sys.argv[2]) # single square width (west-east)
SS_H = float(sys.argv[3]) # single square height (north-south)
HSIZE = float(sys.argv[4])
TRANS_THRESH = float(sys.argv[5])

#fraction on coordinates used for aggregation
#i.e. aggregates amongst 1/FRAC_LAT degrees of latitude
FRAC_LAT = 1./math.degrees(SS_H/RADIUS)
#1600.0 is a good value
FRAC_LON = 1./math.degrees(SS_W/(RADIUS*math.cos(math.radians(DATA[0][0]))))
#800.0 is a good value

#the grid dictionary
MMAPP = {}
# key (lat,lon)
# value (cnt, sum, sumq)

print "Aggregating data...",
sys.stdout.flush()
for row in DATA:
    k = (round(float(row[0])*FRAC_LAT), round(float(row[1])*FRAC_LON))
    if k in MMAPP:
        MMAPP[k][0] += 1.
        MMAPP[k][1] += float(row[2])
        MMAPP[k][2] += float(row[2])*float(row[2])
    else:
        MMAPP[k] = [1, float(row[2]), float(row[2])*float(row[2])]
print "done!"

#max and min used later for palette
MA_V = max([x[1]/x[0] for x in MMAPP.values()])
MI_V = min([x[1]/x[0] for x in MMAPP.values()])

print "BC values"
print "    max: " + str(MA_V)
print "    min: " + str(MI_V)

#max and min of coverage used later for opacity
MA_C = max([math.log(x[0]) for x in MMAPP.values()])
MI_C = min([math.log(x[0]) for x in MMAPP.values()])
print "Coverage (log #measure)"
print "    max: " + str(MA_C)
print "    min: " + str(MI_C)

LATMIN = min([x[0] for x in MMAPP.keys()])
LATMAX = max([x[0] for x in MMAPP.keys()])
LONMIN = min([x[1] for x in MMAPP.keys()])
LONMAX = max([x[1] for x in MMAPP.keys()])
print "Area"
print "    lat min: " + str(LATMIN)
print "    lat max: " + str(LATMAX)
print "    lon min: " + str(LONMIN)
print "    lon max: " + str(LONMAX)

#depending on the input approx size and on the real
#size of the aggregation grid, calulate the size
#in pixel of each aggregated dot and thus
#the image size
DOTSIZE = int(HSIZE/(LONMAX-LONMIN + 1))
SQUARES = (LATMAX - LATMIN + 1)*(LONMAX-LONMIN + 1)
HSIZE = DOTSIZE*(LONMAX - LONMIN + 1)
VSIZE = DOTSIZE*(LATMAX - LATMIN + 1)

# setup scattered data value
X = [(a[1]-LONMIN + .5)*DOTSIZE for a in MMAPP.keys()]
Y = [(a[0]-LATMIN + .5)*DOTSIZE for a in MMAPP.keys()]
Z = [a[1]/a[0] for a in MMAPP.values()]

# setup scattered data coverage
X_C = []
Y_C = []
Z_C = []
INX = 0.5
INY = 0.5

for llat in range(int(LATMIN), int(LATMAX+1)):
    for llon in range(int(LONMIN), int(LONMAX+1)):
        X_C_CAND = INX*DOTSIZE
        Y_C_CAND = INY*DOTSIZE
        if (llat, llon) in MMAPP:
            X_C.append(X_C_CAND)
            Y_C.append(Y_C_CAND)
            Z_C.append(math.log(MMAPP[(llat, llon)][0]))
        elif random.random() < TRANS_THRESH:
            #the number of dot not covered can be uselessy huge
            #only a random subset will be taken
            X_C.append(X_C_CAND)
            Y_C.append(Y_C_CAND)
            Z_C.append(0)
        INX += 1
    INY += 1 
    INX = 0.5

print "Map"
print "    all SQUARES = ", int((LONMAX - LONMIN + 1)*(LATMAX - LATMIN + 1))
print "    rgb SQUARES = ", len(MMAPP)
print "    alpha SQUARES = ", len(Z_C)
print "    pic size (px) = ", int(HSIZE), "x", int(VSIZE), "=", int(HSIZE*VSIZE)

# Interpolation function for rgb and alpha
print "Calculating RBF (rgb)...",
RBF = Rbf(X, Y, Z, epsilon=DOTSIZE)
sys.stdout.flush()
print " done!"

print "Calculated RBF (alpha)...",
RBF_A = Rbf(X_C, Y_C, Z_C, epsilon=DOTSIZE)
sys.stdout.flush()
print " done!"

#image creation
BGR = Image.new("RGBA", (int(HSIZE), int(VSIZE)), (255, 255, 255, 0))
PIXELS = BGR.load()
print "Image initialized!"

#color-opacity function
def qcol((val, mini, maxi), (val_a, mini_a, maxi_a)):
    """function for rgba calculation"""
    rho = (val - mini)/(maxi-mini)
    if rho < 0.:
        rho = 0.
    elif rho > 1.:
        rho = 1.
    red = int(255*math.pow(rho*2, .4)) if rho < .5 else 255
    green = int(255-255*math.pow((rho-.5)*2, .4)) if rho > .5 else 255
    blue = 0
    alpha = int(255*(val_a-mini_a)/(maxi_a-mini_a))
    return (red, green, blue, alpha)

#image writing
for i in range(BGR.size[0]):    # for every pixel:
    sys.stdout.write("Image writing: %d%%\r" % (100*i/BGR.size[0]))
    sys.stdout.flush()
    for j in range(BGR.size[1]):
        RBF_v = RBF(i, BGR.size[1]-j-1)
        alpha_v = RBF_A(i, BGR.size[1]-j-1)
        PIXELS[i, j] = qcol((RBF_v, MI_V, MA_V), (alpha_v, MI_C, MA_C))
BGR.save("bgr.png")
print "Image saved!"

#kml creation
MYKML = KML.kml()
DOC = KML.Document()
MYGO = KML.GroundOverlay(
  KML.color('ffffffff'),
  KML.Icon(KML.href("bgr.png")),
  KML.LatLonBox(
    KML.north(str(LATMAX/FRAC_LAT)),
    KML.south(str(LATMIN/FRAC_LAT)),
    KML.east(str(LONMAX/FRAC_LON)),
    KML.west(str(LONMIN/FRAC_LON)),
    KML.rotation('0')
  )
)
DOC.append(MYGO)
MYKML.append(DOC)

#kml writing
OUTFILE = file('test.kml', 'w')
OUTFILE.write(etree.tostring(MYKML, pretty_print=True))
