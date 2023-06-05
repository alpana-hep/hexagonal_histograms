#!/usr/bin/env python2
import numpy as np
import math
import ROOT
ROOT.gROOT.SetBatch(True)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-n', help="number of cells", type=int, default=9999)
args = parser.parse_args()

UNTIL_THIS_NUMBER = args.n 

# load shared library
ROOT.gInterpreter.ProcessLine('#include "include/HGCalCell.h"')
ROOT.gSystem.Load("./build/libHGCalCell.so")

# input parameters
waferSize = 60
nFine, nCoarse = 0, 10 #222
typeFine, typeCoarse = 0, 1
placementIndex = 0

# coordinate of a hexagonal cell
s3 = math.sqrt(3)
nCorner = 7 # one corner is repeated
hexagon_base = {
#	'x': [2, 1, -1, -2, -1, 1, 2],
#	'y': [0, -1*s3, -1*s3, 0, s3, s3, 0]
	'x': [0, 1*s3, 1*s3, 0, -1*s3, -1*s3, 0],
 	'y': [2, 1, -1, -2, -1, 1, 2]
}

# translation & rotation
global_correction_x, global_correction_y = 2*s3, -5 
global_theta = 5.*math.pi/6. # 150 degree
cos_theta, sin_theta = math.cos(global_theta), math.sin(global_theta)

# load geometry data
fin = open("./data/WaferCellMapTrg.txt", 'r')
contents = fin.readlines()[:223]
fin.close()

# loop over all the cells
counter = 0 # how many cells are drawn
collections = {} # key, value = SiCell, graph
cell_helper = ROOT.HGCalCell(waferSize, nFine, nCoarse)
for i, line in enumerate(contents):
	if i==0: continue # omit heading
	if counter==UNTIL_THIS_NUMBER : break # manually control

	density, _, roc, half, seq, rocpin, sicell, _, _, iu, iv = tuple(line.split()[:11])
	if(iu=="-1" and iv=="-1"): continue # ignore (-1,-1)
	if("CALIB" in rocpin): continue # ignore calibration bins
	if(density == "HD"): break # keep only first set of LD

	coor = cell_helper.cellUV2XY1(int(iu), int(iv), 0, typeCoarse)
	x, y = coor[0], coor[1]

	# evaluate (r, phi) and apply rotation
	r = math.sqrt(pow(x,2)+pow(y,2))
	cos_phi, sin_phi = x/r, y/r
	xprime = r*(cos_phi*cos_theta + sin_phi*sin_theta) + global_correction_x
	yprime = r*(sin_phi*cos_theta - cos_phi*sin_theta) + global_correction_y
	x, y = xprime, yprime

	# create a hexagon
	hexagon = {}
	hexagon['x'] = [ element + x for element in hexagon_base['x'] ]
	hexagon['y'] = [ element + y for element in hexagon_base['y'] ]

	graph = ROOT.TGraph(nCorner, np.array(hexagon['x']), np.array(hexagon['y']))
	graph.SetTitle("")
	graph.GetXaxis().SetTitle("x (arb. unit)")
	graph.GetYaxis().SetTitle("y (arb. unit)")
	
	graph.SetMaximum(200)
	graph.SetMinimum(-200)
	graph.GetXaxis().SetLimits(-200, 200)
	
	graph.SetName("hex_%d" % int(sicell))
	collections[int(sicell)] = graph
	counter+=1

	print "counter=%d, i=%d, (iu,iv) = (%d,%d), (x,y) = (%.2f, %.2f)" % (counter, i, int(iu), int(iv), x, y)

# store graphs in order of sicell
fout = ROOT.TFile("./data/hexagons.root", "RECREATE")
for sicell, graph in collections.items():
	#print "sicell = %d" % sicell
	graph.Write()
fout.Close()

# execute root macro for TH2Poly
import subprocess
subprocess.call("root -l -b -q th2poly.C", shell=True)
