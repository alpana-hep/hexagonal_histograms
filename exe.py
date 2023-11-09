#!/usr/bin/env python3
import numpy as np
import math
import ROOT
ROOT.gROOT.SetBatch(True)
print ('argument')
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-n', help="number of cells", type=int, default=9999)
args = parser.parse_args()
print ('after argument')
UNTIL_THIS_NUMBER = args.n 
print ('loading shared libraries ')

# load shared library
ROOT.gInterpreter.ProcessLine('#include "include/HGCalCell.h"')
ROOT.gSystem.Load("./build/libHGCalCell.so")
print ('after loading shared libraries')
# input parameters
arbUnit_to_cm = 17./24.
waferSize = 60 * arbUnit_to_cm
nFine, nCoarse = 0, 10 #222
typeFine, typeCoarse = 0, 1
placementIndex = 0

import toolbox.geometry as tg
# define type id
# coordinate of a hexagonal cell
s3 = tg.s3
base = tg.base

# LD wafer pentagons
LD_special_polygonal_cells = tg.LD_special_polygonal_cells
LD_special_polygonal_cells_all = tg.LD_special_polygonal_cells_all

# translation & rotation
global_correction_x, global_correction_y = 2*s3, -5 
global_theta = 5.*math.pi/6. # 150 degree
cos_global_theta, sin_global_theta = math.cos(global_theta), math.sin(global_theta)

# for auxiliary boundary lines on the wafer map
import json
dict_my_coordinate_data = {} # key = sicell, value = dict_polygon_coordinates

def get_polygon(sicell, type_polygon, nCorner, x, y, isCM=False, isNC=False):
	polygon_base = base[type_polygon]
	polygon = {}

	resize_factor = 1.0
	if isCM:
		resize_factor = 0.6
	elif isNC:
		resize_factor = 0.4

	polygon['x'] = [ element*arbUnit_to_cm*resize_factor + x for element in polygon_base['x'] ]
	polygon['y'] = [ element*arbUnit_to_cm*resize_factor + y for element in polygon_base['y'] ]

	#if not isCM:
	#	polygon['x'] = [ element*factor + x for element in polygon_base['x'] ]
	#	polygon['y'] = [ element*factor + y for element in polygon_base['y'] ]

	#else: # apply rotation to CM cell coordinates
	#	idxCM = sicell - 198 - 1
	#	theta = 2*math.pi/3. * (idxCM//4)
	#	cos_theta = math.cos(theta)
	#	sin_theta = math.sin(theta)

	#	r = 2.
	#	polygon['x'] = [ element*cos_theta + math.sqrt(pow(r,2)-pow(element,2))*sin_theta + x for element in polygon_base['x'] ]
	#	polygon['y'] = [ element*cos_theta - math.sqrt(pow(r,2)-pow(element,2))*sin_theta + y for element in polygon_base['y'] ]
	#	print idxCM, idxCM//4, nCorner, len(polygon['x'])

	dict_my_coordinate_data[sicell] = polygon

	graph = ROOT.TGraph(nCorner+1, np.array(polygon['x']), np.array(polygon['y']))
	graph.SetTitle("")
	graph.GetXaxis().SetTitle("x (arb. unit)")
	graph.GetYaxis().SetTitle("y (arb. unit)")
	
	graph.SetMaximum(200)
	graph.SetMinimum(-200)
	graph.GetXaxis().SetLimits(-200, 200)
	return graph

# utility
counter = 0 # how many cells are drawn
collections = {} # key, value = SiCell, graph
cell_helper = ROOT.HGCalCell(waferSize, nFine, nCoarse)

# load geometry data
fin = open("./data/WaferCellMapTrg.txt", 'r')
contents = fin.readlines()[:223]
fin.close()

# loop over all the cells
for i, line in enumerate(contents):
	if i==0: continue # omit heading
	if counter==UNTIL_THIS_NUMBER : break # manually control

	density, _, roc, halfroc, seq, rocpin, sicell, _, _, iu, iv, t = tuple([str(ele) if "LD" in ele or "CALIB" in ele else int(ele) for ele in line.strip().split()])
	if(iu==-1 and iv==-1): continue # ignore (-1,-1)
	if(density == "HD"): break # keep only first set of LD

	globalId = 78*roc + 39*halfroc + seq

	# # print globalId vs HGCROC pin
	# print("{{{0},{1}}},").format(globalId, rocpin)

	# # print globalId vs padId
	# print("{{{0},{1}}},").format(globalId, sicell)
	
	coor = cell_helper.cellUV2XY1(int(iu), int(iv), 0, typeCoarse)
	x, y = coor[0], coor[1]

	# evaluate (r, phi) and apply rotation
	r = math.sqrt(pow(x,2)+pow(y,2))
	cos_phi, sin_phi = x/r, y/r
	xprime = r*(cos_phi*cos_global_theta + sin_phi*sin_global_theta) + global_correction_x
	yprime = r*(sin_phi*cos_global_theta - cos_phi*sin_global_theta) + global_correction_y
	x, y = xprime, yprime

	# create a hexagon
	type_polygon, nCorner = tg.type_hexagon, 6
	
	if sicell in LD_special_polygonal_cells[tg.type_hexagon_corner1]:
		type_polygon, nCorner = tg.type_hexagon_corner1, 6
	elif sicell in LD_special_polygonal_cells[tg.type_hexagon_corner2]:
		type_polygon, nCorner = tg.type_hexagon_corner2, 6
	elif sicell in LD_special_polygonal_cells[tg.type_hexagon_corner3]:
		type_polygon, nCorner = tg.type_hexagon_corner3, 6
	elif sicell in LD_special_polygonal_cells[tg.type_hexagon_corner4]:
		type_polygon, nCorner = tg.type_hexagon_corner4, 6
	elif sicell in LD_special_polygonal_cells[tg.type_hexagon_corner5]:
		type_polygon, nCorner = tg.type_hexagon_corner5, 6
	elif sicell in LD_special_polygonal_cells[tg.type_hexagon_corner6]:
		type_polygon, nCorner = tg.type_hexagon_corner6, 6

	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_side1]:
		type_polygon, nCorner = tg.type_pentagon_side1, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_side2]:
		type_polygon, nCorner = tg.type_pentagon_side2, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_side3]:
		type_polygon, nCorner = tg.type_pentagon_side3, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_side4]:
		type_polygon, nCorner = tg.type_pentagon_side4, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_side5]:
		type_polygon, nCorner = tg.type_pentagon_side5, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_side6]:
		type_polygon, nCorner = tg.type_pentagon_side6, 5

	elif sicell in tg.hollow_cells:
		type_polygon, nCorner = tg.type_hollow, 14 
	elif(isinstance(rocpin, str)): # "CALIB"
		type_polygon, nCorner = tg.type_hexagon_small, 6 
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_corner1]:
		type_polygon, nCorner = tg.type_pentagon_corner1, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_corner2]:
		type_polygon, nCorner = tg.type_pentagon_corner2, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_corner3]:
		type_polygon, nCorner = tg.type_pentagon_corner3, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_corner4]:
		type_polygon, nCorner = tg.type_pentagon_corner4, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_corner5]:
		type_polygon, nCorner = tg.type_pentagon_corner5, 5
	elif sicell in LD_special_polygonal_cells[tg.type_pentagon_corner6]:
		type_polygon, nCorner = tg.type_pentagon_corner6, 5

	graph = get_polygon(sicell, type_polygon, nCorner, x, y)	
	graph.SetName("hex_%d" % sicell)
	#collections[sicell] = graph
	collections[globalId] = graph
	counter+=1

	#print "counter=%d, i=%d, (iu,iv) = (%d,%d), (x,y) = (%.2f, %.2f)" % (counter, i, int(iu), int(iv), x, y)

# Add additional cells for CM channels
CMIds = [37, 38, 76, 77, 115, 116, 154, 155, 193, 194, 232, 233]
for idxCM in range(12):
	if idxCM%2==0: # CM0
		type_polygon, nCorner = tg.type_regular_pentagon, 5
	else: # CM1
		type_polygon, nCorner = tg.type_square, 4

	# assign coordinates
	correction_fine_tune_y_coordinate = 2. if idxCM//4 == 0 else 0.
	x = tg.Coordinates_CM_channels['x'][idxCM%4]*arbUnit_to_cm
	y = tg.Coordinates_CM_channels['y'][idxCM%4]*arbUnit_to_cm - correction_fine_tune_y_coordinate*arbUnit_to_cm
	theta = 2*math.pi/3. * (idxCM//4) - math.pi/3.
	cos_theta = math.cos(theta)
	sin_theta = math.sin(theta)

	# evaluate (r, phi) and apply rotation
	r = math.sqrt(pow(x,2)+pow(y,2))
	cos_phi, sin_phi = x/r, y/r
	xprime = r*(cos_phi*cos_theta + sin_phi*sin_theta)
	yprime = r*(sin_phi*cos_theta - cos_phi*sin_theta)
	x, y = xprime, yprime

	sicell = 198+1+idxCM # artificial sicell for CM channels
	globalId = CMIds[idxCM]

	graph = get_polygon(sicell, type_polygon, nCorner, x, y, True)	
	graph.SetName("hex_cm_%d" % idxCM)
	#collections[sicell] = graph
	collections[globalId] = graph
	counter+=1

# Add additional cells for NC channels
NonConnIds = [8, 17, 19, 28, 47, 56, 58, 67, 86, 95, 97, 106, 125, 134, 136, 145, 164, 173, 175, 184, 203, 212, 214, 223]
for idxNC in range(24):
	#type_polygon, nCorner = tg.type_triangle, 3
	type_polygon, nCorner = tg.type_circle, 12 
	x = tg.Coordinates_NC_channels['x'][idxNC%8]*arbUnit_to_cm
	y = tg.Coordinates_NC_channels['y'][idxNC%8]*arbUnit_to_cm
	theta = 2*math.pi/3. * (idxNC//8) - math.pi/3.
	cos_theta = math.cos(theta)
	sin_theta = math.sin(theta)

	# evaluate (r, phi) and apply rotation
	r = math.sqrt(pow(x,2)+pow(y,2))
	cos_phi, sin_phi = x/r, y/r
	xprime = r*(cos_phi*cos_theta + sin_phi*sin_theta)
	yprime = r*(sin_phi*cos_theta - cos_phi*sin_theta)
	x, y = xprime, yprime

	sicell = 198+1+12+idxNC # artificial sicell for NC channels
	globalId = NonConnIds[idxNC]

	graph = get_polygon(sicell, type_polygon, nCorner, x, y, False, True)	
	graph.SetName("hex_nc_%d" % idxNC)
	collections[globalId] = graph
	counter+=1

# store graphs in order of key (sicell or globalId)
fout = ROOT.TFile("./data/hexagons1.root", "RECREATE")
for key, graph in collections.items():
	graph.Write()
	#if key>=78: break
fout.Close()

# store coordinates
with open("data/output_my_coordinate_data.json", 'w') as f:
	json.dump(dict_my_coordinate_data, f, indent=4)

#--------------------------------------------------
# external execute
#--------------------------------------------------
import subprocess

def exe(command):
	print ("\n>>> executing command, ", command)
	subprocess.call(command, shell=True)

#execute python script for coordinate queries
exe("./toolbox/coordinate_loader.py")

## chi2 

                                 
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_chi2\", \"Module_0_run_1691487584_chi2.pdf\", 26, \"false\",\"0\",\"1691487584\",6)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_chi2\", \"Module_1_run_1691487584_chi2.pdf\", 26, \"false\",\"1\",\"1691487584\",6)'")

                                                 
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_chi2\", \"Module_0_run_1691518350_chi2.pdf\", 26, \"false\",\"0\",\"1691518350\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\",  \"Module_1_run_1691518350_chi2\", \"Module_0_run_1691518350_chi2.pdf\", 26, \"false\",\"1\",\"1691518350\",6)'")
# x
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_chi2\", \"Module_0_run_1691489947_chi2.pdf\", 26, \"false\",\"0\",\"1691489947\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\",  \"Module_1_run_1691489947_chi2\", \"Module_0_run_1691489947_chi2.pdf\", 26, \"false\",\"1\",\"1691489947\",6)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_chi2\", \"Module_0_run_1691491176_chi2.pdf\", 26, \"false\",\"0\",\"1691491176\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\",  \"Module_1_run_1691491176_chi2\", \"Module_0_run_1691491176_chi2.pdf\", 26, \"false\",\"1\",\"1691491176\",6)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_0_run_1691497635_chi2\", \"Module_0_run_1691497635_chi2.pdf\", 26, \"false\",\"0\",\"1691497635\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\",  \"Module_1_run_1691497635_chi2\", \"Module_0_run_1691497635_chi2.pdf\", 26, \"false\",\"1\",\"1691497635\",6)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_chi2\", \"Module_0_run_1691503538_chi2.pdf\", 26, \"false\",\"0\",\"1691503538\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\",  \"Module_1_run_1691503538_chi2\", \"Module_0_run_1691503538_chi2.pdf\", 26, \"false\",\"1\",\"1691503538\",6)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_chi2\", \"Module_0_run_1691505681_chi2.pdf\", 26, \"false\",\"0\",\"1691505681\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\",  \"Module_1_run_1691505681_chi2\", \"Module_0_run_1691505681_chi2.pdf\", 26, \"false\",\"1\",\"1691505681\",6)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_chi2\", \"Module_0_run_1691500205_chi2.pdf\", 26, \"false\",\"0\",\"1691500205\",6)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\",  \"Module_1_run_1691500205_chi2\", \"Module_0_run_1691500205_chi2.pdf\", 26, \"false\",\"1\",\"1691500205\",6)'")


### noise/area


exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_noiserPerarea\", \"Module_0_run_1691487584_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691487584\",7)'")

exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_noiserPerarea\", \"Module_1_run_1691487584_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691487584\",7)'")


exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_noiserPerarea\", \"Module_0_run_1691518350_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691518350\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\",  \"Module_1_run_1691518350_noiserPerarea\", \"Module_0_run_1691518350_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691518350\",7)'")

exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_noiserPerarea\", \"Module_0_run_1691489947_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691489947\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\",  \"Module_1_run_1691489947_noiserPerarea\", \"Module_0_run_1691489947_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691489947\",7)'")


exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_noiserPerarea\", \"Module_0_run_1691491176_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691491176\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\",  \"Module_1_run_1691491176_noiserPerarea\", \"Module_0_run_1691491176_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691491176\",7)'")

exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_0_run_1691497635_noiserPerarea\", \"Module_0_run_1691497635_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691497635\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\",  \"Module_1_run_1691497635_noiserPerarea\", \"Module_0_run_1691497635_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691497635\",7)'")

exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_noiserPerarea\", \"Module_0_run_1691503538_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691503538\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\",  \"Module_1_run_1691503538_noiserPerarea\", \"Module_0_run_1691503538_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691503538\",7)'")


exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_noiserPerarea\", \"Module_0_run_1691505681_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691505681\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\",  \"Module_1_run_1691505681_noiserPerarea\", \"Module_0_run_1691505681_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691505681\",7)'")

exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_noiserPerarea\", \"Module_0_run_1691500205_noiserPerarea.pdf\", 26, \"false\",\"0\",\"1691500205\",7)'")
exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\",  \"Module_1_run_1691500205_noiserPerarea\", \"Module_0_run_1691500205_noiserPerarea.pdf\", 26, \"false\",\"1\",\"1691500205\",7)'")




# #Calib_Files_Aug2023/level0_calib_params_run1691489947.txt
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_pedestal_DQMValid\", \"Module_0_run_1691487584_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691487584\",3)'")
                                     
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_pedestal_DQMValid\", \"Module_1_run_1691487584_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691487584\",3)'")
                                                                                                       
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691487584_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691487584\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691487584_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691487584\",4)'")

                                                                                                       
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691487584_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691487584\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691487584_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691487584\",5)'")



# # # execute root macro for TH2Poly
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_pedestal\", \"Module_0_run_1691518350_pedestal.pdf\", 26, \"false\",\"0\",\"1691518350\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\",  \"Module_1_run_1691518350_pedestal\", \"Module_0_run_1691518350_pedestal.pdf\", 26, \"false\",\"1\",\"1691518350\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_pedestal_GaussMean\", \"Module_0_run_1691518350_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691518350\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_1_run_1691518350_pedestal_GaussMean\", \"Module_0_run_1691518350_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691518350\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\",  \"Module_0_run_1691518350_Noise_GaussSigma\", \"Module_0_run_1691518350_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691518350\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\",  \"Module_1_run_1691518350_Noise_GaussSigma\", \"Module_0_run_1691518350_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691518350\",2)'")
                                                                                          
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_pedestal_DQMValid\", \"Module_0_run_1691518350_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691518350\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_1_run_1691518350_pedestal_DQMValid\", \"Module_1_run_1691518350_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691518350\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691518350_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691518350\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_1_run_1691518350_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691518350_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691518350\",4)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_0_run_1691518350_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691518350_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691518350\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691518350_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691518350.txt\", \"Module_1_run_1691518350_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691518350_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691518350\",5)'")

# # ## run2                           
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_pedestal\", \"Module_0_run_1691487584_pedestal.pdf\", 26, \"false\",\"0\",\"1691487584\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_pedestal\", \"Module_0_run_1691487584_pedestal.pdf\", 26, \"false\",\"1\",\"1691487584\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\",\"Module_0_run_1691487584_pedestal_GaussMean\", \"Module_0_run_1691487584_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691487584\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_1_run_1691487584_pedestal_GaussMean\", \"Module_0_run_1691487584_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691487584\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\", \"Module_0_run_1691487584_Noise_GaussSigma\", \"Module_0_run_1691487584_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691487584\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691487584_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691487584.txt\",\"Module_1_run_1691487584_Noise_GaussSigma\", \"Module_0_run_1691487584_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691487584\",2)'")

# ## run3                       
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_pedestal\", \"Module_0_run_1691489947_pedestal.pdf\", 26, \"false\",\"0\",\"1691489947\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_1_run_1691489947_pedestal\", \"Module_0_run_1691489947_pedestal.pdf\", 26, \"false\",\"1\",\"1691489947\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_pedestal_GaussMean\", \"Module_0_run_1691489947_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691489947\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_1_run_1691489947_pedestal_GaussMean\", \"Module_0_run_1691489947_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691489947\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_Noise_GaussSigma\", \"Module_0_run_1691489947_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691489947\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_1_run_1691489947_Noise_GaussSigma\", \"Module_0_run_1691489947_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691489947\",2)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_pedestal_DQMValid\", \"Module_0_run_1691489947_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691489947\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_1_run_1691489947_pedestal_DQMValid\", \"Module_1_run_1691489947_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691489947\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691489947_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691489947\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_1_run_1691489947_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691489947_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691489947\",4)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_0_run_1691489947_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691489947_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691489947\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691489947_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691489947.txt\", \"Module_1_run_1691489947_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691489947_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691489947\",5)'")




# # #run4
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_pedestal\", \"Module_0_run_1691491176_pedestal.pdf\", 26, \"false\",\"0\",\"1691491176\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\",\"Module_1_run_1691491176_pedestal\", \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_pedestal.pdf\", 26, \"false\",\"1\",\"1691491176\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_pedestal_GaussMean\", \"Module_0_run_1691491176_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691491176\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_1_run_1691491176_pedestal_GaussMean\", \"Module_0_run_1691491176_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691491176\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_Noise_GaussSigma\", \"Module_0_run_1691491176_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691491176\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_1_run_1691491176_Noise_GaussSigma\", \"Module_0_run_1691491176_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691491176\",2)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_pedestal_DQMValid\", \"Module_0_run_1691491176_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691491176\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_1_run_1691491176_pedestal_DQMValid\", \"Module_1_run_1691491176_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691491176\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691491176_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691491176\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_1_run_1691491176_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691491176_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691491176\",4)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_0_run_1691491176_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691491176_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691491176\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691491176_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691491176.txt\", \"Module_1_run_1691491176_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691491176_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691491176\",5)'")


# # ##run5

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\",  \"Module_0_run_1691497635_pedestal\", \"Module_0_run_1691497635_pedestal.pdf\", 26, \"false\",\"0\",\"1691497635\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_1_run_1691497635_pedestal\", \"Module_0_run_1691497635_pedestal.pdf\", 26, \"false\",\"1\",\"1691497635\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_0_run_1691497635_pedestal_GaussMean\", \"Module_0_run_1691497635_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691497635\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\",  \"Module_1_run_1691497635_pedestal_GaussMean\", \"Module_0_run_1691497635_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691497635\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\",  \"Module_0_run_1691497635_Noise_GaussSigma\", \"Module_0_run_1691497635_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691497635\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\",  \"Module_1_run_1691497635_Noise_GaussSigma\", \"Module_0_run_1691497635_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691497635\",2)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_0_run_1691497635_pedestal_DQMValid\", \"Module_0_run_1691497635_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691497635\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_1_run_1691497635_pedestal_DQMValid\", \"Module_1_run_1691497635_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691497635\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_0_run_1691497635_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691497635_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691497635\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_1_run_1691497635_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691497635_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691497635\",4)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_0_run_1691497635_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691497635_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691497635\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691497635_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691497635.txt\", \"Module_1_run_1691497635_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691497635_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691497635\",5)'")


# # ## run6

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\",  \"Module_0_run_1691500205_pedestal\", \"Module_0_run_1691500205_pedestal.pdf\", 26, \"false\",\"0\",\"1691500205\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_1_run_1691500205_pedestal\", \"Module_0_run_1691500205_pedestal.pdf\", 26, \"false\",\"1\",\"1691500205\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_pedestal_GaussMean\", \"Module_0_run_1691500205_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691500205\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_1_run_1691500205_pedestal_GaussMean\", \"Module_0_run_1691500205_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691500205\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_Noise_GaussSigma\", \"Module_0_run_1691500205_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691500205\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_1_run_1691500205_Noise_GaussSigma\", \"Module_0_run_1691500205_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691500205\",2)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_pedestal_DQMValid\", \"Module_0_run_1691500205_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691500205\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_1_run_1691500205_pedestal_DQMValid\", \"Module_1_run_1691500205_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691500205\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691500205_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691500205\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_1_run_1691500205_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691500205_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691500205\",4)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_0_run_1691500205_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691500205_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691500205\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691500205_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691500205.txt\", \"Module_1_run_1691500205_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691500205_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691500205\",5)'")




# # ## run7

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_pedestal\", \"Module_0_run_1691503538_pedestal.pdf\", 26, \"false\",\"0\",\"1691503538\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_1_run_1691503538_pedestal\", \"Module_0_run_1691503538_pedestal.pdf\", 26, \"false\",\"1\",\"1691503538\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_pedestal_GaussMean\", \"Module_0_run_1691503538_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691503538\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_1_run_1691503538_pedestal_GaussMean\", \"Module_0_run_1691503538_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691503538\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_Noise_GaussSigma\", \"Module_0_run_1691503538_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691503538\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_1_run_1691503538_Noise_GaussSigma\", \"Module_0_run_1691503538_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691503538\",2)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_pedestal_DQMValid\", \"Module_0_run_1691503538_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691503538\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_1_run_1691503538_pedestal_DQMValid\", \"Module_1_run_1691503538_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691503538\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691503538_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"1691503538\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_1_run_1691503538_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691503538_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"1691503538\",4)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_0_run_1691503538_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691503538_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"1691503538\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691503538_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691503538.txt\", \"Module_1_run_1691503538_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691503538_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"1691503538\",5)'")



# # ##run8

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_pedestal\", \"Module_0_run_1691505681_pedestal.pdf\", 26, \"false\",\"0\",\"1691505681\",0)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_1_run_1691505681_pedestal\", \"Module_0_run_1691505681_pedestal.pdf\", 26, \"false\",\"1\",\"1691505681\",0)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_pedestal_GaussMean\", \"Module_0_run_1691505681_pedestal_GaussMean.pdf\", 26, \"false\",\"0\",\"1691505681\",1)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\", \"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_1_run_1691505681_pedestal_GaussMean\", \"Module_0_run_1691505681_GaussMeanpedestal.pdf\", 26, \"false\",\"1\",\"1691505681\",1)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_Noise_GaussSigma\", \"Module_0_run_1691505681_Noise_GaussSigma.pdf\", 26, \"false\",\"0\",\"1691505681\",2)'")
# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",  \"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_1_run_1691505681_Noise_GaussSigma\", \"Module_0_run_1691505681_Noise_GaussSigma.pdf\", 26, \"false\",\"1\",\"1691505681\",2)'")


# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_pedestal_DQMValid\", \"Module_0_run_1691505681_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"run_1691505681\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_1_run_1691505681_pedestal_DQMValid\", \"Module_1_run_1691505681_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"run_1691505681\",3)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_pedestal_GaussMean_DQMValid\", \"Module_0_run_1691505681_pedestal_DQMValid.pdf\", 26, \"false\",\"0\",\"run_1691505681\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_1_run_1691505681_pedestal_GaussMean_DQMValid\", \"Module_1_run_1691505681_pedestal_DQMValid.pdf\", 26, \"false\",\"1\",\"run_1691505681\",4)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_0_run_1691505681_Noise_GaussSigma_DQMValid\", \"Module_0_run_1691505681_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"0\",\"run_1691505681\",5)'")

# exe("root -l -b -q th2poly.C'(\"./data/hexagons.root\", \"./run_1691505681_pedestal.txt\",\"./Calib_Files_Aug2023/level0_calib_params_run1691505681.txt\", \"Module_1_run_1691505681_Noise_GaussSigma_DQMValid\", \"Module_1_run_1691505681_Noise_GaussSigma_DQMValid.pdf\", 26, \"false\",\"1\",\"run_1691505681\",5)'")



# ## run9
# ## no run beyond this :) ##
