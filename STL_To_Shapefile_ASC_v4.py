# -*- coding: utf-8 -*-



#This code was rewritten by Nicolas Schmid, it runs with python 3. You need to install different packages.

#I recommand using anaconda to create virtual environnements and install the packages



#the only parameter which you should change is the grid size. It give the number of meters between two points of your

#grid on the plot. Recommended between 15 and 100



############Set the grid size ####################

gridsize = 1



sigma = 2.5*50/gridsize #the smoothing parameters depend on the grid size, you can change this but be carefull!

truncate=3.0*50/gridsize

##################################################



#install anaconda: https://www.anaconda.com

"""

conda create --name MY_PYCASTER_ENV

conda activate MY_PYCASTER_ENV

"""



#Then install the following packages:

"""

conda install matplotlib

conda install vtk

pip install numpy-stl  

conda install shapely

conda install fiona

pip install scipy

"""



#Somehow it couldn't find rayCaster in pycaster. My solution was to save pycaster.py (the package) in the same folder as this script and name it mypycaster.py

from tkinter import *

from tkinter import filedialog

import numpy as np

from stl import mesh

import matplotlib.pyplot as plt

import matplotlib.mlab as ml

import vtk

import logging as log

import math

from shapely import geometry

import os

import fiona

import scipy as sp

import scipy.ndimage





class rayCaster(object):

    def __init__(self, mesh):

        self.mesh = mesh  # set the 'mesh'

        self.caster = None

        self._initCaster()  # create a caster

    @classmethod



    def fromSTL(cls, filenameSTL, scale=1.0):

        readerSTL = vtk.vtkSTLReader()  # create a 'vtkSTLReader' object

        readerSTL.SetFileName(

            filenameSTL)  # set the .stl filename in the reader

        readerSTL.Update()  # 'update' the reader i.e. read the .stl file

        polydata = readerSTL.GetOutput()

        # If there are no points in 'vtkPolyData' something went wrong

        if polydata.GetNumberOfPoints() == 0:

            raise ValueError(

                "No point data could be loaded from '" + filenameSTL)

            return None

        # Create a new 'rayCaster' with the mesh loaded from the .stl file

        rT = cls(polydata)

        #if a non-unit 'scale' is given then scale the polydata

        if scale != 1.0:

            rT.scaleMesh(scale)

        #return the 'rayCaster' object

        return rT



    def scaleMesh(self, scale):

        #create a 'vtkTransform' object

        transform = vtk.vtkTransform()

        #create a 'vtkTransformPolyDataFilter' object

        transformFilter = vtk.vtkTransformPolyDataFilter()

        #set the scale transformation

        transform.Scale(scale, scale, scale)  #assuming uniform scale

        #set the 'mesh' as the input of the 'vtkTransformPolyDataFilter'

        transformFilter.SetInput(self.mesh)

        #set 'transform' as the transformation of the 'vtkTransformPolyDataFilter'

        transformFilter.SetTransform(transform)

        #update to apply the transformation

        transformFilter.Update()

        #replace the existing 'mesh' with the scaled one

        self.mesh = transformFilter.GetOutput()

        #update the caster since the 'mesh' changed

        self._updateCaster()



    def _initCaster(self):

        #create a 'vtkOBBTree' object

        self.caster = vtk.vtkOBBTree()

        #set the 'mesh' as the caster's dataset

        self.caster.SetDataSet(self.mesh)

        #build a caster locator

        self.caster.BuildLocator()



    def _updateCaster(self):

        self.caster.SetDataSet(self.mesh)

        self.caster.BuildLocator()



    def castRay(self, pointRaySource, pointRayTarget):

        #create a 'vtkPoints' object to store the intersection points

        pointsVTKintersection = vtk.vtkPoints()

        #perform ray-casting (intersect a line with the mesh)

        code = self.caster.IntersectWithLine(pointRaySource,

                                             pointRayTarget,

                                             pointsVTKintersection, None)



        # Interpret the 'code'. If "0" is returned then no intersections points

        # were found so return an empty list

        if code == 0:

            log.info(

                "No intersection points found for 'pointRaySource': " + str(

                    pointRaySource) + " and 'pointRayTarget': " + str(

                    pointRayTarget))

            return []

        # If code == -1 then 'pointRaySource' lies outside the surface

        elif code == -1:

            log.info("The point 'pointRaySource': " + str(

                pointRaySource) + "lies inside the surface")

        #get the actual data of the intersection points (the point tuples)

        pointsVTKIntersectionData = pointsVTKintersection.GetData()

        #get the number of tuples

        noPointsVTKIntersection = pointsVTKIntersectionData.GetNumberOfTuples()

        #create an empty list that will contain all list objects

        pointsIntersection = []

        # Convert the intersection points to a list of list objects.

        for idx in range(noPointsVTKIntersection):

            _tup = pointsVTKIntersectionData.GetTuple3(idx)

            pointsIntersection.append(_tup)

        #return the list of list objects

        return pointsIntersection



def smooth_with_nan(my_grid,nan_mask,sigma,truncate):

    #this takes a 2d grid containing NaNs as input and applies a gaussian filter on it to make the data smoother

    #reference: https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python

    U = my_grid

    V=U.copy()

    V[np.isnan(U)]=0

    VV=sp.ndimage.gaussian_filter(V,sigma=sigma,truncate=truncate)

    W=0*U.copy()+1

    W[np.isnan(U)]=0

    WW=sp.ndimage.gaussian_filter(W,sigma=sigma,truncate=truncate)

    Z=VV/WW

    Z[nan_mask]=np.nan

    return Z



##############    Tkinter stuff   ##########################################

stlfile = filedialog.askopenfile(title="Select your .stl file")# /DATA



##############    Ray Casting   ##########################################

Mycaster = rayCaster.fromSTL(stlfile.name, scale=1)

my_mesh =mesh.Mesh.from_file(stlfile.name)

data = np.array(my_mesh)

xcoord = []; ycoord = []; zcoord = []

for el in data:

    coord = el[0], el[1], el[2]

    xcoord.append(el[0])

    ycoord.append(el[1])

#make new regular grid using gridsize

xcorner = int(round(min(xcoord),-1))

ycorner = int(round(min(ycoord),-1))

nbcellx = int(round((max(xcoord)-min(xcoord))/gridsize,0))

nbcelly = int(round((max(ycoord)-min(ycoord))/gridsize,0))

xi = xcorner + np.arange(0,nbcellx)*gridsize

yi = ycorner + np.arange(0,nbcelly)*gridsize

[colx, coly] = np.meshgrid(xi, yi)

final_grid_top = np.empty((nbcelly,nbcellx,3))

final_grid_top[:] = np.nan

final_grid_bottom = np.empty((nbcelly,nbcellx,3))

final_grid_bottom[:]=np.nan

for i in range(len(xi)):

    for j in range(len(yi)):

        PointsIntersection = Mycaster.castRay([xi[i],yi[j],10000],[xi[i],yi[j],-10000]) 

        if PointsIntersection:

            final_grid_top[j][i] = PointsIntersection[0] #point intersection on top

            final_grid_bottom[j][i] = PointsIntersection[-1] #point intersection bottom

thickness_grid = final_grid_top[:,:,2]-final_grid_bottom[:,:,2]

mask_thickness = np.isnan(thickness_grid)

thickness_grid[thickness_grid == 0] = np.nan #replace 0 with nan

altitude_grid = final_grid_top[:,:,2]

mask_altitude = np.isnan(altitude_grid)





############ integrate volume of 3D shape with the thickness grid#################################################

volm3=0

for row in thickness_grid:

    for element in row:

        if element!=np.nan and not math.isnan(element):

            volm3+=element*gridsize*gridsize

print("vol in m3:",volm3)

print("vol in km3",volm3/1000000000)



################# get levels for altitude plot ###############################

zmax=np.nanmax(altitude_grid)

print("zmax",zmax)

zmin=np.nanmin(altitude_grid)

print("zmin",zmin)

levelsZ = []

height = zmin-zmin%5

while height < zmax:

    levelsZ.append(height)

    height +=100

if (zmax-1 >height-100):

    levelsZ.append(zmax-1)

levelsZ.append(height)

#levels for thickness plot

levelsZZ = [0,10,50,150,10000]



############# make the grids smoother with a gaussian filter

end_thickness=smooth_with_nan(thickness_grid,mask_thickness,sigma=sigma,truncate=truncate)

end_altitude=smooth_with_nan(altitude_grid,mask_altitude,sigma=sigma,truncate=truncate)



""" plots to check what the smoothing does on our data



plt.figure()

plt.imshow(altitude_grid)

plt.colorbar()

plt.title("altitude before smoothing")

plt.figure()



plt.imshow(end_altitude)

plt.colorbar()

plt.title("altitude after smoothing")



plt.figure()

plt.imshow(thickness_grid)

plt.colorbar()

plt.title("thickness before smoothing")

plt.figure()



plt.imshow(end_thickness)

plt.colorbar()

plt.title("thickness after smoothing")

plt.figure()

"""



##############  Plot ##########################################



#Altitude plot 

CF = plt.contourf(colx,coly,end_altitude, levels = levelsZ ) #Poly

CS = plt.contour(colx,coly,end_altitude, levels=levelsZ , colors='black',linewidths=0.8)# line



plt.colorbar(CF)

plt.title('Altitude plot       Resolution : '+str(gridsize)); plt.xlabel('Coord X (m)'); plt.ylabel('Coord Y (m)')

plt.figure()

#thickness plot

CZF = plt.contourf(colx, coly,end_thickness, levels = levelsZZ,colors=('red','green','blue','orange')) #Poly

CZ = plt.contour(colx,coly,end_thickness, levels=levelsZZ, colors='black',linewidths=0.8)#thickness



plt.title('Thickness plot        Resolution : '+str(gridsize)); plt.xlabel('Coord X (m)'); plt.ylabel('Coord Y (m)')

plt.colorbar(CZF)

plt.show()





##############  Save as shapefile ##########################################



#to store the polygons I followed the following example

# https://chrishavlin.wordpress.com/2019/08/20/shapefiles-in-python-converting-contours-to-shapes/



#create the directories

separator='/'

output_file_name = stlfile.name[:-4] +'_output_res'+str(gridsize)

file_name = stlfile.name[:-4].split("/")[-1]

if not os.path.exists(output_file_name):

    os.makedirs(output_file_name)



################store altitude polygon#################

lvl_lookup = dict(zip(CF.collections, CF.levels))

PolyList=[]

for col in CF.collections:

    z=lvl_lookup[col] # the value of this level

    for contour_path in col.get_paths():

        # create the polygon for this level

        for ncp,cp in enumerate(contour_path.to_polygons()):

            lons = cp[:,0]

            lats = cp[:,1]

            new_shape = geometry.Polygon([(i[0], i[1]) for i in zip(lons,lats)])

            if ncp == 0:                

                poly = new_shape # first shape

            else:

                poly = poly.difference(new_shape) # Remove the holes

            PolyList.append({'poly':poly,'props':{'z': z}})



schema = {'crs': 'epsg:2056','geometry': 'Polygon','properties': {'z': 'float'}}

with fiona.collection(output_file_name+separator+file_name+'_height_poly.shp', "w", "ESRI Shapefile", schema) as output:

    for p in PolyList:

        output.write({'properties': p['props'],

            'geometry': geometry.mapping(p['poly'])})



################store thickness polygon#################

lvl_lookup = dict(zip(CZF.collections, CZF.levels))

PolyList=[]

for col in CZF.collections:

    z=lvl_lookup[col] # the value of this level

    for contour_path in col.get_paths():

        # create the polygon for this level

        for ncp,cp in enumerate(contour_path.to_polygons()):

            lons = cp[:,0]

            lats = cp[:,1]

            new_shape = geometry.Polygon([(i[0], i[1]) for i in zip(lons,lats)])

            if ncp == 0:                

                poly = new_shape # first shape

            else:

                poly = poly.difference(new_shape) # Remove the holes

            PolyList.append({'poly':poly,'props':{'z': z}})



schema = {'geometry': 'Polygon','properties': {'z': 'float'}}

with fiona.collection(output_file_name+separator+file_name+'_thick_poly.shp', "w", "ESRI Shapefile", schema) as output:

    for p in PolyList:

        output.write({'properties': p['props'],

            'geometry': geometry.mapping(p['poly'])})





################store altitude line#################

lvl_lookup = dict(zip(CS.collections, CS.levels))

LineList = []

for col in CS.collections:

    z = lvl_lookup[col]  # The value of this level

    for contour_path in col.get_paths():

        # Create the LineString for this level

            # for each segment of that section

        line = contour_path.vertices.tolist()

        LineList.append({'line': line, 'props': {'z': z}})



schema = {'geometry': 'LineString', 'properties': {'z': 'float'}}

with fiona.open(output_file_name+separator+file_name+'_height_line.shp', "w", "ESRI Shapefile", schema) as output:

    for line in LineList:

        geom = geometry.LineString(line['line'])

        output.write({'properties': line['props'], 'geometry': geometry.mapping(geom)})





################store altitude line#################

lvl_lookup = dict(zip(CZ.collections, CZ.levels))

LineList = []

for col in CZ.collections:

    z = lvl_lookup[col]  # The value of this level

    for contour_path in col.get_paths():

        # Create the LineString for this level

        line = contour_path.vertices.tolist()

        LineList.append({'line': line, 'props': {'z': z}})

schema = {'geometry': 'LineString', 'properties': {'z': 'float'}}

with fiona.open(output_file_name+separator+file_name+'_thick_line.shp', "w", "ESRI Shapefile", schema) as output:

    for line in LineList:

        geom = geometry.LineString(line['line'])

        output.write({'properties': line['props'], 'geometry': geometry.mapping(geom)})



########### save ascii file #################################################



info_str = "ncols\t\t"+str(nbcellx) +'\nnrows\t\t'+str(nbcelly)+'\nxllcorner\t' +str(xcorner) +'\nyllcorner\t' + str(ycorner) + '\ncellsize\t' +str(gridsize) + '\nnodata value\t-9999.000'

new_array=np.nan_to_num(final_grid_top[:,:,2],nan=-9999)

final_array = np.empty(new_array.shape)

for i in range(nbcelly):

    final_array[nbcelly-1-i,:] =new_array[i,:]

np.savetxt(output_file_name + separator +file_name+ '_altitude.asc',final_array,delimiter=' ',fmt="%.3f",header=info_str,comments='')

















