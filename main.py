# main.py

import cadquery
import gmsh
from dolfinx.io import gmshio
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import sys

from cqplugin.export_step import export_step_ap214

# Scaled variables
L = 1.0; W = 0.2

# create cadquery model
cad_model = cadquery.Workplane("XY").move(L/2).box(L, W, W)
cad_exchange = cad_model.toOCC()
#cadquery.exporters.export(cad_model,"test.step")

gmsh.initialize()
gmsh.model.occ.importShapesNativePointer(cad_exchange._address())
gmsh.model.occ.synchronize()

# check for duplicates // sometimes cadquery doesnt flattens the model tree 
gmsh.model.occ.removeAllDuplicates()
# transfinite 1,2,3,4 and 5,6,7,8 to 9 progression 1
for curve in [1,2,3,4,5,6,7,8]:
    gmsh.model.mesh.setTransfiniteCurve(curve,9) 
# transfinite 9,10,11,12 to 41 progression 1 
for curve in [9,10,11,12]:
    gmsh.model.mesh.setTransfiniteCurve(curve,41)
# transfinite surface 1,2,3,4,5,6
for surface in [1,2,3,4,5,6]:
    gmsh.model.mesh.setTransfiniteSurface(surface)
# recombine 1,2,3,4,5,6
for surface in [1,2,3,4,5,6]:
    gmsh.model.mesh.setRecombine(2,surface)
# transfinite volume[1] = [1,2,4,3,5,6,8,7]
for volume in [1]:
    gmsh.model.mesh.setTransfiniteVolume(volume,[1,2,4,3,5,6,8,7])

# simulation
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma
gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=3)


gmsh.finalize()

# analysis
# saving
