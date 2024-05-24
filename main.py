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
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

# create cadquery model
cad_model = cadquery.Workplane("xy").box(L, W, W)
cad_exchange = cad_model.toOCC()

# mesh it 
gmsh.initialize()
gmsh.model.occ.importShapesNativePointer(cad_exchange._address())
gmsh.model.occ.synchronize()


# simulation
gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)
gmsh.finalize()

# analysis
# saving
