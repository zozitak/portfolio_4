# main.py

import cadquery
import gmsh
from dolfinx.io import gmshio
from dolfinx import mesh, fem, plot, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import numpy as np
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

gmsh.model.mesh.generate(3)

# simulation
mu = 1
rho = 1
delta = W / L
gamma = 0.4 * delta**2
beta = 1.25
lambda_ = beta
g = gamma

gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD
domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=3)
gmsh.finalize()
V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

def clamped_boundary(x):
    return np.isclose(x[0], 0)

fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)

u_D = np.array([0, 0, 0], dtype=default_scalar_type)
bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets), V)

T = fem.Constant(domain, default_scalar_type((0, 0, 0)))

ds = ufl.Measure("ds", domain=domain)

def epsilon(u):
    return ufl.sym(ufl.grad(u))  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, default_scalar_type((0, 0, -rho * g)))
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

# analysis
#Von Mises Stress
s = sigma(uh) - 1. / 3 * ufl.tr(sigma(uh)) * ufl.Identity(len(uh))
von_Mises = ufl.sqrt(3. / 2 * ufl.inner(s, s))

V_von_mises = fem.functionspace(domain, ("DG", 0))
stress_expr = fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
stresses = fem.Function(V_von_mises)
stresses.interpolate(stress_expr)

uh_magnitude = ufl.sqrt(ufl.dot(uh, uh))
print('min/max uh:',
      uh_magnitude.vector().array().min(),
      uh_magnitude.vector().array().max())

# saving
with io.XDMFFile(domain.comm, "displacements.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    uh.name = "displacements"
    xdmf.write_function(uh)

with io.XDMFFile(domain.comm, "von_mises_stress.xdmf", "w") as file:
    file.write_mesh(domain)
    uh.name = "von_mises_stress"
    file.write_function(stresses)