# Motivated by this example:
# https://jorgensd.github.io/dolfinx-tutorial/chapter1/fundamentals.html 
# ... but hopefully even simpler!

from mpi4py import MPI
from dolfinx import mesh
from dolfinx.fem import FunctionSpace
from dolfinx import fem
import numpy as np
import ufl
from petsc4py.PETSc import ScalarType
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

# Define the domain to be the unit square.
nx, ny = 8,8
domain = mesh.create_unit_square(MPI.COMM_WORLD, nx, ny, mesh.CellType.quadrilateral)

# Create a Continuous Galerkin function space
V = FunctionSpace(domain, ("CG", 1))

# Define a Dirichlet Boundary condition 
uD = fem.Function(V)
uD.interpolate(lambda x: 1 + x[0]**2 + 2 * x[1]**2)

# Create facet to cell connectivity required to determine boundary facets
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)
boundary_dofs = fem.locate_dofs_topological(V, fdim, boundary_facets)
bc = fem.dirichletbc(uD, boundary_dofs)

# Define the trial functions and test functions
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# Define the terms in the variational formulation
f = fem.Constant(domain, ScalarType(-6))
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = ufl.inner(f, v) * ufl.dx

# Formulate the problem and solve it
problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

# interpolate the solution onto a regular grid
interpolated_solution = interp2d ( domain.geometry.x[:,0], domain.geometry.x[:,1],  np.real(uh.x.array ) )
x = np.linspace(0,1,nx)
y = np.linspace(0,1,ny)
X,Y = np.meshgrid(x,y)
gridded_interpolated_solution = interpolated_solution(x,y).reshape(nx,ny)

# Plot using matplotlib
fig,ax=plt.subplots()
plt.pcolor(x,y,gridded_interpolated_solution)
ax.set_aspect('equal')