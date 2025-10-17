from ufl import exp, dot, sqrt
from dolfinx.fem import Function, functionspace, Expression
from foam2dolfinx import OpenFOAMReader
from dolfinx.io import gmshio
import dolfinx
import festim as F
import numpy as np
import h_transport_materials as htm
from openfoam_to_festim import read_openfoam_data
from dolfinx import cpp as _cpp
from mpi4py import MPI

p, u, mesh, nut = read_openfoam_data(
    "OpenFOAM/turbulent-case/case.foam", final_time=3.0
)
if nut is None:
    raise ValueError(
        "Turbulent viscosity field 'nut' is required for Peclet number calculation."
    )

T = 603.15  # K

D_lipb = htm.diffusivities.filter(material="lipb").mean()
D_0 = D_lipb.pre_exp.magnitude
E_D = D_lipb.act_energy.magnitude

V0 = functionspace(mesh, ("DG", 0))
V1 = functionspace(mesh, ("DG", 1))

# Compute magnitude of velocity
v_mag = sqrt(dot(u, u))

# Compute temperature-dependent diffusivity D(x)
D_diff = D_0 * exp(-E_D / (F.k_B * T))
Schmidt = 0.5
D_turb = nut / Schmidt

D_expr = D_diff + D_turb

# evaluate Cell size
tdim = mesh.topology.dim
num_cells = mesh.topology.index_map(tdim).size_local
cells = np.arange(num_cells, dtype=np.int32)
mesh_ = _cpp.mesh.Mesh_float64(
    mesh.comm, mesh.topology._cpp_object, mesh.geometry._cpp_object
)
h = _cpp.mesh.h(mesh_, tdim, cells)

h_as_function = Function(V0)
h_as_function.x.array[:] = h

# evaluate Peclet number
Pe_expr = v_mag * h_as_function / D_expr

Pe_local = Function(V1)
Pe_local.interpolate(Expression(Pe_expr, V1.element.interpolation_points()))

writer = dolfinx.io.VTXWriter(MPI.COMM_WORLD, "Pe_local.bp", Pe_local, "BP5")
writer.write(t=0.0)
