import festim as F  # using festim2
from dolfinx.io import gmshio
from dolfinx import plot
from mpi4py import MPI
import numpy as np
import pyvista

# LOAD AND READ GMSH MESH

model_rank = 0
mesh, cell_tags, facet_tags = gmshio.read_from_msh(
    "meshing/probe_breeder.msh", MPI.COMM_WORLD, model_rank, gdim=3
)

print(f"Cell tags: {np.unique(cell_tags.values)}")
print(f"Facet tags: {np.unique(facet_tags.values)}")

# DEFINE & INITIALIZE MODEL

my_model = F.HydrogenTransportProblem()

my_model.mesh = F.Mesh(mesh)

material = F.Material(D_0=1, E_D=0)

# SET DOMAINS

vol = F.VolumeSubdomain(id=1, material=material)

inlet = F.SurfaceSubdomain(id=1)
outlet = F.SurfaceSubdomain(id=2)
wall = F.SurfaceSubdomain(id=3)
probe = F.SurfaceSubdomain(id=4)

# pass the meshtags to the model directly
my_model.facet_meshtags = facet_tags
my_model.volume_meshtags = cell_tags

my_model.subdomains = [inlet, outlet, wall, probe, vol]

H = F.Species("H")
my_model.species = [H]

# SET TEMP AND BOUNDARY CONDITIONS

my_model.temperature = 400

my_model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=inlet, value=1, species=H),
    F.FixedConcentrationBC(subdomain=probe, value=0, species=H),
    F.FixedConcentrationBC(subdomain=outlet, value=0.0, species=H),
    F.ParticleFluxBC(subdomain=wall, value=0.0, species=H),
]

# SETTINGS AND EXPORTS

my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

results_folder = "festim_results"

my_model.exports = [F.VTXSpeciesExport(f"{results_folder}/H.bp", field=H)]

# INITIALISE AND RUN

my_model.initialise()
my_model.run()

# VISUALIZE WITH PYVISTA

hydrogen_concentration = H.solution

topology, cell_types, geometry = plot.vtk_mesh(hydrogen_concentration.function_space)
u_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
u_grid.point_data["c"] = hydrogen_concentration.x.array.real
u_grid.set_active_scalars("c")
u_plotter = pyvista.Plotter()
u_plotter.add_mesh(u_grid, show_edges=True, opacity=0.5)

if not pyvista.OFF_SCREEN:
    u_plotter.show()
else:
    figure = u_plotter.screenshot("concentration.png")
