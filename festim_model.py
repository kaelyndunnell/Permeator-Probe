import festim as F  # using festim2
from dolfinx.io import gmshio
from dolfinx import plot
from mpi4py import MPI
import numpy as np
import pyvista
from openfoam_to_festim import read_openfoam_data

# from meshing.cad_to_gmsh import inlet_marker, outlet_marker, wall_marker, probe_marker TODO: clean this up

# markers for gmsh TODO: do not make this repetitive

inlet_marker = 1
outlet_marker = 2
wall_marker = 3
probe_marker = 4

p, u = read_openfoam_data(final_time=100)

print("Building FESTIM model...")

# LOAD AND READ GMSH MESH

model_rank = 0
mesh, cell_tags, facet_tags = gmshio.read_from_msh(
    "OpenFOAM/probe-case/probe_breeder.msh", MPI.COMM_WORLD, model_rank, gdim=3
)

print(f"Cell tags: {np.unique(cell_tags.values)}")
print(f"Facet tags: {np.unique(facet_tags.values)}")

# DEFINE & INITIALIZE MODEL

my_model = F.HydrogenTransportProblem()

my_model.mesh = F.Mesh(mesh)

material = F.Material(D_0=1, E_D=0)

# SET DOMAINS

vol = F.VolumeSubdomain(id=1, material=material)

# use same tags as gmsh markers
inlet = F.SurfaceSubdomain(id=inlet_marker)
outlet = F.SurfaceSubdomain(id=outlet_marker)
wall = F.SurfaceSubdomain(id=wall_marker)
probe = F.SurfaceSubdomain(id=probe_marker)

# pass the meshtags to the model directly
my_model.facet_meshtags = facet_tags
my_model.volume_meshtags = cell_tags

my_model.subdomains = [inlet, outlet, wall, probe, vol]

H = F.Species("H")
my_model.species = [H]

# SET TEMP AND BOUNDARY CONDITIONS

my_model.temperature = 603.15  # K

my_model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=inlet, value=1, species=H),
    F.FixedConcentrationBC(subdomain=probe, value=0, species=H),
    F.FixedConcentrationBC(subdomain=outlet, value=0.0, species=H),
    F.ParticleFluxBC(subdomain=wall, value=0.0, species=H),
]

advection = F.AdvectionTerm(velocity=u, subdomain=vol, species=H)
my_model.advection_terms = [advection]

# SETTINGS AND EXPORTS

dt = F.Stepsize(
    initial_value=10, growth_factor=1.1, cutback_factor=0.9, target_nb_iterations=5
)

my_model.settings = F.Settings(
    atol=1e-10,
    rtol=1e-10,
    transient=True,
    final_time=200,
    stepsize=dt,
)

results_folder = "festim_results"

my_model.exports = [F.VTXSpeciesExport(filename=f"{results_folder}/H.bp", field=H)]

# INITIALISE AND RUN

my_model.initialise()
my_model.run()

# # VISUALIZE WITH PYVISTA

# hydrogen_concentration = H.solution

# topology, cell_types, geometry = plot.vtk_mesh(hydrogen_concentration.function_space)
# u_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
# u_grid.point_data["c"] = hydrogen_concentration.x.array.real
# u_grid.set_active_scalars("c")
# u_plotter = pyvista.Plotter()
# u_plotter.add_mesh(u_grid, show_edges=True, opacity=0.5)

# if not pyvista.OFF_SCREEN:
#     u_plotter.show()
# else:
#     figure = u_plotter.screenshot("concentration.png")
