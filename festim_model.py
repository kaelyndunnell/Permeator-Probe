import festim as F  # using festim2
import numpy as np
from dolfinx import fem
import ufl
from dolfinx import cpp as _cpp
from openfoam_to_festim import read_openfoam_data
from dolfinx.log import set_log_level, LogLevel


def evaluate_stabalisation_term(mesh, u, delta):
    """See more at https://www.comsol.com/blogs/understanding-stabilization-methods"""

    # evaluate Cell size
    tdim = mesh.topology.dim
    num_cells = mesh.topology.index_map(tdim).size_local
    cells = np.arange(num_cells, dtype=np.int32)
    mesh_ = _cpp.mesh.Mesh_float64(
        mesh.comm, mesh.topology._cpp_object, mesh.geometry._cpp_object
    )
    h = _cpp.mesh.h(mesh_, tdim, cells)
    V0 = fem.functionspace(mesh, ("DG", 0))
    h_as_function = fem.Function(V0)
    h_as_function.x.array[:] = h

    # Compute magnitude of velocity
    v_mag = ufl.sqrt(ufl.dot(u, u))

    D_art = delta * v_mag * h_as_function

    return D_art


# markers for gmsh TODO: do not make this repetitive
inlet_marker = 1
outlet_marker = 2
wall_marker = 3
probe_marker = 4

p, u, mesh, nut, facet_meshtags, volume_meshtags = read_openfoam_data(
    "OpenFOAM/turbulent-case/probe.foam", final_time=300
)

# DEFINE & INITIALIZE MODEL

print("Building FESTIM model...")

my_model = F.HydrogenTransportProblem()

my_model.mesh = F.Mesh(mesh)
my_model.facet_meshtags = facet_meshtags
my_model.volume_meshtags = volume_meshtags

breeder_temperature = 603.15  # K
D_0_PbLi = 4.03e-08  # m2/s
E_D_PbLi = 0.20  # eV

D_diff = D_0_PbLi * ufl.exp(-E_D_PbLi / (F.k_B * breeder_temperature))
D_art = evaluate_stabalisation_term(mesh=mesh, u=u, delta=0.1)

D_expr = D_diff + D_art
V = fem.functionspace(mesh, ("CG", 1))
D_pbli = fem.Function(V)
D_pbli.interpolate(fem.Expression(D_expr, V.element.interpolation_points()))
material = F.Material(D=D_pbli)

# SET DOMAINS

vol = F.VolumeSubdomain(id=1, material=material)

# use same tags as gmsh markers
inlet = F.SurfaceSubdomain(id=inlet_marker)
outlet = F.SurfaceSubdomain(id=outlet_marker)
wall = F.SurfaceSubdomain(id=wall_marker)
probe = F.SurfaceSubdomain(id=probe_marker)

my_model.subdomains = [inlet, outlet, wall, probe, vol]

H = F.Species("H")
my_model.species = [H]

# SET TEMP AND BOUNDARY CONDITIONS

my_model.temperature = breeder_temperature  # K

my_model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=inlet, value=1e16, species=H),
    F.FixedConcentrationBC(subdomain=probe, value=0, species=H),
]

advection = F.AdvectionTerm(velocity=u, subdomain=vol, species=H)
my_model.advection_terms = [advection]

# SETTINGS AND EXPORTS

dt = F.Stepsize(
    initial_value=1, growth_factor=1.05, cutback_factor=0.9, target_nb_iterations=5
)

my_model.settings = F.Settings(
    atol=1e04,
    rtol=1e-10,
    transient=True,
    final_time=200,
    stepsize=dt,
)

results_folder = "festim_results"

my_model.exports = [F.VTXSpeciesExport(filename=f"{results_folder}/H.bp", field=H)]

# INITIALISE AND RUN

my_model.initialise()
set_log_level(LogLevel.INFO)
my_model.run()
