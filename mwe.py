import gmsh
import dolfinx
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
import festim as F
from dolfinx import fem
import ufl
from basix.ufl import element
import numpy as np
from dolfinx import cpp as _cpp
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


gmsh.initialize()
gmsh.model.add("mwe")

r_inner = 0.1  # fluid radius (m)
r_tube = 0.105  # tube wall radius (m)
length = 0.4  # cylinder length (m)

factory = gmsh.model.occ

fluid = factory.addRectangle(0, 0, 0, length, r_inner)
tube = factory.addRectangle(0, 0, 0, length, r_tube)
walls, interface = factory.cut(
    [(2, tube)], [(2, fluid)], removeObject=True, removeTool=False
)

factory.synchronize()
interface_curves = [s[1] for s in interface[0]]
interface_tag = gmsh.model.addPhysicalGroup(1, interface_curves)
gmsh.model.setPhysicalName(1, interface_tag, "interface")

fluid_tag = gmsh.model.addPhysicalGroup(2, [fluid])
gmsh.model.setPhysicalName(2, fluid_tag, "fluid")

wall_tag = gmsh.model.addPhysicalGroup(2, [walls[0][1]])
gmsh.model.setPhysicalName(2, wall_tag, "wall")

surfaces = gmsh.model.getEntities(dim=1)

inlet_surfaces = []
outlet_surfaces = []

for s in surfaces:
    com = gmsh.model.occ.getCenterOfMass(s[0], s[1])
    # inlet at x=0, outlet at x=length
    if abs(com[0]) < 1e-6:
        inlet_surfaces.append(s[1])
    elif abs(com[0] - length) < 1e-6:
        outlet_surfaces.append(s[1])

if inlet_surfaces:
    inlet_tag = gmsh.model.addPhysicalGroup(1, inlet_surfaces)
    gmsh.model.setPhysicalName(1, inlet_tag, "inlet")

if outlet_surfaces:
    outlet_tag = gmsh.model.addPhysicalGroup(1, outlet_surfaces)
    gmsh.model.setPhysicalName(1, outlet_tag, "outlet")

gmsh.model.mesh.generate(2)

gmsh.write("mwe.msh")

mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=2)
my_mesh = mesh_data.mesh

gmsh.finalize()

my_model = F.HydrogenTransportProblemDiscontinuous()
my_model.mesh = F.Mesh(my_mesh)
my_model.facet_meshtags = mesh_data.facet_tags
my_model.volume_meshtags = mesh_data.cell_tags


el = element("Lagrange", my_mesh.topology.cell_name(), 2, shape=(my_mesh.geometry.dim,))


V = dolfinx.fem.functionspace(my_model.mesh.mesh, el)

velocity = dolfinx.fem.Function(V)

velocity.interpolate(lambda x: (-100 * x[1] * (x[1] - 1), np.full_like(x[0], 0.0)))


D_diff = 2 * ufl.exp(-1 / (F.k_B * 500))

# add stabilization term for diffusion
D_art = evaluate_stabalisation_term(mesh=my_mesh, u=velocity, delta=0.1)

D_expr = D_diff + D_art
V = fem.functionspace(my_mesh, ("CG", 1))
D_fluid = fem.Function(V)
D_fluid.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

dummy_fluid = F.Material(D=D_fluid, K_S_0=1, E_K_S=1)
dummy_tube = F.Material(D_0=2, E_D=1, K_S_0=2, E_K_S=2)

inlet = F.SurfaceSubdomain(id=inlet_tag)
outlet = F.SurfaceSubdomain(id=outlet_tag)
fluid = F.VolumeSubdomain(id=fluid_tag, material=dummy_fluid)
tube = F.VolumeSubdomain(id=wall_tag, material=dummy_tube)

my_model.subdomains = [inlet, outlet, fluid, tube]

my_model.surface_to_volume = {inlet: fluid, outlet: fluid, interface_tag: fluid}
H = F.Species("H", subdomains=my_model.volume_subdomains)
my_model.species = [H]
my_model.interfaces = [
    F.Interface(id=interface_tag, subdomains=[fluid, tube], penalty_term=100)
]

my_model.temperature = 500

advection_term = F.AdvectionTerm(
    velocity=velocity,
    subdomain=fluid,
    species=H,
)

my_model.advection_terms = [advection_term]

my_model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=inlet, value=1, species=H),
]
my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

concentration_field_fluid = F.VTXSpeciesExport(
    filename=f"H_fluid.bp", field=H, subdomain=fluid
)
concentration_field_tube = F.VTXSpeciesExport(
    filename=f"H_tube.bp", field=H, subdomain=tube
)

my_model.exports = [
    concentration_field_fluid,
    concentration_field_tube,
]


my_model.initialise()
set_log_level(LogLevel.INFO)
my_model.run()
