import gmsh
from mpi4py import MPI
from dolfinx import fem
import numpy as np
import festim as F
from basix.ufl import element
from dolfinx.io import VTXWriter, XDMFFile, gmsh as gmshio
import ufl
from dolfinx import cpp as _cpp
from dolfinx.log import set_log_level, LogLevel
import h_transport_materials as htm


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
    v_mag_cell = np.zeros(num_cells)
    dofs_per_cell = V0.dofmap.index_map_bs  # Should be 1 for DG0

    # Use dolfinx function to integrate |u|^2 over each cell
    u_sq_expr = ufl.dot(u, u)
    u_sq_fn = fem.Function(V0)
    u_sq_fn.interpolate(fem.Expression(u_sq_expr, V0.element.interpolation_points))
    v_mag_cell[:] = np.sqrt(u_sq_fn.x.array[:])  # now per cell

    # Stabilization term per cell
    D_art_array = delta * v_mag_cell * h_as_function.x.array

    # Return as DG0 function
    D_art_fn = fem.Function(V0)
    D_art_fn.x.array[:] = D_art_array

    return D_art_fn


gmsh.initialize()
gmsh.model.add("mwe")

r_inner = 0.1
r_tube = 0.11
length = 0.4

factory = gmsh.model.occ

fluid = factory.addCylinder(0, 0, 0, length, 0, 0, r_inner)
tube = factory.addCylinder(0, 0, 0, length, 0, 0, r_tube)
walls, interface = factory.cut(
    [(3, tube)], [(3, fluid)], removeObject=True, removeTool=False
)

factory.synchronize()
interface_surfaces = [1]
interface_tag = gmsh.model.addPhysicalGroup(
    2, interface_surfaces
)  # these are tagged wrong
gmsh.model.setPhysicalName(2, interface_tag, "interface")

outlet = gmsh.model.addPhysicalGroup(2, [2])
gmsh.model.setPhysicalName(2, outlet, "outlet")

inlet = gmsh.model.addPhysicalGroup(2, [3])
gmsh.model.setPhysicalName(2, inlet, "inlet")

wall_surf = gmsh.model.addPhysicalGroup(2, [4, 5, 6])
gmsh.model.setPhysicalName(2, wall_surf, "wall_surf")

fluid_tag = gmsh.model.addPhysicalGroup(3, [fluid])
gmsh.model.setPhysicalName(3, fluid_tag, "fluid")

wall_tag = gmsh.model.addPhysicalGroup(3, [walls[0][1]])
gmsh.model.setPhysicalName(3, wall_tag, "wall")

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
gmsh.model.mesh.generate(3)
gmsh.write("mwe.msh")

mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)
my_mesh = mesh_data.mesh
gmsh.finalize()

my_model = F.HydrogenTransportProblemDiscontinuous()
my_model.mesh = F.Mesh(my_mesh)
my_model.facet_meshtags = mesh_data.facet_tags
my_model.volume_meshtags = mesh_data.cell_tags

el = element("Lagrange", my_mesh.topology.cell_name(), 1, shape=(my_mesh.geometry.dim,))
V = fem.functionspace(my_mesh, el)

velocity = fem.Function(V)
fluid_cells = mesh_data.cell_tags.find(fluid_tag)
velocity.interpolate(
    lambda x: (
        -5e-2
        * (1.0 - (x[0] / length) ** 2)
        * (1.0 - (x[1] ** 2 + x[2] ** 2) / r_inner**2),
        np.zeros_like(x[0]),
        np.zeros_like(x[0]),
    ),
    cells0=fluid_cells,
)

D_0_PbLi = 4.03e-08  # m2/s
E_D_PbLi = 0.2021  # eV

D_diff = D_0_PbLi * ufl.exp(-E_D_PbLi / (F.k_B * 500))

# add stabilization term for diffusion
D_art = evaluate_stabalisation_term(mesh=my_mesh, u=velocity, delta=0.1)
D_expr = D_diff + D_art
V = fem.functionspace(my_mesh, ("CG", 1))
D_fluid = fem.Function(V)
D_fluid.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

# my_writer_2 = VTXWriter(MPI.COMM_WORLD, "D_field.bp", D_pbli, "BP5")
# my_writer_2.write(t=0)

breeder_material = F.Material(
    D=D_fluid, K_S_0=1.43e23, E_K_S=0.13
)  # https://theses.hal.science/tel-04906459v1

# probe material parameters -- alpha-Fe
htm_D_Fe = htm.diffusivities.filter(material="Fe").mean()
htm_S_Fe = htm.solubilities.filter(material="Fe").mean()

iron = F.Material(
    D_0=htm_D_Fe.pre_exp.magnitude,
    E_D=htm_D_Fe.act_energy.magnitude,
    K_S_0=htm_S_Fe.pre_exp.magnitude,
    E_K_S=htm_S_Fe.act_energy.magnitude,
)

# D_diff = 1e-4 * ufl.exp(-0.2 / (F.k_B * 500))
# D_art = evaluate_stabalisation_term(mesh=my_mesh, u=velocity, delta=0.1)
# D_expr = D_diff + D_art
# V_CG = fem.functionspace(my_mesh, ("CG", 1))
# D_fluid = fem.Function(V_CG)
# D_fluid.interpolate(fem.Expression(D_expr, V_CG.element.interpolation_points))

my_writer = VTXWriter(MPI.COMM_WORLD, "velocity_field.bp", velocity, "BP5")
my_writer.write(t=0)
my_writer_2 = VTXWriter(MPI.COMM_WORLD, "D_field.bp", D_fluid, "BP5")
my_writer_2.write(t=0)

# with XDMFFile(MPI.COMM_WORLD, "facet_tags.xdmf", "w") as xdmf:
#     xdmf.write_mesh(my_mesh)
#     xdmf.write_meshtags(mesh_data.facet_tags, my_mesh.geometry)
# with XDMFFile(MPI.COMM_WORLD, "volume_tags.xdmf", "w") as xdmf:
#     xdmf.write_mesh(my_mesh)
#     xdmf.write_meshtags(mesh_data.cell_tags, my_mesh.geometry)


# dummy_fluid = F.Material(D=D_fluid, K_S_0=1, E_K_S=0)
dummy_tube = F.Material(D_0=3.87e-8, E_D=0.04, K_S_0=3.07e23, E_K_S=0.279)

inlet = F.SurfaceSubdomain(id=inlet)
outlet = F.SurfaceSubdomain(id=outlet)
vacuum = F.SurfaceSubdomain(id=wall_surf)

fluid_sd = F.VolumeSubdomain(
    id=fluid_tag,
    material=breeder_material,
)
tube_sd = F.VolumeSubdomain(
    id=wall_tag,
    material=dummy_tube,
)

my_model.subdomains = [inlet, outlet, vacuum, fluid_sd, tube_sd]
my_model.surface_to_volume = {inlet: fluid_sd, outlet: fluid_sd, vacuum: tube_sd}

H = F.Species("H", subdomains=my_model.volume_subdomains)
my_model.species = [H]
my_model.interfaces = [
    F.Interface(id=interface_tag, subdomains=[fluid_sd, tube_sd], penalty_term=10000)
]
my_model.temperature = 500

advection_term = F.AdvectionTerm(velocity=velocity, subdomain=fluid_sd, species=H)
my_model.advection_terms = [advection_term]

my_model.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=inlet, value=1, species=H),
    F.FixedConcentrationBC(subdomain=vacuum, value=0, species=H),
]

my_model.settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)

concentration_field_fluid = F.VTXSpeciesExport(
    filename="H_fluid.bp", field=H, subdomain=fluid_sd
)
concentration_field_tube = F.VTXSpeciesExport(
    filename="H_tube.bp", field=H, subdomain=tube_sd
)
my_model.exports = [concentration_field_fluid, concentration_field_tube]

my_model.initialise()
set_log_level(LogLevel.INFO)
my_model.run()
