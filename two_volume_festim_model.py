import festim as F  # using festim2
import numpy as np
from scifem import assemble_scalar
from dolfinx import fem
from dolfinx.io import VTXWriter, XDMFFile
import ufl
from dolfinx import cpp as _cpp
from openfoam_to_festim import read_openfoam_data
from dolfinx.log import set_log_level, LogLevel
import matplotlib.pyplot as plt
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
from basix.ufl import element
import h_transport_materials as htm


class SurfaceAdvectionFlux(F.SurfaceFlux):
    """Computes the advection flux of a field on a given surface

    Args:
        field (festim.Species): species for which the surface flux is computed
        surface (festim.SurfaceSubdomain1D): surface subdomain
        filename (str, optional): name of the file to which the surface flux is exported

    Attributes:
        see `festim.SurfaceFlux`
    """

    def __init__(self, field, surface, filename, velocity_field):

        super().__init__(field=field, surface=surface, filename=filename)
        self.velocity_field = velocity_field

    @property
    def title(self):
        return f"{self.field.name} advection flux surface {self.surface.id}"

    def compute(self, u, ds: ufl.Measure, entity_maps=None):
        """Computes the value of the advection flux at the surface

        Args:
            u: field for which the flux is computed
            velocity (float): magnitude of velocity in m/s
            ds: surface measure of the model
            entity_maps: entity maps relating parent mesh and submesh
        """

        # obtain mesh normal from field
        # if case multispecies, solution is an index, use sub_function_space
        if isinstance(u, ufl.indexed.Indexed):
            mesh = self.field.sub_function_space.mesh
        else:
            mesh = u.function_space.mesh
        n = ufl.FacetNormal(mesh)

        surface_flux = assemble_scalar(
            fem.form(
                -self.D * ufl.dot(ufl.grad(u), n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )
        advective_flux = assemble_scalar(
            fem.form(
                u * ufl.inner(self.velocity_field, n) * ds(self.surface.id),
                entity_maps=entity_maps,
            )
        )

        self.value = surface_flux + advective_flux
        self.data.append(self.value)


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


def build_festim_model(
    openfoam_data_file,
    openfoam_final_time,
    breeder_temperature,
    delta,
    results_folder,
):

    # READ OPENFOAM MESH

    p, u, openfoam_mesh, nut, facet_meshtags, volume_meshtags = read_openfoam_data(
        openfoam_data_file, final_time=openfoam_final_time
    )

    # READ GMSH MESH
    model_rank = 0
    festim_mesh_data = gmshio.read_from_msh(
        "meshing/festim_mesh.msh", MPI.COMM_WORLD, model_rank, gdim=3
    )
    festim_mesh = festim_mesh_data.mesh

    # DEFINE & INITIALIZE MODEL

    print("Building FESTIM model...")

    my_model = F.HydrogenTransportProblemDiscontinuous()

    my_model.mesh = F.Mesh(festim_mesh)
    my_model.facet_meshtags = festim_mesh_data.facet_tags
    my_model.volume_meshtags = festim_mesh_data.cell_tags

    # interpolate OpenFOAM velocity field onto FESTIM mesh
    el = element(
        "Lagrange",
        festim_mesh.topology.cell_name(),
        1,
        shape=(festim_mesh.geometry.dim,),
    )
    V_openfoam = fem.functionspace(openfoam_mesh, el)
    V_festim = fem.functionspace(festim_mesh, el)

    u_openfoam = fem.Function(V_openfoam)
    u_openfoam.interpolate(
        u
    )  # u is a fem.function.Function! paraview visualization of u_openfoam is accurate with this formulation
    festim_velocity = fem.Function(V_festim)

    festim_cells = festim_mesh_data.cell_tags.find(2)  # breeder cells to interpolate to

    interpolation_data = fem.create_interpolation_data(
        V_to=V_festim, V_from=V_openfoam, cells=festim_cells
    )

    festim_velocity.interpolate_nonmatching(
        u_openfoam, cells=festim_cells, interpolation_data=interpolation_data
    )

    # my_writer = VTXWriter(MPI.COMM_WORLD, "velocity_field.bp", festim_velocity, "BP5")
    # my_writer.write(t=0)

    D_0_PbLi = 4.03e-08  # m2/s
    E_D_PbLi = 0.2021  # eV

    D_diff = D_0_PbLi * ufl.exp(-E_D_PbLi / (F.k_B * breeder_temperature))

    # add stabilization term for diffusion
    D_art = evaluate_stabalisation_term(
        mesh=festim_mesh, u=festim_velocity, delta=delta
    )
    D_expr = D_diff + D_art
    V = fem.functionspace(festim_mesh, ("CG", 1))
    D_pbli = fem.Function(V)
    D_pbli.interpolate(fem.Expression(D_expr, V.element.interpolation_points))

    # my_writer_2 = VTXWriter(MPI.COMM_WORLD, "D_field.bp", D_pbli, "BP5")
    # my_writer_2.write(t=0)

    breeder_material = F.Material(
        D=D_pbli, K_S_0=1.43e23, E_K_S=0.13
    )  # https://theses.hal.science/tel-04906459v1

    u = htm.ureg
    solubility_alpha_iron = htm.Solubility(
        S_0=0.51 * u.mol * u.m**-3 * u.Pa**-0.5,
        E_S=27 * u.kJ * u.mol**-1,
        source="10.1016/S0022-3115(96)00670-8",
    )

    # probe material parameters -- alpha-Fe
    alpha_iron = F.Material(
        D_0=3.87e-8,
        E_D=0.04,
        K_S_0=solubility_alpha_iron.pre_exp.magnitude,
        E_K_S=solubility_alpha_iron.act_energy.magnitude,
    )  # from doi: 10.1016/S0022-3115(98)00038-5

    # SET DOMAINS

    # use same tags as gmsh markers
    probe_marker = 1
    breeder_marker = 2

    probe = F.VolumeSubdomain(id=probe_marker, material=alpha_iron)
    breeder = F.VolumeSubdomain(id=breeder_marker, material=breeder_material)

    interface_marker = 3
    inlet_marker = 4
    outlet_marker = 6
    wall_marker = 5
    vacuum_marker = 7

    inlet = F.SurfaceSubdomain(id=inlet_marker)
    outlet = F.SurfaceSubdomain(id=outlet_marker)
    wall = F.SurfaceSubdomain(id=wall_marker)
    vacuum = F.SurfaceSubdomain(id=vacuum_marker)

    my_model.subdomains = [
        probe,
        breeder,
        inlet,
        outlet,
        wall,
        vacuum,
    ]

    my_model.surface_to_volume = (
        {  # anything defines BC for needs to be here (and exports)
            inlet: breeder,
            outlet: breeder,
            vacuum: probe,
        }
    )

    H = F.Species("H", subdomains=my_model.volume_subdomains)
    my_model.species = [H]

    my_model.interfaces = [
        F.Interface(
            id=interface_marker,
            subdomains=[breeder, probe],
            penalty_term=1e25,
        ),
    ]

    # SET TEMP AND BOUNDARY CONDITIONS

    my_model.temperature = breeder_temperature  # K

    advection = F.AdvectionTerm(velocity=festim_velocity, subdomain=breeder, species=H)
    my_model.advection_terms = [advection]

    N_A = 6.022e23
    c_in = (
        1.5e-2 * N_A
    )  # atms/m3 , inspired by tritium concentration (mol/m3) of OB loop from Utili 2023

    print(f"Inlet Concentration is {c_in} #/m3.")

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, value=c_in, species=H),
        F.FixedConcentrationBC(subdomain=vacuum, value=0, species=H),
    ]

    # SETTINGS

    my_model.settings = F.Settings(
        atol=1e-10,
        rtol=1e-10,
        transient=False,
    )

    # EXPORTS

    # outlet_advective_flux = SurfaceAdvectionFlux(
    #     field=H,
    #     surface=outlet,
    #     filename=f"{results_folder}/outlet_advective_flux.csv",
    #     velocity_field=festim_velocity,
    # )
    permeation_flux = F.SurfaceFlux(
        field=H, surface=vacuum, filename=f"{results_folder}/permeation_flux.csv"
    )

    concentration_field_breeder = F.VTXSpeciesExport(
        filename=f"{results_folder}/H_breeder.bp", field=H, subdomain=breeder
    )
    concentration_field_probe = F.VTXSpeciesExport(
        filename=f"{results_folder}/H_probe.bp", field=H, subdomain=probe
    )

    my_model.exports = [
        # outlet_advective_flux,
        permeation_flux,
        concentration_field_breeder,
        concentration_field_probe,
    ]

    return my_model


if __name__ == "__main__":

    my_model = build_festim_model(
        openfoam_data_file="OpenFOAM/laminar-case/probe.foam",
        openfoam_final_time=300,
        breeder_temperature=603.15,
        delta=0.1,
        results_folder="festim_results",
    )

    # INITIALISE AND RUN
    my_model.initialise()
    set_log_level(LogLevel.INFO)
    my_model.run()
