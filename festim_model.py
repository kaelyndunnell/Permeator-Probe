import festim as F  # using festim2
import numpy as np
from scifem import assemble_scalar
from dolfinx import fem
import ufl
from dolfinx import cpp as _cpp
from openfoam_to_festim import read_openfoam_data
from dolfinx.log import set_log_level, LogLevel
import matplotlib.pyplot as plt


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
    festim_final_time,
    steady=True,
):

    # markers for gmsh TODO: do not make this repetitive
    inlet_marker = 1
    outlet_marker = 2
    wall_marker = 3
    probe_marker = 4

    # READ MESH

    p, u, mesh, nut, facet_meshtags, volume_meshtags = read_openfoam_data(
        openfoam_data_file, final_time=openfoam_final_time
    )

    # DEFINE & INITIALIZE MODEL

    print("Building FESTIM model...")

    my_model = F.HydrogenTransportProblem()

    my_model.mesh = F.Mesh(mesh)
    my_model.facet_meshtags = facet_meshtags
    my_model.volume_meshtags = volume_meshtags

    D_0_PbLi = 4.03e-08  # m2/s
    E_D_PbLi = 0.2021  # eV

    D_diff = D_0_PbLi * ufl.exp(-E_D_PbLi / (F.k_B * breeder_temperature))

    # add stabilization term for diffusion
    D_art = evaluate_stabalisation_term(mesh=mesh, u=u, delta=delta)

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

    N_A = 6.022e23
    c_in = (
        1.5e-2 * N_A
    )  # atms/m3 , inspired by tritium concentration (mol/m3) of OB loop from Utili 2023

    my_model.boundary_conditions = [
        F.FixedConcentrationBC(subdomain=inlet, value=c_in, species=H),
        F.FixedConcentrationBC(subdomain=probe, value=0, species=H),
    ]

    advection = F.AdvectionTerm(velocity=u, subdomain=vol, species=H)
    my_model.advection_terms = [advection]

    # SETTINGS

    if not steady:
        dt = F.Stepsize(
            initial_value=1,
            growth_factor=1.05,
            cutback_factor=0.9,
            target_nb_iterations=5,
        )
        my_model.settings = F.Settings(
            atol=1e04,
            rtol=1e-10,
            transient=True,
            final_time=festim_final_time,
            stepsize=dt,
        )

    else:
        my_model.settings = F.Settings(
            atol=1e04,
            rtol=1e-10,
            transient=False,
        )

    # EXPORTS

    outlet_surface_flux = SurfaceAdvectionFlux(
        field=H,
        surface=outlet,
        filename=f"{results_folder}/outlet_surface_flux.csv",
        velocity_field=u,
    )
    probe_flux = F.SurfaceFlux(
        field=H, surface=probe, filename=f"{results_folder}/probe_surface_flux.csv"
    )

    inventory = F.TotalVolume(
        field=H, volume=vol, filename=f"{results_folder}/inventory.csv"
    )

    concentration_field = F.VTXSpeciesExport(filename=f"{results_folder}/H.bp", field=H)

    my_model.exports = [
        outlet_surface_flux,
        probe_flux,
        inventory,
        concentration_field,
    ]

    return my_model


if __name__ == "__main__":

    my_model = build_festim_model(
        openfoam_data_file="OpenFOAM/turbulent-case/probe.foam",
        openfoam_final_time=300,
        breeder_temperature=603.15,
        delta=0.1,
        results_folder="festim_results",
        festim_final_time=300,  # ignored for steady state
        steady=False,
    )

    # INITIALISE AND RUN
    my_model.initialise()
    my_model.run()
