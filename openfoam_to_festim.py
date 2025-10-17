import festim as F
from foam2dolfinx import OpenFOAMReader
from dolfinx import fem
from festim.helpers import nmm_interpolate
from dolfinx.io import VTXWriter
from mpi4py import MPI


def read_openfoam_data(file_name, final_time):
    """
    Read OpenFOAM data from a file and return the pressure, velocity, and viscosity (if it exists) fields.
    """
    print("Reading OpenFOAM data...")
    openfoam_reader = OpenFOAMReader(filename=file_name, cell_type=10)
    p = openfoam_reader.create_dolfinx_function(t=final_time, name="p")
    u = openfoam_reader.create_dolfinx_function(t=final_time, name="U")
    mesh = openfoam_reader.dolfinx_meshes_dict["default"]

    try:
        # read turbulent viscosity if it exists
        nut = openfoam_reader.create_dolfinx_function(t=final_time, name="nut")
        print("'nut' field found and read.")
    except Exception:
        nut = None
        print("no 'nut' field found.")

    return p, u, mesh, nut


def export_openfoam_data(p, u):
    """
    Export OpenFOAM data to VTX files.
    """

    print("Exporting OpenFOAM data")
    writer_p = VTXWriter(
        MPI.COMM_WORLD,
        "OpenFOAM/pressure.bp",
        p,
        "BP5",
    )
    writer_u = VTXWriter(
        MPI.COMM_WORLD,
        "OpenFOAM/velocity.bp",
        u,
        "BP5",
    )

    writer_p.write(t=0)
    writer_u.write(t=0)


if __name__ == "__main__":
    # read openfoam data
    p, u = read_openfoam_data("OpenFOAM/probe-case/case.foam", final_time=100)

    export_openfoam_data(p, u)
