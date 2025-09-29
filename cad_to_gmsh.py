import gmsh
import os
from dolfinx import plot
import pyvista
import numpy as np
from dolfinx.io import gmshio
from mpi4py import MPI

########################################
###### CREATE MESH FROM CAD MODEL ######
########################################

# LOAD CAD AND INITIALIZE MESH

gmsh.initialize()

gmsh.model.add("candido_probe")

cad_file_path = "breeder.step"
gmsh.model.occ.importShapes(cad_file_path)

gmsh.model.occ.synchronize()

##### TAG & NAME PHYSICAL GROUPS #####

# VOLUMES
volumes = gmsh.model.getEntities(
    dim=3
)  # gets volumes using 3. to get surfaces, you would use 2, etc.

# only exporting the breeder volume for the simulation! so have one volume
# tagging this volume as tag 1
vol_marker = 1

# assign volumes with gmsh
# gmsh.model.addPhysicalGroup(dimension, [cad tags to include in this group], gmsh tag)
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], vol_marker)

# set name of volume groups
# gmsh.model.setPhysicalName(dimension, gmsh tag, group name)
gmsh.model.setPhysicalName(volumes[0][0], vol_marker, "breeder")


# SURFACES
surfaces = gmsh.model.occ.getEntities(dim=2)  # now getting surfaces using dimension = 2

# get surface ids from opening CAD (.step) file in gmsh api
# choose tools > visibility and select the surface to see
wall_tag = 1
outlet_tag = 2
inlet_tag = 3
probe_tags = [
    4,
    5,
    # 6,
]  # 4 and 5 are the caps (patches at ends) of the inner capsule (which is NOT included in the surface of the capsule (tag 6)!)

# markers for gmsh
wall_marker = 1
outlet_marker = 2
inlet_marker = 3
probe_marker = 4

# assign surfaces with gmsh
gmsh.model.addPhysicalGroup(surfaces[0][0], [wall_tag], wall_marker)
gmsh.model.setPhysicalName(surfaces[0][0], wall_marker, "wall")

gmsh.model.addPhysicalGroup(surfaces[0][0], [outlet_tag], outlet_marker)
gmsh.model.setPhysicalName(surfaces[0][0], outlet_marker, "outlet")

gmsh.model.addPhysicalGroup(surfaces[0][0], [inlet_tag], inlet_marker)
gmsh.model.setPhysicalName(surfaces[0][0], inlet_marker, "inlet")

gmsh.model.addPhysicalGroup(surfaces[0][0], probe_tags, probe_marker)
gmsh.model.setPhysicalName(surfaces[0][0], probe_marker, "probe")


##### MESH SIZE & REFINEMENT #####

# refine mesh near the probe
distance = gmsh.model.mesh.field.add(
    "Distance"
)  # create mesh size field calculating distance from specified model entities to other parts of geometry
# distance variable stores the id of this field
gmsh.model.mesh.field.setNumbers(
    distance, "FacesList", probe_tags
)  # calculating distance from probe surfaces
# returns calculated distance for each mesh point in model to specified entities (here, probe surfaces)

probe_minor_radius = 0.1  # (d/2) in cm, from CAD
resolution = probe_minor_radius / 2

threshold = gmsh.model.mesh.field.add(
    "Threshold"
)  # creates new mesh size field to control mesh element sizes
gmsh.model.mesh.field.setNumber(
    threshold, "IField", distance
)  # bases threshold input field off of field id of distance field (stored in distance variable)
gmsh.model.mesh.field.setNumber(
    threshold, "LcMin", resolution
)  # minimum characteristic length
gmsh.model.mesh.field.setNumber(
    threshold, "LcMax", 20 * resolution
)  # maximum characteristic length
gmsh.model.mesh.field.setNumber(
    threshold, "DistMin", 0.5 * probe_minor_radius
)  # minimum characteristic distance
gmsh.model.mesh.field.setNumber(
    threshold, "DistMax", probe_minor_radius
)  # maximum characteristic distance

# set minimum mesh field size via probe refinement
minimum = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold])
gmsh.model.mesh.field.setAsBackgroundMesh(
    minimum
)  # sets this as default background mesh size

# # to only show surfaces (not volumes) in pyvista
# gmsh.option.setNumber("Mesh.SurfaceFaces", 1)  # show 2D faces
# gmsh.option.setNumber("Mesh.VolumeEdges", 0)  # hide 3D element edges

##### GENERATE MESH #####

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)  # 3D mesh
gmsh.write("probe_breeder.msh")
gmsh.finalize()

########################################
############ VISUALIZE MESH ############
########################################

# read mesh, volume tags, and surface tags
mesh, volume_tags, surface_tags = gmshio.read_from_msh(
    "probe_breeder.msh", MPI.COMM_WORLD, 0, gdim=3
)

print(f"Volume tags: {np.unique(volume_tags.values)}")
print(f"Surface tags: {np.unique(surface_tags.values)}")

# initialize pyvista for visualization
pyvista.start_xvfb()
pyvista.set_jupyter_backend("html")

# TODO comment rest of this

tdim = mesh.topology.dim

mesh.topology.create_connectivity(tdim, tdim)
topology, cell_types, geometry = plot.vtk_mesh(mesh, tdim)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

plotter = pyvista.Plotter()
plotter.add_mesh(
    grid, show_edges=True, opacity=0.5
)  # change show_edges to False if want to visualize just surfaces


if not pyvista.OFF_SCREEN:
    plotter.show()
else:
    figure = plotter.screenshot("mesh.png")
