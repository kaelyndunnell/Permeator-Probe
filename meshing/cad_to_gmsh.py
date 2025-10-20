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

gmsh.option.setString(
    "Geometry.OCCTargetUnit", "M"
)  # make sure gmsh reads .step file in meters

gmsh.model.add("candido_probe")

cad_file_path = "meshing/breeder_with_probe_surface.step"
gmsh.model.occ.importShapes(cad_file_path)

gmsh.model.occ.synchronize()


##### TAG & NAME PHYSICAL GROUPS #####

# VOLUMES
volumes = gmsh.model.getEntities(
    dim=3
)  # gets volumes using 3. to get surfaces, you would use 2, etc.

# exported breeder and probe volume, tagged according to visibility in gmsh
# treating tag and marker as the same here
vacuum_vol_tag = 1
probe_vol_tag = 2
breeder_vol_tag = 3

# assign volumes with gmsh
# gmsh.model.addPhysicalGroup(dimension, [cad tags to include in this group], gmsh tag, name of group)
gmsh.model.addPhysicalGroup(
    volumes[0][0], [vacuum_vol_tag], vacuum_vol_tag, name="vacuum"
)
gmsh.model.addPhysicalGroup(volumes[0][0], [probe_vol_tag], probe_vol_tag, name="probe")
gmsh.model.addPhysicalGroup(
    volumes[0][0], [breeder_vol_tag], breeder_vol_tag, name="breeder"
)


# SURFACES
surfaces = gmsh.model.occ.getEntities(dim=2)  # now getting surfaces using dimension = 2

# get surface ids from opening CAD (.step) file in gmsh api
# choose tools > visibility and select the surface to see
inlet_tag = 8
outlet_tag = 9
wall_tag = 10
inside_probe_surf = [4, 11, 12]
outside_probe_surf = [5, 6, 7, 13, 14, 15]
vacuum_surf_tags = [1, 2, 3]

# 5 and 13 are the same surface, 6 and 7 are the outer loops of 14 and 15

# markers for gmsh
inlet_marker = 1
outlet_marker = 2
wall_marker = 3
probe_marker = 4
vacuum_marker = 5

# assign surfaces with gmsh
gmsh.model.addPhysicalGroup(surfaces[0][0], [wall_tag], wall_marker, name="wall")
gmsh.model.addPhysicalGroup(surfaces[0][0], [outlet_tag], outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(surfaces[0][0], [inlet_tag], inlet_marker, name="inlet")
gmsh.model.addPhysicalGroup(
    surfaces[0][0], outside_probe_surf, probe_marker, name="probe"
)
gmsh.model.addPhysicalGroup(
    surfaces[0][0], inside_probe_surf + vacuum_surf_tags, vacuum_marker, name="vacuum"
)


##### MESH SIZE & REFINEMENT #####

# probe_minor_radius = 0.1  # (d/2) in cm, from CAD
# tube_radius = 13  # in cm, from CAD

gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)  # high refinement for probe

gmsh.model.occ.synchronize()

# refinement for tube with boundary layer refinement
inlet_outlet_wall = [
    inlet_marker,
    outlet_marker,
    wall_marker,
]  # set distance field near inlet, outlet, and wall surfaces

distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", inlet_outlet_wall)

# use threshold field to refine near surfaces
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMin", 0.005
)  # smallest mesh size near surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMax", 0.10
)  # mesh size far from surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMin", 0.015
)  # distance where within which mesh is fully refined
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMax", 0.02
)  # distance where mesh transitions to coarse size

# set threshold field as background field (doesn't impact probe surface)
gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)


##### GENERATE MESH #####

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)  # 3D mesh

gmsh.fltk.run()  # comment out if want to run without GUI

gmsh.write("meshing/probe_breeder_with_boundary.msh")
gmsh.finalize()

########################################
############ VISUALIZE MESH ############
########################################

# # read mesh, volume tags, and surface tags
# mesh, volume_tags, surface_tags = gmshio.read_from_msh(
#     "meshing/probe_breeder_with_boundary.msh", MPI.COMM_WORLD, 0, gdim=3
# )

# print(f"Volume tags: {np.unique(volume_tags.values)}")
# print(f"Surface tags: {np.unique(surface_tags.values)}")

# # initialize pyvista for visualization
# pyvista.start_xvfb()
# pyvista.set_jupyter_backend("html")

# tdim = mesh.topology.dim  # gets dimensions of the mesh cells / elements

# mesh.topology.create_connectivity(
#     tdim, tdim
# )  # generates and stores connectivity between mesh entities of all dimensions (tdim)

# topology, cell_types, geometry = plot.vtk_mesh(
#     mesh, tdim
# )  # preps mesh for plotting in pyvista
# grid = pyvista.UnstructuredGrid(
#     topology, cell_types, geometry
# )  # makes the mesh in pyvista

# plotter = pyvista.Plotter()  # initialize plotter, then plot below
# plotter.add_mesh(
#     grid, show_edges=True, opacity=0.5
# )  # change show_edges to False if want to visualize just surfaces


# if not pyvista.OFF_SCREEN:
#     plotter.show()
# else:
#     figure = plotter.screenshot("mesh.png")
