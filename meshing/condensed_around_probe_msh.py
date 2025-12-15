import gmsh
import os
from dolfinx import plot
import pyvista
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmsh as gmshio

########################################
###### CREATE MESH FROM CAD MODEL ######
########################################

# LOAD CAD AND INITIALIZE MESH

gmsh.initialize()

gmsh.option.setString(
    "Geometry.OCCTargetUnit", "M"
)  # make sure gmsh reads .step file in meters

gmsh.model.add("candido_probe")

cad_file_path = "meshing/condensed_around_probe.step"
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
inlet_tag = 3
outlet_tag = 1
wall_tag = 2
probe_tags = [4, 5, 6, 7]
# 4 and 5 are the caps (patches at ends) of the inner capsule (which is NOT included in the surface of the capsule (tag 6)!)

# markers for gmsh
inlet_marker = 1
outlet_marker = 2
wall_marker = 3
probe_marker = 4

# assign surfaces with gmsh
gmsh.model.addPhysicalGroup(surfaces[0][0], [wall_tag], wall_marker, name="wall")
gmsh.model.addPhysicalGroup(surfaces[0][0], [outlet_tag], outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(surfaces[0][0], [inlet_tag], inlet_marker, name="inlet")
gmsh.model.addPhysicalGroup(surfaces[0][0], probe_tags, probe_marker, name="probe")


##### MESH SIZE & REFINEMENT #####

gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)  # high refinement for probe

gmsh.model.occ.synchronize()

# refinement for tube
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
    threshold_field, "SizeMin", 0.003
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

gmsh.model.occ.synchronize()

##### GENERATE MESH #####

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)  # 3D mesh

gmsh.fltk.run()  # comment out if want to run without GUI

gmsh.write("meshing/condensed.msh")
gmsh.finalize()
