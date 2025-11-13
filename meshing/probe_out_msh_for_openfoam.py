import gmsh
import os
from dolfinx import plot
import pyvista
import numpy as np
from mpi4py import MPI

########################################
###### CREATE MESH FROM CAD MODEL ######
########################################

# LOAD CAD AND INITIALIZE MESH

gmsh.initialize()

gmsh.option.setString(
    "Geometry.OCCTargetUnit", "M"
)  # make sure gmsh reads .step file in meters

gmsh.model.add("probe_out")

cad_file_path = "meshing/probe_out.step"
gmsh.model.occ.importShapes(cad_file_path)

gmsh.model.occ.synchronize()


##### TAG & NAME PHYSICAL GROUPS #####

# VOLUMES
volumes = gmsh.model.getEntities(
    dim=3
)  # gets volumes using 3. to get surfaces, you would use 2, etc.

vol_marker = 1

gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], vol_marker)
gmsh.model.setPhysicalName(volumes[0][0], vol_marker, "breeder")


# SURFACES
surfaces = gmsh.model.occ.getEntities(dim=2)  # now getting surfaces using dimension = 2

inlet_tag = 5
outlet_tag = 4
wall_tags = [1, 2, 3]
probe_tags = [6, 7, 8, 9, 10]

# markers for gmsh
inlet_marker = 1
outlet_marker = 2
wall_marker = 3
probe_marker = 4

# assign surfaces with gmsh
gmsh.model.addPhysicalGroup(surfaces[0][0], wall_tags, wall_marker, name="wall")
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

gmsh.write("meshing/probe_out_openfoam.msh")
gmsh.finalize()
