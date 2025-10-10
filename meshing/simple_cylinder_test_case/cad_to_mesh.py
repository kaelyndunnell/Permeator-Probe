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

gmsh.model.add("cylinder")

cad_file_path = "meshing/simple_cylinder_test_case/simple_cylinder.step"
gmsh.model.occ.importShapes(cad_file_path)

gmsh.model.occ.synchronize()

##### TAG & NAME PHYSICAL GROUPS #####

# VOLUMES
volumes = gmsh.model.getEntities(
    dim=3
)  # gets volumes using 3. to get surfaces, you would use 2, etc.

# tagging this volume as tag 1
vol_marker = 1

# assign volumes with gmsh
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], vol_marker, name="fluid")

# SURFACES
surfaces = gmsh.model.occ.getEntities(dim=2)  # now getting surfaces using dimension = 2

# get surface ids from opening CAD (.step) file in gmsh api
inlet_tag = 2
outlet_tag = 3
wall_tag = 1

# markers for gmsh
inlet_marker = 1
outlet_marker = 2
wall_marker = 3

# assign surfaces with gmsh
gmsh.model.addPhysicalGroup(surfaces[0][0], [wall_tag], wall_marker, name="wall")
gmsh.model.addPhysicalGroup(surfaces[0][0], [outlet_tag], outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(surfaces[0][0], [inlet_tag], inlet_marker, name="inlet")

##### MESH SIZE & REFINEMENT #####

gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 30)

gmsh.model.occ.synchronize()

##### GENERATE MESH #####

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)  # 3D mesh

gmsh.fltk.run()  # comment out if want to run without GUI

gmsh.write("meshing/simple_cylinder_test_case/cylinder.msh")
gmsh.finalize()
