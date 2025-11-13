import gmsh

###############################################
###### CREATE FESTIM MESH FROM CAD MODEL ######
###############################################

# LOAD CAD AND INITIALIZE MESH

gmsh.initialize()
gmsh.option.setString(
    "Geometry.OCCTargetUnit", "M"
)  # make sure gmsh reads .step file in meters
gmsh.model.add("festim_mesh")

cad_file_path = "meshing/two_volumes.step"

entities = gmsh.model.occ.importShapes(cad_file_path)
gmsh.model.occ.synchronize()

# EXTRACT ALL VOLUMES
volumes = [e for e in gmsh.model.occ.getEntities() if e[0] == 3]

print(f"Extracted {len(volumes)} raw volumes from CAD.")

# volumes info for debugging
for v in volumes:
    com = gmsh.model.occ.getCenterOfMass(v[0], v[1])
    bbox = gmsh.model.getBoundingBox(v[0], v[1])

##### FRAGMENT VOLUMES & GENERATE SHARED SURFACES #####
print("Fragmenting volumes to define interfaces...")
gmsh.model.occ.fragment(volumes, [])
gmsh.model.occ.synchronize()

# FINAL VOLUMES AFTER FRAGMENT
final_volumes = gmsh.model.getEntities(dim=3)
print(f"Final number of volumes: {len(final_volumes)}")

##### TAG & NAME PHYSICAL GROUPS #####
# need to open mesh in gmsh gui to determine the proper tagging as below

# volumes
probe_marker = 1
breeder_marker = 2

gmsh.model.addPhysicalGroup(3, [1, 2, 3], probe_marker, name=f"probe")
gmsh.model.addPhysicalGroup(3, [4], breeder_marker, name=f"breeder")

# interface surfaces
surfaces = gmsh.model.getEntities(dim=2)
interface_surfaces = [1, 5, 6, 8, 9]
interface_marker = 3
gmsh.model.addPhysicalGroup(2, interface_surfaces, interface_marker, name="interfaces")

# other surfaces
inlet = 11
inlet_marker = 4

walls = 12
walls_marker = 5

outlet = 13
outlet_marker = 6

vacuum_surfaces = [4, 7, 10]
vacuum_marker = 7

gmsh.model.addPhysicalGroup(2, [inlet], inlet_marker, name="inlet")
gmsh.model.addPhysicalGroup(2, [outlet], outlet_marker, name="outlet")
gmsh.model.addPhysicalGroup(2, [walls], walls_marker, name="walls")
gmsh.model.addPhysicalGroup(2, vacuum_surfaces, vacuum_marker, name="vacuum")


##### MESH SIZE & REFINEMENT #####
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)

gmsh.model.occ.synchronize()

# refinement for tube with boundary layer refinement
inlet_outlet_walls = [
    inlet,
    outlet,
    walls,
]  # set distance field near inlet, outlet, and wall surfaces

distance_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", inlet_outlet_walls)

# use threshold field to refine near surfaces
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMin", 0.01
)  # smallest mesh size near surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "SizeMax", 0.2
)  # mesh size far from surfaces
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMin", 0.05
)  # distance where within which mesh is fully refined
gmsh.model.mesh.field.setNumber(
    threshold_field, "DistMax", 0.1
)  # distance where mesh transitions to coarse size

# set threshold field as background field (doesn't impact probe surface)
gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

##### SYNC & GENERATE MESH #####
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)

##### SAVE MESH #####
output_file = "meshing/festim_mesh.msh"
gmsh.write(output_file)
gmsh.finalize()
