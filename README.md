# Permeator Probe for Tritium Measurement

This repository develops a helical permeator probe design for simulating tritium measurement in fusion fuel cycle components. The permeator probe is based on the geometry from [HyPer-QuarCh II: A laboratory-scale device for hydrogen isotopes permeation experiments](https://doi.org/10.1016/j.fusengdes.2021.112920) by L. Candido et al (2021). Using OpenFOAM for CFD simulations as coupled to FESTIM for tritium transport, the permeator probe can be modeled and tested as a potential tritium measurement device.  


## How to Run: 

Clone the repository: 

 ```
git clone https://github.com/kaelyndunnell/Permeator-Probe
cd Permeator-Probe 
```

Run this command to create a new environment with the right dependencies (e.g. dolfinx, OpenFOAM):

```
conda env create -f environment.yml 
```

Then, activate the environment: 

```
conda activate permeator-probe-env
```

## Workflow: 

Create mesh: 

 ```
 python meshing/cad_to_gmsh.py
 ```

 Convert mesh to OpenFOAM format: 

 ```
 mv meshing/probe_breeder.msh OpenFOAM/probe-case
 cd OpenFOAM/probe-case
 gmshToFoam probe_breeder.msh
 ```

 Check the mesh with OpenFOAM: 

 ```
 checkMesh
 ```

 > **_NOTE:_**  Ensure the bounding box of the mesh is in proper units (i.e., our mesh is in centimeters, so the bounding box should be on the order of 1e-2). If it isn't, run 
 `transformPoints "scale=(0.001 0.001 0.001)"` in the command line.

 Finally, run the OpenFOAM simulation:

 ```
 foamRun -solver incompressibleFluid
 ```

