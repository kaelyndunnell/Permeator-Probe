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

> **_NOTE:_**  The `environment.yml` file does not include OpenFOAM installation, which will need to be done independently. Refer to [OpenFOAMv13 installation instructions](https://openfoam.org/download/) for your operating system. 


## Workflow: 

To run the probe-in OpenFOAM simulation, first create the mesh: 

 ```
 python meshing/cad_to_gmsh_for_openfoam.py
 ```

 Convert mesh to OpenFOAM format: 

 ```
 mv meshing/probe_in_openfoam_mesh.msh OpenFOAM/laminar-case
 cd OpenFOAM/laminar-case
 gmshToFoam probe_in_openfoam_mesh.msh
 ```

 > **_NOTE:_**  Ensure the bounding box of the mesh is in the proper units -- with our geometry, it should be `(-0.067 -0.1 -0.067) (0.067 0.1 0.067)`. If it isn't, use the `transformPoints` command accordingly.

 Check the mesh: 

 ```
 checkMesh
 ```
 
 Finally, run the OpenFOAM simulation:

 ```
 foamRun -solver incompressibleFluid
 ```

When the OpenFOAM simulation is finished, the results can be fed into the FESTIM simulation by running: 

```
python two_volume_festim_model.py
```

The results can be read from the `festim_results` folder that is created from the FESTIM simulation. 