#!/bin/bash

for i in $(seq 1 4); do # flow rates, same as in parametrization model

    cd OpenFOAM/flow_rate_parametrization/probe_case_${i}_kEpsilon

    gmshToFoam openfoam_mesh.msh
    checkMesh
    simpleFoam

    cd .. 
    gmshToFoam openfoam_mesh.msh
    checkMesh
    mv OpenFOAM/flow_rate_parametrization/probe_case_${i}_kEpsilon OpenFOAM/flow_rate_parametrization/probe_case_${i}_kOmega
    mapFields -consistent OpenFOAM/flow_rate_parametrization/probe_case_${i}_kEpsilon
    potentialFoam -writep -writePhi
    simpleFoam 
