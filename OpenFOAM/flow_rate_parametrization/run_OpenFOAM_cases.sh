#!/bin/bash

for i in $(seq 1 4); do # flow rates, same as in parametrization model

    cd probe_case_${i}_kEpsilon

    gmshToFoam openfoam_mesh.msh
    checkMesh

    cd ..
    python3 change_patch_names.py --case $i --turbulence kEpsilon
    cd probe_case_${i}_kEpsilon

    simpleFoam

    echo "Ran Probe Case $i kEpsilon"

    cd .. 
    cd probe_case_${i}_kOmega
    gmshToFoam openfoam_mesh.msh
    checkMesh

    cd ..
    python3 change_patch_names.py --case $i --turbulence kOmega
    cd probe_case_${i}_kOmega

    mv probe_case_${i}_kEpsilon probe_case_${i}_kOmega
    mapFields -consistent probe_case_${i}_kEpsilon
    potentialFoam -writep -writePhi
    simpleFoam 

    echo "Running Probe Case $i kOmega"

done