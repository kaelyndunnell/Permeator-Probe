#!/bin/bash

for i in $(seq 2 4); do # flow rates, same as in parametrization model

    cd probe_case_${i}_kEpsilon

    gmshToFoam probe_in_openfoam_mesh.msh
    checkMesh

    cd ..
    python3 change_patch_names.py --case $i --turbulence kEpsilon
    cd probe_case_${i}_kEpsilon
    potentialFoam -writep -writePhi

    simpleFoam

    echo "Ran Probe Case $i kEpsilon"

    cd .. 
    cd probe_case_${i}_kOmega
    gmshToFoam openfoam_mesh.msh
    checkMesh

    cd ..
    python3 change_patch_names.py --case $i --turbulence kOmega
    # TODO: create python script to change start time for kOmega depending on endtime for kEpsilon
    cd probe_case_${i}_kOmega
    
    simpleFoam 

    echo "Running Probe Case $i kOmega"

done