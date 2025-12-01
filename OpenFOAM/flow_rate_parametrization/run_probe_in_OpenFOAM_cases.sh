#!/bin/bash

for i in $(seq 2 4); do # flow rates, same as in parametrization model

    cd probe_in/probe_case_${i}_kEpsilon

    gmshToFoam probe_in_openfoam_mesh.msh
    checkMesh

    cd ..
    cd ..
    python3 change_patch_names.py --case $i --turbulence kEpsilon --geometry probe_in
    cd probe_in/probe_case_${i}_kEpsilon
    potentialFoam -writep -writePhi

    simpleFoam

    echo "Ran Probe Case $i kEpsilon"

    cd .. 
    cd ..
    cd probe_in/probe_case_${i}_kOmega
    gmshToFoam probe_in_openfoam_mesh.msh
    checkMesh

    cd ..
    cd ..
    python3 change_patch_names.py --case $i --turbulence kOmega --geometry probe_in
    # TODO: create python script to change start time for kOmega depending on endtime for kEpsilon
    cd probe_in/probe_case_${i}_kOmega
    
    simpleFoam 

    echo "Running Probe Case $i kOmega"

done