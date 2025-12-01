#!/bin/bash

for i in $(seq 2 4); do # flow rates, same as in parametrization model

    cd probe_out/probe_case_${i}_kEpsilon

    gmshToFoam probe_out_openfoam.msh
    checkMesh

    cd ..
    cd ..
    python3 change_patch_names.py --case $i --turbulence kEpsilon --geometry probe_out
    cd probe_out/probe_case_${i}_kEpsilon
    potentialFoam -writep -writePhi

    simpleFoam

    echo "Ran Probe Case $i kEpsilon"

    cd ..
    cd .. 
    cd probe_out/probe_case_${i}_kOmegaSST
    gmshToFoam probe_out_openfoam.msh
    checkMesh

    cd ..
    cd ..
    python3 change_patch_names.py --case $i --turbulence kOmegaSST --geometry probe_out
    cp -r probe_out/probe_case_${i}_kEpsilon probe_out/probe_case_${i}_kOmegaSST
    cd probe_out/probe_case_${i}_kOmegaSST
    
    mapFields -consistent probe_case_${i}_kEpsilon
    potentialFoam -writep -writePhi
    simpleFoam 

    echo "Ran Probe Case $i kOmegaSST."

done