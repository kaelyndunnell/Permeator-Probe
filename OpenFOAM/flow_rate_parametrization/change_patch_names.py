import os
import shutil
import argparse


def replace_specific_variables(filename, targets):
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    out = []
    inside = None

    for line in lines:
        stripped = line.strip()
        if stripped in targets:
            inside = stripped
        if stripped == "}":
            inside = None
        if inside and "type" in stripped and "patch" in stripped:
            line = line.replace("patch", "wall")
        out.append(line)

    with open(filename, "w", encoding="utf-8") as f:
        f.writelines(out)


# CHANGE 'PATCH' TO 'WALL' FOR WALLS AND PROBE IN CONSTANT/POLYMESH/BOUNDARY

parser = argparse.ArgumentParser(description="Determines which case to run in.")

parser.add_argument("--case", type=str, help="Number of OpenFOAM case being run.")
parser.add_argument(
    "--turbulence", type=str, help="Turbulence type, kEpsilon or kOmega."
)

args = parser.parse_args()

case_number = args.case
model = args.turbulence

replace_specific_variables(
    filename="probe_case_"
    + str(case_number)
    + "_"
    + str(model)
    + "/constant/polyMesh/boundary",
    targets={"wall", "probe"},
)
