from fluid_parameters import (
    calculate_inlet_velocity,
    calculate_reynolds_number,
    calculate_initial_k,
    calculate_initial_epsilon,
    calculate_initial_omega,
)
from LiPb_properties import calculate_LiPb_kinematic_viscosity
import numpy as np
import festim as F
import os
import shutil


def change_variable_in_openfoam_file(filename, old_value, new_value):
    with open(filename, "r", encoding="utf-8") as f:
        content = f.read()

    content = content.replace(str(old_value), str(new_value))

    with open(filename, "w") as f:
        f.write(content)


# CONSTANT BREEDER PARAMETERS
breeder = "LiPb"
inlet_diameter = 0.13  # m from CAD
breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019
k_b = F.k_B  # eV/K, boltzmann constant
E_D = 19500 * 1.0364e-5  # = 0.202098
LiPb_diffusivity = 4.03e-8 * np.exp(
    -E_D / (k_b * breeder_temperature)
)  # m2/s ; from Utili 2023, 1 J/mol = 1.0364E-5eV
kinematic_viscosity = calculate_LiPb_kinematic_viscosity(
    breeder_temperature, LiPb_density, breeder
)

# FLOW RATE PARAMETRIZATION
flow_rates = [1, 2, 3, 4]  # kg/s ; from Utili 2023

for rate in flow_rates:
    inlet_velocity = calculate_inlet_velocity(
        rate, inlet_diameter, LiPb_density, breeder
    )
    Re = calculate_reynolds_number(
        inlet_velocity,
        inlet_diameter,
        kinematic_viscosity,
        breeder,
        suppress_print=True,
    )
    k = calculate_initial_k(inlet_velocity)
    # epsilon = calculate_initial_epsilon(k, characteristic_length=inlet_diameter)
    omega = calculate_initial_omega(k, inlet_diameter)

    # TODO: add suppress print function

    openfoam_folder = f"OpenFOAM/flow_rate_parametrization/probe_case_{rate}_kg_s-1"
    os.makedirs(openfoam_folder, exist_ok=True)

    shutil.copytree(
        "OpenFOAM/kOmega-case/0/", openfoam_folder + "/0"
    )  # p, nut files are the same as benchmark kOmega case
    shutil.copytree("OpenFOAM/kOmega-case/system/", openfoam_folder + "/system")
    shutil.copytree("OpenFOAM/kOmega-case/constant/", openfoam_folder + "/constant")

    variables_dict = {
        "U": [0.0076, inlet_velocity],
        "k": [2.215e-07, k],
        "omega": [0.00661, omega],
    }

    for name, values in variables_dict.items():
        change_variable_in_openfoam_file(
            filename=openfoam_folder + "/0/" + name,
            old_value=values[0],
            new_value=values[1],
        )

# then need to have some way to run the openfoam simulations directly from this file
