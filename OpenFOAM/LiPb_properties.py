from fluid_parameters import (
    calculate_inlet_velocity,
    calculate_kinematic_viscosity,
    calculate_reynolds_number,
    plot_reynolds_number_vs_inlet_velocity,
    calculate_initial_k,
    calculate_initial_epsilon,
)
import numpy as np
import festim as F

breeder = "LiPb"

flow_rate = 1  # kg/s ; from Utili 2023
inlet_diameter = 0.13  # m from CAD

breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019
tube_diameter = 13e-2  # m, diameter of tube from CAD

k_b = F.k_B  # eV/K, boltzmann constant
E_D = 19500 * 1.0364e-5  # = 0.202098
LiPb_diffusivity = 4.03e-8 * np.exp(
    -E_D / (k_b * breeder_temperature)
)  # m2/s ; from Utili 2023, 1 J/mol = 1.0364E-5eV

inlet_velocity = calculate_inlet_velocity(
    flow_rate, inlet_diameter, LiPb_density, breeder
)

kinematic_viscosity = calculate_kinematic_viscosity(
    breeder_temperature, LiPb_density, breeder
)

Re = calculate_reynolds_number(
    inlet_velocity, tube_diameter, kinematic_viscosity, breeder
)
k = calculate_initial_k(inlet_velocity)
epsilon = calculate_initial_epsilon(k, characteristic_length=tube_diameter)

print(f"Initial turbulence kinetic energy for {breeder}: {k} m2/s2")
print(f"Initial turbulence dissipation rate for {breeder}: {epsilon} m2/s3")

plot_reynolds_number_vs_inlet_velocity(
    tube_diameter, kinematic_viscosity, breeder_temperature, breeder
)
