import numpy as np

# generate inlet velocity for OpenFOAM simulation using flow rate
flow_rate = 1  # kg/s ; from Utili 2023
inlet_diameter = 0.13  # m from CAD
inlet_area = np.pi * (inlet_diameter / 2) ** 2  # m^2

breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019

inlet_velocity = flow_rate * LiPb_density ** (-1) * inlet_area ** (-1)  # m/s
print(f"Inlet velocity for LiPb flow rate of {flow_rate}kg/s is {inlet_velocity}m/s.")


# calculate the kinematic viscosity of a fluid at a given temp
LiPb_dynamic_viscosity = (
    0.0061091
    - 2.2574e-5 * breeder_temperature
    + 3.766e-8 * breeder_temperature**2
    - 2.2887e-11 * breeder_temperature**3
)  # Pa.s = kg s / m s2; equation from Martelli 2019

kinematic_viscosity = LiPb_dynamic_viscosity / LiPb_density  # m2/s

print(
    f"Kinematic viscosity of LiPb at {breeder_temperature}K is {kinematic_viscosity}m2/s."
)

# script to calculate the reynolds number of a given fluid and thus if turbulence is needed
fluid = "LiPb"  # fluid in tube

characteristic_length = 13e-2  # m, diameter of tube

reynolds_number = (inlet_velocity * characteristic_length) / kinematic_viscosity

print(
    f"Reynolds number for {fluid} is {reynolds_number}. Reynolds numbers above 3500 are generally considered to indicate turbulent flow."
)
