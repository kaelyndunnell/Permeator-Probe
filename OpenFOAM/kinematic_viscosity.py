# calculate the kinematic viscosity of a fluid at a given temp

import numpy as np
from inlet_velocity import LiPb_density, breeder_temperature

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
