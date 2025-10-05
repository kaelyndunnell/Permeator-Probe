# script to calculate the reynolds number of a given fluid and thus if turbulence is needed

import numpy as np
from inlet_velocity import inlet_velocity
from kinematic_viscosity import kinematic_viscosity

fluid = "LiPb"  # fluid in tube

characteristic_length = 13e-2  # m, diameter of tube

reynolds_number = (inlet_velocity * characteristic_length) / kinematic_viscosity

print(
    f"Reynolds number for {fluid} is {reynolds_number}. Reynolds numbers above 3500 are generally considered to indicate turbulent flow."
)
