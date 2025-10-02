# generate inlet velocity for OpenFOAM simulation using flow rate

import numpy as np

flow_rate = 1  # kg/s ; from Utili 2023
inlet_diameter = 0.13  # m from CAD
inlet_area = np.pi * (inlet_diameter / 2) ** 2  # m^2

breeder_temperature = 603.15  # K from Utili 2023
LiPb_density = (
    10520.35 - 1.19051 * breeder_temperature
)  # kg/m3 ; equation from Martelli 2019

inlet_velocity = flow_rate * LiPb_density ** (-1) * inlet_area ** (-1)  # m/s
print(f"Inlet velocity for LiPb flow rate of {flow_rate}kg/s is {inlet_velocity}m/s.")
