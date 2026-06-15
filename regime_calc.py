import numpy as np

# example case is 10 mol/m3 probe-in

# flux = 2.767666454849072e18 # 10 mol m-3 probe in
flux = 4.781980280397347e+17 # 10 mol m-3 probe out
# flux = 40398548257366.88 # 10e-3 mol m-3, 4e13 probe in 

probe_thickness = 0.02e-2  # m
D_probe = 9.38e-9
D_breeder = 1e-9
K_breeder = 1.62e21
K_probe = 3.07e23
# diameter = 13.8e-2 # m
diameter = 6e-2 # m probe out compartment diameter
# mass_trans_coeff = (
#     5.1e-3  # https://www.sciencedirect.com/science/article/abs/pii/092037969190064W
# )

c_in = 10 * 6.022e23 # high case
# c_in = 10e-3 * 6.022e23 # low case
# c_b_interface_out = 6e23
# c_b_interface_in = 5.9e24
# c_p_interface_out = 1.1e26
# c_vac_out = 1.2e25
# mass_trans_coeff = abs(
#     D_probe
#     / probe_thickness
#     * (c_p_interface_out - c_vac_out)
#     / (c_in - c_b_interface_out)
# )
# print(f"mass trans coeff is {mass_trans_coeff}")

# calc mass trans coeff

viscosity = 1.1956969781435462e-07  # kinematic
v_in = 7.6e-3  # m/s

Re = v_in * diameter / viscosity
Sc = 0.7 #viscosity / D_breeder  # 0.7
# Sh = 0.3 + (
#     0.62 * Re ** (1 / 2) * Sc ** (1 / 3) * (1 + (0.4 / Sc) ** (2 / 3)) ** (-1 / 4)
# ) * (1 + (Re / 28200) ** (5 / 8))
Sh = 0.023 * Re ** (0.8) * Sc ** (1 / 3)

mass_trans_coeff = Sh * D_breeder / diameter
print(f"mass transfer coefficient is {mass_trans_coeff}")

# pressure_breeder = (c_in / K_breeder) ** 2

K_r = 1.26e-25  # recombo coeff
# h = 0.24e-2  # m

partition_parameter = D_probe / mass_trans_coeff * K_probe / K_breeder / probe_thickness
print(f"partition parameter is {partition_parameter}")

# W = (
#     (partition_parameter + 1) ** 2
#     * flux
#     * D_probe
#     * K_probe
#     * np.sqrt(pressure_breeder)
#     / probe_thickness
#     - (partition_parameter + 1)
# ) ** (-1) / 4

######### UROGORRI 2023 MODEL #########


W_0 = K_r * c_in * probe_thickness / (D_probe * K_breeder / K_probe)

print(f"mixed permeation number is {W_0}")

calculated_flux = (  # assuming LLR from Alberghi
    mass_trans_coeff * c_in * partition_parameter / (partition_parameter )
)

print(f"calculated permeation flux = {calculated_flux:.2e}")