import numpy as np
from fluid_parameters import calculate_inlet_velocity

# reference: https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-epsilon.html


def calculate_initial_k(inlet_velocity):
    """Calculate initial turbulence kinetic energy.

    Parameters
    ----------
    inlet_velocity : float
        Inlet velocity in m/s.

    Returns
    -------
    float
        Initial turbulence kinetic energy in m2/s2.
    """

    # assume initial turbulence is isotropic so that U'2_x = U'2_y = U'2_z
    return (3 / 2) * (0.05 * inlet_velocity) ** 2  # 5% of inlet velocity


def calculate_initial_epsilon(k, characteristic_length):
    """Calculate initial turbulence dissipation rate.

    Parameters
    ----------
    k : float
        Initial turbulence kinetic energy in m2/s2.
    characteristic_length : float
        Characteristic length in m.

    Returns
    -------
    float
        Initial turbulence dissipation rate in m2/s3.
    """
    return 0.09 ** (3 / 4) * k ** (3 / 2) / characteristic_length


if __name__ == "__main__":
    breeder = "LiPb"

    flow_rate = 1  # kg/s ; from Utili 2023
    inlet_diameter = 0.13  # m from CAD

    breeder_temperature = 603.15  # K from Utili 2023
    LiPb_density = (
        10520.35 - 1.19051 * breeder_temperature
    )  # kg/m3 ; equation from Martelli 2019

    tube_diameter = 13e-2  # m, diameter of tube from CAD

    inlet_velocity = calculate_inlet_velocity(
        flow_rate, inlet_diameter, LiPb_density, breeder
    )  # m/s

    k = calculate_initial_k(inlet_velocity)
    epsilon = calculate_initial_epsilon(k, characteristic_length=tube_diameter)

    print(f"Initial turbulence kinetic energy for {breeder}: {k} m2/s2")
    print(f"Initial turbulence dissipation rate for {breeder}: {epsilon} m2/s3")
