import numpy as np
import matplotlib.pyplot as plt


def calculate_inlet_velocity(flow_rate, inlet_diameter, breeder_density, breeder):
    """Calculate the inlet velocity of fluid breeder at a given flow rate, inlet diameter, breeder density, and temperature.
    Used for OpenFOAM simulation.

    Parameters
    ----------
    flow_rate : float
        Flow rate in kg/s.
    inlet_diameter : float
        Inlet diameter in m.
    breder_density : float
        Breeder fluid density in kg/m3.
    breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Inlet velocity in m/s.
    """

    inlet_area = np.pi * (inlet_diameter / 2) ** 2  # m^2

    inlet_velocity = flow_rate * breeder_density ** (-1) * inlet_area ** (-1)  # m/s
    print(
        f"Inlet velocity for {breeder} flow rate of {flow_rate}kg/s is {inlet_velocity}m/s."
    )

    return inlet_velocity


def calculate_kinematic_viscosity(breeder_temperature, breeder_density, breeder):
    """Calculate the kinematic viscosity of a fluid breeder at a given temperature and density.
    Used for OpenFOAM simulation.

    Parameters
    ----------
    breeder_temperature : float
        Breeder temperature in K.
    breeder_density : float
        Breeder density in kg/m3.
     breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Kinematic viscosity in m2/s.
    """
    breeder_dynamic_viscosity = (
        0.0061091
        - 2.2574e-5 * breeder_temperature
        + 3.766e-8 * breeder_temperature**2
        - 2.2887e-11 * breeder_temperature**3
    )  # Pa.s = kg s / m s2; equation from Martelli 2019

    kinematic_viscosity = breeder_dynamic_viscosity / breeder_density  # m2/s

    print(
        f"Kinematic viscosity of {breeder} at {breeder_temperature}K is {kinematic_viscosity}m2/s."
    )

    return kinematic_viscosity


def calculate_reynolds_number(
    inlet_velocity,
    characteristic_length,
    kinematic_viscosity,
    breeder,
    suppress_print=False,
):
    """Calculate the reynolds number of a fluid breeder given its inlet velocity, characteristic length, and kinematic viscosity.
    Used to determine if turbulence is needed for OpenFOAM simulation.

    Parameters
    ----------
    inlet_velocity : float
        Inlet velocity in m/s.
    characteristic_length : float
        Characteristic length in m.
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
     breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Reynolds number (dimensionless).
    """

    reynolds_number = (inlet_velocity * characteristic_length) / kinematic_viscosity

    if not suppress_print:

        print(f"Reynolds number for {breeder} is {reynolds_number}.")

        if reynolds_number > 3500:
            print(f"Flow is turbulent.")
        else:
            print(f"Flow is laminar.")

    return reynolds_number


def calculate_peclet_number(inlet_velocity, characteristic_length, diffusivity):
    """Calculate the Peclet number of a fluid breeder given its inlet velocity, characteristic length, and diffusivity.
    Used to determine the relative importance of advection and diffusion in mass transport.

    Parameters    ----------
    inlet_velocity : float
        Inlet velocity in m/s.
    characteristic_length : float
        Characteristic length in m.
    diffusivity : float
        Diffusivity in m2/s.

    Returns
    -------
    float
        Peclet number (dimensionless).
    """
    peclet_number = diffusivity / (inlet_velocity * characteristic_length)
    print(f"Peclet number is {peclet_number} for a diffusivity of {diffusivity}m2/s.")
    return peclet_number


def plot_reynolds_number_vs_inlet_velocity(
    characteristic_length, kinematic_viscosity, breeder_temperature, breeder
):
    """Plot the reynolds number of a fluid breeder with varying inlet velocities.
    Assumes constant characteristic length and kinematic viscosity.

    Parameters
    ----------
    characteristic_length : float
        Characteristic length in m.
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
    breeder_temperature : float
        Breeder temperature in K.
     breeder : str
        Breeder fluid name.
    """
    inlet_velocities = np.linspace(0, 1e-2, 10000)  # m/s
    Re_numbers = []

    for inlet_velocity in inlet_velocities:  # m/s
        Re_numbers.append(
            calculate_reynolds_number(
                inlet_velocity,
                characteristic_length,
                kinematic_viscosity,
                breeder,
                suppress_print=True,
            )
        )

    plt.plot(inlet_velocities, Re_numbers, "b-")
    plt.axhline(y=3500, color="r", linestyle="--", label="Turbulence Threshold")
    plt.xlabel("Inlet Velocity (m/s)")
    plt.ylabel("Reynolds Number")
    plt.title(
        f"Reynolds Number vs Inlet Velocity for {breeder} at {breeder_temperature}K."
    )
    plt.legend()
    plt.show()


if __name__ == "__main__":

    breeder = "LiPb"

    flow_rate = 1  # kg/s ; from Utili 2023
    inlet_diameter = 0.13  # m from CAD

    breeder_temperature = 603.15  # K from Utili 2023
    LiPb_density = (
        10520.35 - 1.19051 * breeder_temperature
    )  # kg/m3 ; equation from Martelli 2019
    tube_diameter = 13e-2  # m, diameter of tube from CAD

    k_b = 8.617e-5  # eV/K, boltzmann constant
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

    plot_reynolds_number_vs_inlet_velocity(
        tube_diameter, kinematic_viscosity, breeder_temperature, breeder
    )

    calculate_peclet_number(inlet_velocity, tube_diameter, diffusivity=LiPb_diffusivity)
