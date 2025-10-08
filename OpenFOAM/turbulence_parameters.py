import numpy as np


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
    return (3 / 2) * (0.05 * 1.1e-2) ** 2  # 5% of inlet velocity


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
    return 0.09 ** (3 / 4) * calculate_initial_k() ** (3 / 2) / characteristic_length
