# script to produce mesh convergence study using OpenFOAM simulations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# test three meshes to start

finer_mesh = {
    "Elements": 3332802,
    "Max_Velocity": 1.0e-2,
    "Max_Pressure": 1.1e-4,
    "Final_Pressure": 3.9e-5,
}
starting_mesh = {
    "Elements": 1593795,
    "Max_Velocity": 1.0e-2,
    "Max_Pressure": 1.0e-4,
    "Final_Pressure": 3.9e-5,
}
coarser_mesh = {
    "Elements": 702995,
    "Max_Velocity": 1.0e-2,
    "Max_Pressure": 1.1e-4,
    "Final_Pressure": 3.2e-5,
}
coarsest_mesh = {
    "Elements": 585602,
    "Max_Velocity": 1.0e-2,
    "Max_Pressure": 1.0e-4,
    "Final_Pressure": 3.2e-5,
}


elements = [
    finer_mesh["Elements"],
    starting_mesh["Elements"],
    coarser_mesh["Elements"],
    coarsest_mesh["Elements"],
]
max_velocity = [
    finer_mesh["Max_Velocity"],
    starting_mesh["Max_Velocity"],
    coarser_mesh["Max_Velocity"],
    coarsest_mesh["Max_Velocity"],
]
max_pressure = [
    finer_mesh["Max_Pressure"],
    starting_mesh["Max_Pressure"],
    coarser_mesh["Max_Pressure"],
    coarsest_mesh["Max_Pressure"],
]
final_pressure = [
    finer_mesh["Final_Pressure"],
    starting_mesh["Final_Pressure"],
    coarser_mesh["Final_Pressure"],
    coarsest_mesh["Final_Pressure"],
]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

ax1.plot(elements, max_velocity)
ax2.plot(elements, max_pressure)
ax3.plot(elements, final_pressure)

formatter = ticker.ScalarFormatter(useOffset=False, useMathText=True)
formatter.set_powerlimits((0, 0))  # sci fi notation for y axis
ax1.yaxis.set_major_formatter(formatter)
ax2.yaxis.set_major_formatter(formatter)
ax3.yaxis.set_major_formatter(formatter)

ax1.set_ylabel("Max Velocity (m/s)")
ax1.set_title("Mesh Convergence Study")
ax2.set_ylabel("Max Pressure (m2/s2)")
ax2.set_xlabel("Number of elements")
ax3.set_ylabel("Final Pressure (m2/s2)")

plt.tight_layout()
plt.show()
