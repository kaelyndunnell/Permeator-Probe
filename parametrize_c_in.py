from dolfinx.log import set_log_level, LogLevel
from two_volume_festim_model import build_festim_model


N_A = 6.022e23

# mol_m3_c_in = [5e-9,5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,5e-1,5,5e1]  # mol/m3, corresponds to about 3e15 to 3e25 atms/m3
mol_m3_c_in = [6.5e-5, 1e2]  # mol/m3, corresponds to 1e19 and 6e25 #/m3

for conc in mol_m3_c_in:
    c_in = conc * N_A  # atms/m3

    print(f"Running simulation for inlet Concentration of {c_in} #/m3.")

    my_model = build_festim_model(
        openfoam_data_file="OpenFOAM/laminar-case/probe.foam",
        openfoam_final_time=300,
        breeder_temperature=603.15,
        delta=0.1,
        c_in=c_in,
        results_folder=f"festim_results_cin_{conc}_mol_m3",
        insulated=True,
        visualize_fields=True,
    )

    # INITIALISE AND RUN
    my_model.initialise()
    set_log_level(LogLevel.INFO)
    my_model.run()
