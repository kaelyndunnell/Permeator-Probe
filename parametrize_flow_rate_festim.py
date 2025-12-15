from two_volume_festim_probe_in import build_festim_model

N_A = 6.022e23
c_in = (
    1.5e-2 * N_A
)  # atms/m3 , inspired by tritium concentration (mol/m3) of OB loop from Utili 2023

openfoam_data_files = {
    2: (
        "OpenFOAM/flow_rate_parametrization/probe_in/probe_case_2_kOmega/case.foam",
        704,
    ),
    3: (
        "OpenFOAM/flow_rate_parametrization/probe_in/probe_case_3_kOmega/case.foam",
        753,
    ),
    4: (
        "OpenFOAM/flow_rate_parametrization/probe_in/probe_case_4_kOmega/case.foam",
        814,
    ),
}

for flow_rate, (openfoam_data_file, final_time) in openfoam_data_files.items():

    print(f"Running simulation for flow rate of {flow_rate} kg/s.")

    my_model = build_festim_model(
        openfoam_data_file=openfoam_data_file,
        openfoam_final_time=final_time,
        breeder_temperature=603.15,
        delta=0.1,
        c_in=c_in,
        Sc=0.7,
        results_folder=f"OpenFOAM_abacus/festim_results_flow_{flow_rate}_kg_s",
        insulated=True,
        visualize_fields=True,
    )

    # INITIALISE AND RUN
    my_model.initialise()
    my_model.run()

    print(f"Simulation for flow rate {flow_rate} kg/s complete.")
