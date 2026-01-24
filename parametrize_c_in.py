# from two_volume_festim_probe_out import build_festim_model

from two_volume_festim_probe_in import build_festim_model

N_A = 6.022e23

mol_m3_c_in = [
    5e-5,
    1e-4,
    5e-4,
    1e-3,
    5e-3,
    1e-2,
    1e-1,
    5e-1,
    1,
    5,
    1e1,
    5e1,
]  # mol/m3, corresponds to about 3e15 to 3e25 atms/m3

for conc in mol_m3_c_in:
    c_in = conc * N_A  # atms/m3

    print(f"Running simulation for inlet concentration of {c_in} #/m3.")

    # PROBE IN MODEL
    my_model = build_festim_model(
        openfoam_data_file="OpenFOAM/benchmark_cases_LiPb/kOmega-case/case.foam",
        openfoam_final_time=208,
        festim_mesh_file="meshing/probe_in_festim_mesh.msh",
        breeder_temperature=603.15,
        delta=10,
        c_in=c_in,
        Sc=0.7,
        results_folder=f"surfaceRecomboBC/festim_abacus_probe_in/festim_results_cin_{conc}_mol_m3",
        insulated=True,
        visualize_fields=False,
    )

    # PROBE OUT MODEL
    # my_model = build_festim_model(
    #     openfoam_data_file="OpenFOAM/benchmark_cases_LiPb/probe_out_kOmegaSST/case.foam",
    #     openfoam_final_time=2457,
    #     festim_mesh_file="meshing/probe_out_festim_mesh.msh",
    #     breeder_temperature=603.15,
    #     delta=0.1,
    #     c_in=c_in,
    #     Sc=0.7,  # seems to be default in OpenFOAM, find a reference to back up
    #     results_folder=f"surfaceRecomboBC/festim_abacus_probe_out/festim_results_cin_{conc}_mol_m3",
    #     insulated=True,
    #     visualize_fields=False,
    # )

    # INITIALISE AND RUN
    my_model.initialise()
    my_model.run()

    print(f"Simulation for inlet concentration {c_in} #/m3 complete.")
