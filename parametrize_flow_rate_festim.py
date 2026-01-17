# from two_volume_festim_probe_in import build_festim_model
from two_volume_festim_probe_out import build_festim_model

N_A = 6.022e23
c_in = (
    1e01 * N_A
)  # atms/m3 , inspired by tritium concentration (mol/m3) of OB loop from Utili 2023

# PROBE OUT DATA FILES
openfoam_data_files = {
    2: (
        "OpenFOAM/flow_rate_parametrization/probe_out/probe_case_2_kOmegaSST/case.foam",
        5108,
    ),
    3: (
        "OpenFOAM/flow_rate_parametrization/probe_out/probe_case_3_kOmegaSST/case.foam",
        4269,
    ),
    4: (
        "OpenFOAM/flow_rate_parametrization/probe_out/probe_case_4_kOmegaSST/case.foam",
        810,
    ),
}

# # PROBE IN DATA FILES
# openfoam_data_files = {
#     2: (
#         "OpenFOAM/flow_rate_parametrization/probe_in/probe_case_2_kOmega/probe.foam",
#         4812,
#     ),
#     3: (
#         "OpenFOAM/flow_rate_parametrization/probe_in/probe_case_3_kOmega/probe.foam",
#         5707,
#     ),
#     4: (
#         "OpenFOAM/flow_rate_parametrization/probe_in/probe_case_4_kOmega/case.foam",
#         5686,
#     ),
# }

for flow_rate, (openfoam_data_file, final_time) in openfoam_data_files.items():
    print(f"Running simulation for flow rate of {flow_rate} kg/s.")

    # # PROBE IN MODEL
    # my_model = build_festim_model(
    #     openfoam_data_file=openfoam_data_file,
    #     openfoam_final_time=final_time,
    #     festim_mesh_file="meshing/probe_in_festim_mesh.msh",
    #     breeder_temperature=603.15,
    #     delta=10,
    #     c_in=c_in,
    #     Sc=0.7,  # seems to be default in OpenFOAM, 10.1007/s10652-005-5656-9
    #     results_folder=f"surfaceRecomboBC/OpenFOAM_abacus_probe_in/festim_results_flow_{flow_rate}_kg_s",
    #     insulated=True,
    #     visualize_fields=True,
    # )

    # PROBE OUT MODEL
    my_model = build_festim_model(
        openfoam_data_file=openfoam_data_file,
        openfoam_final_time=final_time,
        festim_mesh_file="meshing/probe_out_festim_mesh.msh",
        breeder_temperature=603.15,
        delta=0.1,
        c_in=c_in,
        Sc=0.7,
        results_folder=f"surfaceRecomboBC/OpenFOAM_abacus_probe_out/festim_results_flow_{flow_rate}_kg_s",
        insulated=True,
        visualize_fields=True,
    )

    # INITIALISE AND RUN
    my_model.initialise()
    my_model.run()

    print(f"Simulation for flow rate {flow_rate} kg/s complete.")
