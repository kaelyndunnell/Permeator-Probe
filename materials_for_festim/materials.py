import h_transport_materials as htm


alpha_Fe_recombination = htm.recombination_coeffs.filter(material=htm.IRON)[0]

# print(alpha_Fe_recombination)

u = htm.ureg
solubility_alpha_iron = htm.Solubility(
    S_0=0.51 * u.mol * u.m**-3 * u.Pa**-0.5,
    E_S=27 * u.kJ * u.mol**-1,
    source="10.1016/S0022-3115(96)00670-8",
)

diffusivity_alpha_iron = (
    htm.diffusivities.filter(material=htm.IRON)
    .filter(exclude=True, isotope="H")
    .filter(exclude=True, isotope="D")
    .mean()
)
lipb_solubility = (
    htm.solubilities.filter(material=htm.LIPB)
    .filter(exclude=True, isotope="H")
    .filter(exclude=True, isotope="D")
    .mean()
)
lipb_diffusivity = (
    htm.diffusivities.filter(material=htm.LIPB)
    .filter(exclude=True, isotope="H")
    .filter(exclude=True, isotope="D")
    .mean()
)

# print(lipb_diffusivity)
# print(lipb_solubility)
# print(solubility_alpha_iron)
# print(diffusivity_alpha_iron)

u = htm.ureg
recombination_Nb = htm.RecombinationCoeff(
    pre_exp=1.88e-18,  # m4/s
    act_energy=74 * u.kJ * u.mol**-1,
    source="10.1238/Physica.Topical.103a00113",
)

# print(recombination_Nb)
