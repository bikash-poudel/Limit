'''
Created on 06.08.2024
@author: poudel-b
'''

# -*- coding: utf-8 -*-


def delta_to_concentration(delta_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}, M_w=0.018,
                           M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    R_i = delta_to_ratio(delta_i, solute_i, R_ref)  # isotopic ratio of isotopic species i in sample
    c_solutes = R_i * density_H2O * M_i[solute_i] / M_w
    return c_solutes


def delta_to_ratio(delta_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    R_i = (delta_i / 1000 + 1) * R_ref[solute_i]  # isotopic ratio of isotopic species i in sample
    return R_i


def concentration_to_delta(c_solute_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}, M_w=0.018,
                           M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    R_i = concentration_to_ratio(c_solute_i, solute_i, M_w, M_i, density_H2O)  # isotopic ratio of isotopic species i

    delta_i = ratio_to_delta(R_i, solute_i, R_ref)

    return delta_i


def concentration_to_ratio(c_solutes_i, solute_i, M_w=0.018, M_i={"2H": 0.019, "18O": 0.020},
                           density_H2O=1000.0):
    """
    Converts the given concentration of isotope species i  (kg/m**3) into the ratio of isotope species i

    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """

    R_i = c_solutes_i / (density_H2O * M_i[solute_i] / M_w)

    return R_i


def ratio_to_delta(R_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    delta_i = (R_i / R_ref[solute_i] - 1) * 1000
    return delta_i


def delta_testcases(Isotopologue, testcase=1):

    # delta signature in alimentation water and atmosphere

    delta_ali = {"2H": -65, "18O": -8}

    # delta signature in atmospheric vapor
    if testcase in [1, 3]:
        delta_atm = {"2H": -65, "18O": -8}
    elif testcase in [2, 4, 5, 6]:
        delta_atm = {"2H": -112, "18O": -15}

    else:
        raise ValueError

    return delta_ali[Isotopologue], delta_atm[Isotopologue]


def test_case_args(testcase=1, BA=False):

    # Testcases: Mathieu and Bariac (1996)
    ignore = {'ignoredvi': True,
              'ignoredli': True,
              'ignorealphai': True,
              'ignorealphaik': True,
              'Barnes_Allison': BA
              }

    # Condition: Barnes_Allison is only allowed in testcase 6
    if ignore['Barnes_Allison'] and testcase != 6:
        raise ValueError("Barnes_Allison can only be True if testcase is 6")

    if testcase == 1 or testcase == 2:
        return ignore
    elif testcase == 3 or testcase == 4:
        ignore['ignorealphai'] = False
        return ignore
    elif testcase == 5:
        ignore['ignorealphai'] = False
        ignore['ignoredli'] = False
        return ignore
    elif testcase == 6:
        ignore['ignorealphai'] = False
        ignore['ignoredli'] = False
        ignore['ignoredvi'] = False
        ignore['ignorealphaik'] = False

        return ignore
    else:
        raise NotImplementedError

