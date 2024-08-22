
""""
Created on 15.05.2023
@author: poudel-b
"""

# Mathieu and Bariac (1996) likelihood tests from 6.1 Braud_et_al 2005a

# -*- coding: utf-8 -*-

import cmf
import numpy as np

from iso_time import iso_time
from iso_layers import iso_layer
from Boundary_conditions import BoundaryCondition, iso_atmosphere
import solve_iso_transport as iso
from Visualize import Visualize


def main():
    # initialize
    # units in kg - m - days

    # Likelihood Tests (Mathieu and Bariac (1996).)
    ignorealphai = True
    ignorealphaik = True
    ignoredl = True
    ignoredv = True

    my_time = i_time()
    my_layers = layers()
    flux = fluxes(my_time)
    atm = atmosphere()
    B_C = BC(my_time)

    #  Stability Check
    print(check_stability(my_layers, my_time, flux))

    # Run simulations
    C_t = iso.run_1D_model(atm, my_layers, flux, B_C, solutes=['2H'],
                           ignore_alpha_i=ignorealphai,
                           ignore_alpha_i_k=ignorealphaik,
                           ignore_dl_i=ignoredl,
                           ignore_dv_i=ignoredv)

    C2H = C_t['2H'][0:50]
   #C18O = C_t['18O'][0:50]

    # Convert Result into delta notation
    my_vectorized_function = np.vectorize(iso.concentration_to_delta, excluded=['solute_i'])
    C2H_delta = my_vectorized_function(C2H, solute_i='2H')
    #C18O_delta = my_vectorized_function(C18O, solute_i='18O')

    #Visualize
    plot = Visualize(my_layers[0:50], my_time)
    plot.profile(C2H_delta, solute_i='2H', print_time_steps=24*30)
    plot.breakthrough(C2H_delta, solute_i='2H', print_steps=10)


def i_time():

    time_i = iso_time(time_units='hours',
                      delta_time=1,  # must be in given time_units
                      final_time=24*300  # must be in given time units
                      )
    return time_i


def layers():

    # Layer information
    top_layer_thickness = 0.000001   # Top layer of 1e-6 for numerical stability (Braud_et_al-2005)
    soil_depth = 1  # m
    dz = 0.01  # m
    layers = soil_depth / dz  # No of discrete layers

    initial_c_2H_soil = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
    initial_c_18O_soil = iso.delta_to_concentration(delta_i=-8, solute_i="18O")

    my_layers = []
    # Top Layer
    Top_layer = iso_layer(upper_boundary=0,
                          lower_boundary=top_layer_thickness,
                          initial_c_solutes={"2H": initial_c_2H_soil, "18O": initial_c_18O_soil},
                          initial_theta=0.30,  # volumetric water content m3/m3
                          theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                          theta_sat=0.35,  # volumetric water content m3/m3 at saturation
                          initial_T=303,  # soil temperature in Kelvin
                          porosity=0.35,  # porosity of the soil m3/m3
                          tortuosity=0.67)  # tortuosity of the soil m/m))
    my_layers.append(Top_layer)

    upper_boundary = top_layer_thickness
    for layer in range(0, int(layers)):
        lower_boundary = upper_boundary + dz
        new_layer = iso_layer(upper_boundary=upper_boundary,
                 lower_boundary=lower_boundary,
                 initial_c_solutes={"2H": initial_c_2H_soil, "18O": initial_c_18O_soil},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 initial_theta=0.30,  # volumetric water content m3/m3
                 theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                 theta_sat=0.35,  # volumetric water content m3/m3 at saturation
                 initial_T=303,  # soil temperature in Kelvin
                 porosity=0.35,  # porosity of the soil m3/m3
                 tortuosity=0.67)  # tortuosity of the soil m/m)

        my_layers.append(new_layer)
        upper_boundary = lower_boundary

    return my_layers


def fluxes(time):

    # Flux information
    ql = [1/1000/86400*0] * time.time_steps  # m / time_units
    qv = [0] * time.time_steps   # m / time_units

    return [ql, qv]


def atmosphere():

    initial_c_2H_atmosphere = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
    initial_c_18O_atmosphere = iso.delta_to_concentration(delta_i=-15, solute_i="18O")
    atm = iso_atmosphere(initial_c_atmosphere={"2H": initial_c_2H_atmosphere, "18O": initial_c_18O_atmosphere},
                         initial_T_atmosphere=303,  # temperature in the atmosphere in Kelvin
                         initial_Rh_atmosphere=0.20,  # relative humidity of the atmosphere (-)
                         initial_Pa=10 ** 5,  # atmospheric pressure (Pa)
                         initial_wind_speed=2.0,  # wind speed at top of canopy in m/s
                         initial_Precipitation=2.0,  # Precipitation in m3/d
                         initial_c_Precipitation={"2H": 1.0, "18O": 1.0},
                         # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                         )
    return atm


def BC(time):

    # Boundary conditions
    pot_evaportion = 1.005 * 10**-5  # kg / m2 / time

    U_boundary_conc = iso.delta_to_concentration(delta_i=-0, solute_i='2H')
    L_boundary_conc = iso.delta_to_concentration(delta_i=-65, solute_i='2H')

    BC = BoundaryCondition()
    BC.upper_boundary('atmosphere', [pot_evaportion] * time.time_steps)  # neuman[flux] / dirichlet[constant conc] / atmosphere[potential evaporation]
    BC.lower_boundary('dirichlet', [L_boundary_conc] * time.time_steps)   # TODO: need to check lower boundary

    return BC


def check_stability(layers, time, flux):

    for layer in layers:
        stable = True
        for ql, qv in zip(*flux):
            if not (stability(ql, time.dt, layer.thickness) and stability(qv, time.dt, layer.thickness)):
                stable = False
                break
        if stable:
            return 'Stability Checked! Stable'
        else:
            return 'Stability Checked! Unstable'


def stability(q, dt, dz):
    if q * dz / dt <= 0.5:
        return True
    else:
        return False


if __name__ == '__main__':
    main()



