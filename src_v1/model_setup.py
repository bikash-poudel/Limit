
""""
Created on 15.05.2023
@author: poudel-b
"""

# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

from iso_layers import iso_layer
from Boundary_conditions import BoundaryCondition, iso_atmosphere
import solve_iso_transport as iso_solve


def stability(q, dt, dz):
    if q * dz / dt * 24 <= 0.5:
        return True
    else:
        return False


def check_stability(ql, qv, dt, dz):
    stable = True
    for ql, qv in zip(ql, qv):
        if not (stability(ql, dt, dz) and stability(qv, dt, dz)):
            stable = False
            break
    if stable:
        return 'Stability Checked! Stable'
    else:
        return 'Stability Checked! Unstable'


def main():

    # TODO: check for diffusion time scale

    # initialize
    # units in kg - m - days

    # Layer information
    upper_boundary = 0
    soil_depth = 1  # m
    layers = 1000  # No of discrete layers
    dz = soil_depth / layers  # m

    my_layers = []
    for layer in range(layers):

        lower_boundary = upper_boundary + dz
        new_layer = iso_layer(upper_boundary=upper_boundary,
                 lower_boundary=lower_boundary,
                 initial_c_solutes={"2H": 0.1537265, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 initial_theta=0.25,  # volumetric water content m3/m3
                 theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                 theta_sat=0.35,  # volumetric water content m3/m3 at saturation
                 initial_T=283.15,  # soil temperature in Kelvin
                 porosity=0.35,  # porosity of the soil m3/m3
                 tortuosity= 0.67 )  # tortuosity of the soil m/m)

        my_layers.append(new_layer)
        upper_boundary = lower_boundary

    # Time information
    time_final = 300  # days
    dt = 1/86400 * 60 * 60   # 1 hour
    time_steps = int(time_final / dt)

    # Flux information
    ql = [0.000000001] * time_steps  # m / day  [ql / 86400  per sec] dz / dt
    qv = [0.0] * time_steps  # m/day

    #  Stability Check
    print(check_stability(ql, qv, dt, dz))

    # Boundary conditions
    qi_ev = [0.00000000015] * time_steps     # evaporation flux kg / day
    qi_pcp = [0.0] * time_steps  # [0.3] * 60 + [0.2] * 60 + [0.1] * 120  # precipitation flux kg / day (variable flow)
    qi_upper = np.add(qi_pcp, qi_ev).tolist()

    BC = BoundaryCondition()
    BC.upper_boundary('neuman', qi_upper)    # neuman[flux] / dirichlet[constant conc]
    BC.lower_boundary('neuman', [0.0] * time_steps)

    time = [time_steps, dt]
    Q = [ql, qv]

    # Run simulations
    C_t = iso_solve.run_1D_model(atmosphere(), my_layers, time, Q, BC, solutes=['2H'])
    C2H = C_t['2H']
    # C18O = C_t['18O']

    # Convert Result into delta signature
    my_vectorized_function = np.vectorize(iso_solve.concentration_to_delta, excluded=['solute_i'])
    C2H_dlta = my_vectorized_function(C_t['2H'], solute_i='2H')
    C2H_delta = np.round(C2H_dlta, 2)

    # VISUALIZATION
    # Concentration vs time
    fig1, ax = plt.subplots(figsize=(10, 5), dpi=300)
    for index, value_index in enumerate(range(0, len(my_layers), 200)):
        ax.plot(np.arange(0, time_steps) * dt, C2H_delta[value_index],
                label="{:.1f}m".format(my_layers[value_index].center))
    # ax.plot(np.arange(0, time_steps) * dt, C2H[-1], label="lower boundary")

    ax.set_xlabel('time [days]')
    ax.set_ylabel('delta')
    ax.set_title('conc breakthrough [2H]')
    # plt.yscale('log')
    ax.legend()

    # timestamp1 = datetime.now().strftime("%Y%m%d-%H%M%S")
    # fig1.savefig( 'D:\\Isotope transport\\Scripts\\output\\{}.png'.format(timestamp1))
    plt.show()

    # Concentration vs depth

    fig2, ax1 = plt.subplots(figsize=(10, 5), dpi=300)
    for t in range(0, time_steps, 24 * 20):
        ax1.plot(C2H_delta[:, t], -np.arange(0, len(my_layers)) * dz,
                 label='day {}'.format(t / 24 ))

    ax1.set_xlabel('delta')
    ax1.set_ylabel('depth [m]')
    ax1.set_title('conc profiles [2H]')
    # plt.yscale('log')
    ax1.legend()

    # timestamp2 = datetime.now().strftime("%Y%m%d-%H%M%S")
    # fig2.savefig( 'D:\\Isotope transport\\Scripts\\output\\{}.png'.format(timestamp2))
    plt.show()
    return


def atmosphere():

    initial_c_2H_atmosphere = iso_solve.delta_to_concentration(delta_i=-112, solute_i="2H")
    initial_c_18O_atmosphere = iso_solve.delta_to_concentration(delta_i=-8, solute_i="18O")
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




if __name__ == '__main__':
    main()

