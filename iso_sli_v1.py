
""""
Created on 10.09.2024
@author: poudel-b
"""

# -*- coding: utf-8 -*-
import os
import Sli
import matplotlib.pyplot as plt
import numpy as np
from src import *
from src_v1 import Boundary_conditions as BC
from src_v1 import solve_iso_transport as iso


def main():

    sli = get_sli(testcase=6)
    atm = _atm(sli, testcase=1)
    layers = _layers(sli, testcase=1)
    solute = '2H'

    D = []
    for dt in range(1, len(sli.get_in_soil()) - 1):

        print(dt)
        layers = update_layers(layers, sli, dt)

        # Flux information
        ql = sli.qlsig(dt)[1:]  # m / day  [ql / 86400  per sec] dz / dt
        qv = sli.qvsig(dt)[1:]  # m/day

        # Boundary conditions
        qev = sli.qevap(dt)  # evaporation flux kg / day
        qi_pcp = sli.qprec(dt)  # [0.3] * 60 + [0.2] * 60 + [0.1] * 120  # precipitation flux kg / day (variable flow)
        # qi_upper = np.add(qi_pcp, qi_ev).tolist()

        Q = [ql, qv]

        bc = BC.BoundaryCondition()
        bc.upper_boundary('atmosphere', qev)  # neuman[flux] / dirichlet[constant conc]
        bc.lower_boundary('neuman', 0.0)

        delta_t = sli.dt(dt)
        c_t = iso.run_1D_model(atm, layers, Q, bc, delta_t, solute,
                               ignore_alpha_i=True,
                               ignore_alpha_i_k=True,
                               ignore_dl_i=True,
                               ignore_dv_i=True)

        update_concentration(layers, c_t, solute)

        D.append(delta(c_t, solute))

    return D


def get_sli(testcase=1):
    pth = os.getcwd()
    path = os.path.abspath(os.path.join(pth, "", ".."))

    # imports all the variable files /variables folder: testcase-1, sig=1
    sli = Sli.SlI(path + '\_sli_\sli_label3\iso_variables_{}'.format(testcase))

    return sli


def _layers(sli, testcase):
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    dt = 0  # initial states
    lower_boundaries = np.cumsum(sli.dx(dt))
    theta = sli.theta(dt)
    theta_r = sli.thetar(dt)
    theta_sat = sli.thetasat(dt)
    tortuosity = sli.tortuosity(dt)
    T_soil = sli.T_soil0(dt)
    rH = sli.R_humidity(dt)
    psi = sli.matric_pot(dt)

    if testcase in [1, 2, 3, 4, 5, 6]:
        init_c_iso_2H = [iso_delta.delta_to_concentration(-65, '2H')] * len(lower_boundaries)  # 0.15367056287906838
        init_c_iso_18O = [iso_delta.delta_to_concentration(-8, '18O')] * len(lower_boundaries)

    elif testcase in [7, 8]:
        init_c_iso_2H = [iso_delta.delta_to_concentration(0, '2H')] * len(lower_boundaries)  # 0.15367056287906838
        init_c_iso_18O = [iso_delta.delta_to_concentration(0, '18O')] * len(lower_boundaries)
    else:
        raise ValueError

    my_layers = []
    id = 0
    upper_boundary = 0
    for lower_boundary, c2H, c18O, th, tr, tsat, tor, T, r_H, PSI in \
            zip(lower_boundaries, init_c_iso_2H, init_c_iso_18O, theta, theta_r, theta_sat, tortuosity, T_soil, rH,
                psi):
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=upper_boundary,
                                                lower_boundary=lower_boundary,
                                                conc_iso_liquid={"2H": c2H, "18O": c18O},
                                                theta=th,
                                                theta_0=tr,
                                                theta_sat=tsat,
                                                tortuosity=tor,
                                                T=T + Tzero_sli,
                                                rH=r_H,
                                                psi=PSI
                                                )
        id += 1
        upper_boundary = lower_boundary
        my_layers.append(new_layer)
        # c.add_layer(new_layer)

    return my_layers


def _atm(sli, testcase):
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    dt = 0  # initial states
    # atmospheric variables
    patm = 100000  # from SLI_solve
    Tatm = sli.Tatm(dt) + Tzero_sli
    Rh_atm = sli.R_humidity_atm(dt)
    wind_speed = sli.wind_speed(dt)

    if testcase == 1 or testcase == 3:
        c_iso_2H = iso_delta.delta_to_concentration(-65, '2H')  # sli.civa(1) / sli.cva(1)  0.15367056287906838
        c_iso_18O = iso_delta.delta_to_concentration(-8, '18O')
    elif testcase in [2, 4, 5, 6]:
        c_iso_2H = iso_delta.delta_to_concentration(-112, '2H')  # sli.civa(dt) / sli.cva(dt)  0.15367056287906838
        c_iso_18O = iso_delta.delta_to_concentration(-15, '18O')  # sli.civa(dt) / sli.cva(dt)
    elif testcase in [7, 8]:
        c_iso_2H = iso_delta.delta_to_concentration(-100, '2H')  # sli.civa(dt) / sli.cva(dt)  0.15367056287906838
        c_iso_18O = iso_delta.delta_to_concentration(-14, '18O')  # sli.civa(dt) / sli.cva(dt)
    else:
        raise NotImplementedError

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=10,  # canopy height [m] (e.g. 40)
                         d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.1 * 10,
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=1.1,  # leaf area index (e.g. 2.0)
                         extku=1.5)
    return atm


def update_layers(layers, sli, dt):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision
    Temp = np.array(sli.T_soil0(dt)) + Tzero_sli

    for l, theta, T, rh, psi in zip(layers, sli.theta(dt), Temp, sli.R_humidity(dt), sli.matric_pot(dt)):

        # Update layers
        l.theta_t0 = l.theta
        l.theta = theta
        l.T = T
        l.rH = rh
        l.psi = psi

    return layers


def update_concentration(layers, conc, solute):

    for lr, c in zip(layers, conc):
        lr.set_conc_iso_liquid(c, solute)

    return layers


def delta(conc, solute):

    d = []
    for c in conc:
        d.append(iso_delta.concentration_to_delta(c, solute))

    return d


if __name__ == '__main__':
    main()

