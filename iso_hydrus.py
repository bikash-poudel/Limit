'''
Created on 03.06.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

import hydrus
from src import *

import matplotlib.pyplot as plt


def _layers(c, hy, testcase):
    Tzero_sli = 273.16  # [k] 0 celcius in kelvin

    dt = 0  # initial states
    lower_boundaries = hy.depth

    theta_r = hy.thetar
    theta_sat = hy.theta_sat
    tortuosity = [0.67] * len(hy.depth)

    theta = hy.theta(dt)
    T_soil = hy.T(dt)
    psi = hy.head(dt)
    # rH = sli.R_humidity(dt)

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    init_c_iso_2H = [iso_delta.delta_to_concentration(ali_2H, '2H')] * len(lower_boundaries)
    init_c_iso_18O = [iso_delta.delta_to_concentration(ali_18O, '18O')] * len(lower_boundaries)

    id = 0
    upper_boundary = 0
    for lower_boundary, c2H, c18O, th, tr, tsat, tor, T, PSI in \
            zip(lower_boundaries, init_c_iso_2H, init_c_iso_18O, theta, theta_r, theta_sat, tortuosity, T_soil, psi):
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=upper_boundary,
                                                lower_boundary=lower_boundary,
                                                conc_iso_liquid={"2H": c2H, "18O": c18O},
                                                theta=th,
                                                theta_0=tr,
                                                theta_sat=tsat,
                                                tortuosity=tor,
                                                T=T + Tzero_sli,
                                                psi=PSI
                                                )
        id += 1
        upper_boundary = lower_boundary
        c.add_layer(new_layer)

    return c


def _atm(hy, testcase):
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = hy.T_s() + Tzero_sli

    Rh_atm = 0.2  # sli.R_humidity_atm(dt)
    wind_speed = 2.0  # sli.wind_speed(dt)

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=2,  # 10,  # canopy height [m] (e.g. 40)
                         d0=0.67 * 2,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.01,  # 0.01
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=0,  # 1.1,  # 1.1  # leaf area index (e.g. 2.0)
                         extku=0)
    return atm


def update_storages(c, hy, dt):
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Need not update atmosphere
    """"
    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = sli.Tatm(dt) + Tzero_sli
    Rh_atm = sli.R_humidity_atm(dt)
    c_iso_18O = sli.civa(dt) / sli.cva(dt)
    c_iso_2H = 0.15367056287906838
    wind_speed = sli.wind_speed(dt)
    c.update_atmosphere(c_atm={'2H': c_iso_2H, '18O': c_iso_18O}, T=Tatm, Rh=Rh_atm, Pa=patm, wind_speed=wind_speed)
    """

    # Update layers
    theta = hy.theta(dt)
    T_soil = np.array(hy.T(dt)) + Tzero_sli
    rH = None
    psi = hy.head(dt)  # sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c, hy, dt):
    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Surface variables
    T_surface = hy.T_surface(dt) + Tzero_sli  ## TODO: recheck the surface temperature value ?????
    q_ev = hy.q_evap(dt) # sli.qevap(dt)]

    ql_surface = hy.ql_surface(dt)
    qv_surface = hy.qv_surface(dt)
    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=0.0, qv_surface=0.0)

    # soil fluxes
    ql = np.array(hy.q_liquid(dt)[1:])
    qv = np.array(hy.q_vapor(dt)[1:])

    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)


def run_iso(p, hy, **kwargs):
    solutes = ["2H", '18O']
    c = p.get_cells()[0]  # get current cell of project

    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}

    for delta_t, dt, in zip(hy.dt, hy.time[1:]):  # starting from next time step (n+1)

        print(dt, '  days')

        # c_iso["2H"].append(c.conc_2H), c_iso["18O"].append(c.conc_18O)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        # update storage states and boundaries to current time
        update_storages(c, hy, dt), update_boundaries(c, hy, dt)
        print([l.theta for l in c.layers])

        for solute in solutes:

            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=None, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


def iso_setup(hy, testcase=1, **kwargs):
    # Define boundary storages
    atm = _atm(hy, testcase)  # atmosphere

    p = iso_project()  # create a project
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell

    c = p.get_cells()[0]  # get current cell
    _layers(c, hy, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation(),

    return run_iso(p, hy, **kwargs)


path_H1D_files = 'D:\Hydrus\sli_vapor'
path_H1D = 'C:\Program Files (x86)\PC-Progress\Hydrus-1D 4.xx\H1D_CALC.EXE'
#hydrus.run_hydrus(path_H1D, path_H1D_files)

h = hydrus.hydrus(path_H1D_files)
ignore = iso_delta.test_case_args(testcase=1)
delta = iso_setup(hy=h, testcase=1, **ignore)
