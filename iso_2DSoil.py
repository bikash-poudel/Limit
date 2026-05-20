'''
Created on 03.12.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

from _2DSoil import _2D_Soil
from src import *

import matplotlib.pyplot as plt
import time


def plot_moisture(depth, soil, t):

    plt.plot(soil.theta(t)[1:], -depth[1:], label='theta')
    plt.legend()
    plt.show()


def plot_fluxes(depth, soil, t):

    plt.plot(soil.q_l(t)[1:], -depth[1:], label='ql')
    plt.plot(soil.q_v(t)[1:], -depth[1:], label='qv')
    plt.legend()
    plt.show()


def plot_delta(delta, depth, solute, test_cases):

    try:
        for test in test_cases:
            d = delta[test]
            plt.plot(d[solute], -depth[1:], label='testcase:{}'.format(test))

        plt.legend()
        plt.show()

    except ValueError:
        return NotImplementedError


def _atm(testcase):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = 30
    Rh_atm = 0.2
    wind_speed = 2  # m/s

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm + Tzero_sli, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=10,  # canopy height [m] (e.g. 40)
                         d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.1 * 10,
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=1.1,  # leaf area index (e.g. 2.0)
                         extku=1.5)
    return atm


def _layers(c, testcase):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision
    lower_boundaries = np.cumsum([0.05]*100)
    nx = len(lower_boundaries)
    theta, theta_r, theta_sat, tortuosity, T_soil = [0.5] * nx, [0.05] * nx, [0.547] * nx, [0.67] * nx, [25] * nx
    # rH = sli.R_humidity(dt)
    #psi = sli.matric_pot(dt)

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    init_c_iso_2H = [iso_delta.delta_to_concentration(ali_2H, '2H')] * len(lower_boundaries)
    init_c_iso_18O = [iso_delta.delta_to_concentration(ali_18O, '18O')] * len(lower_boundaries)

    id = 0
    upper_boundary = 0
    for lower_boundary, c2H, c18O, th, tr, tsat, tor, T in \
            zip(lower_boundaries, init_c_iso_2H, init_c_iso_18O, theta, theta_r, theta_sat, tortuosity, T_soil):

        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=upper_boundary,
                                                lower_boundary=lower_boundary,
                                                conc_iso_liquid={"2H": c2H, "18O": c18O},
                                                theta=th,
                                                theta_0=tr,
                                                theta_sat=tsat,
                                                tortuosity=tor,
                                                T=T + Tzero_sli,
                                                rH=None
                                                )
        id += 1
        upper_boundary = lower_boundary
        c.add_layer(new_layer)

    return c


def update_storages(c, soil, dt):

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
    theta = soil.theta(dt)
    T_soil = np.array(soil.T(dt)) + Tzero_sli
    rH = None
    psi = soil.head(dt)  # sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c, soil, dt):
    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Surface variables
    T_surface = soil.T(dt)[0] + Tzero_sli  ## TODO: recheck the surface temperature value ?????
    q_ev = soil.qevap(dt)

    ql_surface = 0
    qv_surface = 0
    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

    # soil fluxes
    ql = soil.q_l(dt)[1:]  # p.array()
    qv = soil.q_v(dt)[1:]  # np.array()
    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

    # boundary storage
    c.update_connection_to_aquifer()


def run_iso(p, soil, **kwargs):

    solutes = ["2H", "18O"]
    c = p.get_cells()[0]  # get current cell of project

    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}

    time_steps = len(soil.get_head())
    for dt in range(time_steps-1):
        print(dt / 24, ' days')
        print(c.conc_2H_delta[0])

        c_iso["2H"].append(c.conc_2H), c_iso["18O"].append(c.conc_18O)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        # update storage states and boundaries to current time
        update_storages(c, soil, dt), update_boundaries(c, soil, dt)
        for solute in solutes:

            delta_t = 60*60
            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=None, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


def iso_setup(soil, testcase=1, **kwargs):

    # Define boundary storages
    atm = _atm(testcase)  # atmosphere
    aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # define aquifer as boundary isotope storage

    p = iso_project()  # create a project
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell

    c = p.get_cells()[0]  # get current cell
    _layers(c, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation()
    c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

    return run_iso(p, soil, **kwargs)


def run_testcases(soil, test_cases):

    delta = {}
    for Testcase in test_cases:
        print('Testcase:', Testcase)

        ignore = iso_delta.test_case_args(Testcase)

        delta[Testcase] = {}
        d = iso_setup(soil, testcase=Testcase, **ignore)

        # delta at the end of simulation for each test cases
        delta[Testcase]["2H"] = d["2H"][-1]
        delta[Testcase]["18O"] = d["18O"][-1]

    return delta


path = r"D:\Isotope transport\soil water models\VaporModel_Example_B\VaporModel_Example" \
       r"\VaporIllustrationEx2\03_VaryEvpConsTemp"
_soil = _2D_Soil(path)  # reading input from Soil2D model

tests = [1]
diso = run_testcases(_soil, tests)

depth = np.arange(0, 50.5, 2.5)
#plot_delta(diso, depth,'2H', tests)

