'''
Created on 03.12.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

from MOIST import moist
from src import *

import matplotlib.pyplot as plt


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


def _layers(c, m, testcase):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    lower_boundaries = np.arange(0.1, 1.1, 0.1)
    theta, tmp, head = m.theta(0), m.T(0), m.head(0)  # 0 index for first time step values
    theta_sat, theta_r, tortuosity = 0.35, 0.01, 0.67,

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    init_c_iso_2H = [iso_delta.delta_to_concentration(ali_2H, '2H')] * len(lower_boundaries)
    init_c_iso_18O = [iso_delta.delta_to_concentration(ali_18O, '18O')] * len(lower_boundaries)

    id = 0
    upper_boundary = 0
    for lower_boundary, c2H, c18O, th, T, h in zip(lower_boundaries, init_c_iso_2H, init_c_iso_18O, theta, tmp, head):

        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=upper_boundary,
                                                lower_boundary=lower_boundary,
                                                conc_iso_liquid={"2H": c2H, "18O": c18O},
                                                theta=0.35,
                                                theta_0=theta_r,
                                                theta_sat=theta_sat,
                                                tortuosity=tortuosity,
                                                T=T + Tzero_sli,
                                                rH=None,
                                                psi=h
                                                )
        id += 1
        upper_boundary = lower_boundary
        c.add_layer(new_layer)

    return c


def update_storages(c, m, dt):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Need not update atmosphere
    dt = dt
    # Update layers
    theta = m.theta(dt)
    T_soil = np.array(m.T(dt)) + Tzero_sli
    rH = None
    psi = m.head(dt)  # sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)

    # boundary storage
    h0 = 0  # sli.pond_height(dt)
    c.update_pond(pond_height=h0)

    # c_aquifer = sli.cali(dt)
    c.update_aquifer(c_iso={"2H": 0.0, "18O": 0.0})


def update_boundaries(c, m, dt):
    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Surface variables
    T_surface = m.Ts(dt) + Tzero_sli  ## TODO: recheck the surface temperature value ?????
    q_ev = m.qev(dt)
    ql_surface = m.qls(dt)
    qv_surface = m.qvs(dt)
    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

    # soil fluxes
    ql = m.ql(dt)[1:]  # p.array()
    qv = m.qv(dt)[1:]  # np.array()
    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

    # boundary storage
    c.update_connection_to_aquifer()

    # transpiration
    ql_trans = [0.0] * len(ql)
    c.update_transpiration(q_trans=ql_trans)

    # boundary storage
    c.update_connection_to_aquifer()


def run_iso(p, m, **kwargs):

    solutes = ["2H", "18O"]
    c = p.get_cells()[0]  # get current cell of project

    c_iso_delta = {'2H': [], '18O': []}
    T, count = 0, 0
    for i, dt in enumerate(m.dt):

        delta_t = float(dt)
        T += delta_t
        count += 1

        print('count: ', count, float(T) / 86400, ' days')
        print(c.conc_2H_delta[0], c.conc_18O_delta[0])

        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        # update storage states and boundaries to current time
        update_storages(c, m, i), update_boundaries(c, m, i)
        for solute in solutes:

            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=None, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


def iso_setup(m, testcase=1, **kwargs):

    # Define boundary storages
    atm = _atm(testcase)  # atmosphere
    pd = iso_storages.iso_pond(pond_height=0)  # define pond
    aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # define aquifer as boundary isotope storage

    p = iso_project()  # create a project
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell

    c = p.get_cells()[0]  # get current cell
    _layers(c, m, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation(), c.add_transpiration(), c.add_surface_runoff()
    c.add_precipitation(), c.add_pond(pd), c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

    return run_iso(p, m, **kwargs)


def run_testcases(m, test_cases):

    delta = {}
    for Testcase in test_cases:
        print('Testcase:', Testcase)

        ignore = iso_delta.test_case_args(Testcase)

        delta[Testcase] = {}
        d = iso_setup(m, testcase=Testcase, **ignore)

        # delta at the end of simulation for each test cases
        delta[Testcase]["2H"] = d["2H"][-1]
        delta[Testcase]["18O"] = d["18O"][-1]

    return delta


path = r"D:\Isotope transport\soil water models\MOIST\8397416\thoritical_test\output"
moist = moist(path)
tests = [1]

diso = run_testcases(moist, tests)

#######################
m = moist
t = len(m.dt)

qs = [m.qev(tm) for tm in range(0,t)]
plt.plot([3950]*9, np.arange(0, 2.25e-7, 0.25e-7))
plt.plot(qs)
plt.yscale('log')
plt.show()


ths0 = [m.theta(tm)[0] for tm in range(0,t)]
ths1 = [m.theta(tm)[1] for tm in range(0,t)]
ths2 = [m.theta(tm)[2] for tm in range(0,t)]
d = np.arange(0.05, 0.35, 0.03)
plt.plot([3950]*10, d)
plt.plot(ths0, label='1')
plt.plot(ths1, label='2')
plt.plot(ths2, label='3')
plt.legend()
plt.show()

ts0 = [m.T(tm)[0] for tm in range(0,t)]
ts1 = [m.T(tm)[1] for tm in range(0,t)]
ts2 = [m.T(tm)[2] for tm in range(0,t)]
d = np.arange(27, 30, 0.3)
plt.plot([3950]*10, d)
plt.plot(ts0, label='1')
plt.plot(ts1, label='2')
plt.plot(ts2, label='3')
plt.legend()
plt.show()

ql0 = [m.ql(tm)[0] for tm in range(0,t)]
ql1 = [m.ql(tm)[1] for tm in range(0,t)]
ql2 = [m.ql(tm)[2] for tm in range(0,t)]

d = np.arange(-1.5e-7, 1e-7, 0.2e-7)
plt.plot([3950]*13, d)
plt.plot(ql0, label='1')
plt.plot(ql1, label='2')
plt.plot(ql2, label='3')
plt.legend()
plt.show()

qv0 = [m.qv(tm)[0] for tm in range(0,t)]
qv1 = [m.qv(tm)[1] for tm in range(0,t)]
qv2 = [m.qv(tm)[2] for tm in range(0,t)]

d = np.arange(-1.5e-9, 1e-9, 0.2e-9)
plt.plot([3950]*13, d)
plt.plot(qv0, label='1')
plt.plot(qv1, label='2')
plt.plot(qv2, label='3')
plt.legend()
plt.show()
