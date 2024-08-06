'''
Created on 03.06.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

import os
import numpy as np

import Sli
import iso_project
import iso_storages
import iso_delta

import matplotlib.pyplot as plt


def visualize(delta, sli, Isotopologue):

    depth = np.insert(np.cumsum(np.array(sli.dx(0))), 0, 0)
    centers = -(depth[:-1] + depth[1:]) / 2.0
    d = centers.tolist()[:10]

    for case in delta.keys():

        plt.plot(delta[case][Isotopologue][:10], d, label='testcase_{}'.format(case))
        #plt.plot(delta[2][:10], d, label='testcase_2')

    #plt.xlim([-15, 15])
    plt.xlabel('delta_{}'.format(Isotopologue))
    plt.ylabel('depth [m]')
    plt.title('initial testcases after 250 days')
    plt.legend()
    plt.gca().set_aspect(aspect=150)
    plt.show()


def test_case(testcase=1):

    ignore = {'ignoredvi': True, 'ignoredli': True, 'ignorealphai': True,
              'ignorealphaik': True}  # Testcases: Mathieu and Bariac (1996)

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


def get_sli():
    pth = os.getcwd()
    path = os.path.abspath(os.path.join(pth, "..", ".."))

    # imports all the variable files /variables folder: testcase-1, sig=1
    sli = Sli.SlI(path + '\_sli_\sli_label3\iso_variables_1')

    return sli


def _layers(c, sli):

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

    init_c_iso_2H = [iso_delta.delta_to_concentration(-65, '2H')] * len(lower_boundaries)  # 0.15367056287906838
    init_c_iso_18O = [iso_delta.delta_to_concentration(-8, '18O')] * len(lower_boundaries)  # sli.civa(dt) / sli.cva(dt)

    id = 0
    upper_boundary = 0
    for lower_boundary, c2H, c18O, th, tr, tsat, tor, T, r_H, PSI in \
            zip(lower_boundaries, init_c_iso_2H, init_c_iso_18O, theta, theta_r, theta_sat, tortuosity, T_soil, rH, psi):
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
        c.add_layer(new_layer)

    return c


def _atm(sli, testcase):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    dt = 0  # initial states
    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = sli.Tatm(dt) + Tzero_sli
    Rh_atm = sli.R_humidity_atm(dt)
    wind_speed = sli.wind_speed(dt)

    if testcase == 1 or testcase == 3:
        c_iso_2H = iso_delta.delta_to_concentration(-65, '2H')  # 0.15367056287906838
        c_iso_18O = iso_delta.delta_to_concentration(-8, '18O')  # sli.civa(dt) / sli.cva(dt)
    elif testcase in [2, 4, 5, 6]:
        c_iso_2H = iso_delta.delta_to_concentration(-112, '2H')  # 0.15367056287906838
        c_iso_18O = iso_delta.delta_to_concentration(-15, '18O')  # sli.civa(dt) / sli.cva(dt)
    else:
        raise NotImplementedError

    atm = iso_storages.iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
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


def update_storages(c, sli, dt):

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
    theta = sli.theta(dt)
    T_soil = np.array(sli.T_soil0(dt)) + Tzero_sli
    rH = sli.R_humidity(dt)
    psi = sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)

    # boundary storage
    h0 = sli.pond_height(dt)
    c.update_pond(pond_height=h0)

    c_aquifer = sli.cali(dt)
    c.update_aquifer(c_iso={"2H": 1.0, "18O": c_aquifer})


def update_boundaries(c, sli, dt):
    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # precipitation
    q_prec = sli.qprec(dt)
    c_prec = sli.cprec(dt)
    c.update_precipitation(q_prec=q_prec, c_prec={"2H": 1.0, "18O": c_prec})

    # runoff
    qrunoff = sli.qrunoff(dt)
    c.update_runoff(q_runoff=qrunoff)

    # Surface variables
    T_surface = sli.Ts(dt) + Tzero_sli
    q_ev = sli.qevap(dt)
    ql_surface = sli.ql0(dt)
    qv_surface = sli.qv0(dt)
    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

    # soil fluxes
    ql = sli.qlsig(dt)[1:]
    qv = sli.qvsig(dt)[1:]

    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

    # transpiration
    ql_trans = sli.qex(dt)
    c.update_transpiration(q_trans=ql_trans)

    # boundary storage
    c.update_connection_to_aquifer()


def run_iso(p, sli, **kwargs):

    c = p.get_cells()[0]  # get current cell of project

    c_iso = {'2H': [c.conc_2H], '18O': [c.conc_18O]}
    c_iso_delta = {'2H': [c.conc_2H_delta], '18O': [c.conc_18O_delta]}

    solutes = ["2H", "18O"]

    for dt in range(1, 5000):  #len(sli.get_in_soil())):  # starting from next time step (n+1)
        print(dt)
        # update storage states and boundaries to current time
        update_storages(c, sli, dt), update_boundaries(c, sli, dt)

        for solute in solutes:

            delta_t = sli.dt(dt)
            dc = p.run(Isotopologue=solute, delta_time=delta_t, **kwargs)

            current_conc = c.get_conc_layers(Isotopologue=solute)
            c_t = list(np.array(current_conc) + np.array(dc))
            delta = [iso_delta.concentration_to_delta(c_iso, solute) for c_iso in c_t]

            c_iso[solute].append(c_t), c_iso_delta[solute].append(delta)

            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)

    return c_iso, c_iso_delta


def iso_setup(sli, testcase=None, **kwargs):

    p = iso_project.iso_project()  # create a project

    atm = _atm(sli, testcase)  # atmosphere
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell
    c = p.get_cells()[0]
    _layers(c, sli)  # add soil layers to the new cell
    c.install_connections()  # flux connections between the layers

    pd = iso_storages.iso_pond(pond_height=0)
    c.add_pond(pd)

    # boundary connections
    c.add_evaporation(), c.add_transpiration(), c.add_surface_runoff(), c.add_precipitation()

    aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # aquifer as boundary isotope storage
    c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

    return run_iso(p, sli, **kwargs)


def run_testcases(test_cases, sli):

    delta = {}
    d_solute = {}
    for Testcase in test_cases:

        print('Testcase:', Testcase)
        cases = test_case(testcase=Testcase)

        c, d = iso_setup(sli, testcase=Testcase, **cases)

        # delta at the end of simulation for each test cases
        d_solute["2H"] = d["2H"][-1]
        d_solute["18O"] = d["18O"][-1]
        delta[Testcase] = d_solute

    return delta


def moisture(sli):

    n_timesteps = len(sli.get_scaler())

    ql = np.array(sli.qlsig(n_timesteps - 1)) * 1000 * 86400
    qv = np.array(sli.qvsig(n_timesteps - 1)) * 1000 * 86400
    theta = sli.theta(n_timesteps - 1)

    depth = np.insert(np.cumsum(np.array(sli.dx(0))), 0, 0)
    centers = -(depth[:-1] + depth[1:]) / 2.0
    d = centers.tolist()[:10]

    plt.plot(theta[:10], d)
    plt.xlabel('m3 m-3')
    plt.ylabel('depth m')
    plt.title('moisture at 250-day')
    plt.show()

    plt.plot(ql[:10], d, label='ql')
    plt.plot(qv[:10], d, label='qv')
    plt.xlabel('mm per day')
    plt.ylabel('depth m')
    plt.title('moisture_flux at 250-day')
    plt.legend()
    plt.show()


slI = get_sli()
ignore = test_case(testcase=1)
#delta = run_testcases([1], slI)

conc, d = iso_setup(sli=slI, testcase=1, **ignore)
#visualize(delta=delta, sli=slI, Isotopologue="2H")
#visualize(delta=delta, sli=slI, Isotopologue="18O")

#moisture(slI)



