'''
Created on 03.12.2024
@author: poudel-b
'''
import numpy as np

# -*- coding: utf-8 -*-

from MOIST import moist, moist_v1
from src import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm


def visualize(p, delta, Isotopologue='2H'):

    depth = [-l.center.z for l in p.get_cells()[0].layers]
    plt.figure(figsize=(7, 10))
    for tests in delta.keys():
        plt.plot(delta[tests][Isotopologue][-1][:10], depth[:10], label='test_{}'.format(tests))

    plt.xlabel(r'$\delta$ {}'.format(Isotopologue))
    plt.ylabel('depth $[m]$')
    plt.title('Initial testcases after 250 days')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()


def _atm(testcase):

    Tzero_sli = 273.17  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = 30
    Rh_atm = 0.2
    wind_speed = 0.00001  # m/s

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm + Tzero_sli, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=2.5,  # canopy height [m] (e.g. 40)
                         d0=0,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.001,
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=0,  # leaf area index (e.g. 2.0)
                         extku=0)

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
    # dt = dt
    # Update layers
    theta = m.theta(dt)
    T_soil = np.array(m.T(dt)) + Tzero_sli
    rH = None
    psi = m.head(dt)  # sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)

    # boundary storage
    # h0 = 0  # sli.pond_height(dt)
    # c.update_pond(pond_height=h0)

    # c_aquifer = sli.cali(dt)
    # c.update_aquifer(c_iso={"2H": 0.0, "18O": 0.0})


def update_boundaries(c, m, dt):
    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.17  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Surface variables
    T_surface = 303.17  # m.Ts(dt) + Tzero_sli  ## TODO: recheck the surface temperature value ?????
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
    # c.update_connection_to_aquifer()

    # transpiration
    # ql_trans = [0.0] * len(ql)
    # c.update_transpiration(q_trans=ql_trans)

    # boundary storage
    # c.update_connection_to_aquifer()


def run_iso(p, m, **kwargs):

    solutes = ["2H"]
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

    return p, c_iso_delta


def iso_setup(m, testcase=1, **kwargs):

    # Define boundary storages
    atm = _atm(testcase)  # atmosphere
    # pd = iso_storages.iso_pond(pond_height=0)  # define pond
    # aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # define aquifer as boundary isotope storage

    p = iso_project()  # create a project
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell

    c = p.get_cells()[0]  # get current cell
    _layers(c, m, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation()  # , c.add_transpiration(), c.add_surface_runoff()
    #c.add_precipitation(), c.add_pond(pd), c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

    return run_iso(p, m, **kwargs)


def run_testcases(m, test_cases):

    delta = {}
    for Testcase in test_cases:
        print('Testcase:', Testcase)

        ignore = iso_delta.test_case_args(Testcase)

        delta[Testcase] = {}
        p, d = iso_setup(m, testcase=Testcase, **ignore)

        # delta at the end of simulation for each test cases
        delta[Testcase]["2H"] = d["2H"][-1]
        delta[Testcase]["18O"] = d["18O"][-1]

    return p, delta


path = r"D:\Isotope transport\soil water models\MOIST\8397416\thoritical_test\output"
m = moist_v1(path)
tests = [1]

p_iso, diso = run_testcases(m, tests)

# Iso_delta
visualize(p_iso, diso, Isotopologue='2H')
visualize(p_iso, diso, Isotopologue='18O')


#######################
t = len(m.dt)

f = 86400 * 1000  # ms-1 to mm day-1

delta_t = m.dt
cum_t = np.cumsum(delta_t)[:2340]
days = cum_t / 86400
t_steps = 2340  # len(m.dt)

plt_days = [50, 100, 150, 200, 250]
plt_steps = [np.argmin(np.abs(days - target)) for target in plt_days]

depth = -np.arange(0.1, 1.1, 0.1)
qev = np.array([m.qev(t) * f for t in range(t_steps)])

theta = np.array([m.theta(t) for t in plt_steps])

ql = np.array([m.ql(t)[1:] for t in plt_steps])
qv = np.array([m.qv(t)[1:] for t in plt_steps])
T = np.array([m.T(t) for t in plt_steps])
pot = np.array([m.head(t) for t in plt_steps])

################## theta ###################
for th, l in zip(theta, plt_days):
    plt.plot(th, depth, label=str(l) + ' days')

plt.title('theta')
plt.xlabel('[mm3 / mm3 ]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()

############### ql profile #######################
colors = cm.viridis(np.linspace(0, 1, len(ql)))
for q_l, q_v, l, c in zip(ql, qv, plt_days, colors):
    plt.plot(q_l * f, depth, label=str(l) + ' days', color=c)
    plt.plot(q_v * f, depth, color=c, linestyle='--')

plt.title('ql')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


######################## qv profile ##################
for q_v, l in zip(qv, plt_days):
    plt.plot(q_v * f, depth, label=str(l) + ' days')

plt.title('qv')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()

################# coupled profile ###############
plt.plot((ql[-1]) * f, depth, label='ql')
plt.plot(qv[-1] * f, depth, label='qv')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


################### Temp profile ###########################
for Tmp, l in zip(T, plt_days):
    plt.plot(Tmp, depth, label=str(l) + ' days')

plt.title('Temp')
plt.xlabel('tmp [C]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


###############  potential profile ###############
for p, l in zip(pot, plt_days):
    plt.plot(p, depth, label=str(l) + ' days')

plt.title('matric potential')
plt.xlabel('potential [m]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


####################################################
############# temporal surface fluxes #############
q_ev = np.array([m.qev(t) for t in range(t_steps)]) * f
ql0 = np.array([m.qls(t) for t in range(t_steps)]) * f
qv0 = np.array([m.qvs(t) for t in range(t_steps)]) * f
q_tot = np.array(ql0) + np.array(qv0)

plt.plot(days, q_ev, label='qev')
plt.plot(days, ql0, label='ql0')
plt.plot(days, qv0, label='qv0')
plt.plot(days, q_tot, label='total_fluxes')

plt.title('surface fluxes')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()


############ ql qv Temporal #########################

qv1 = np.array([m.qv(t)[0] for t in range(t_steps)]) * f
qv2 = np.array([m.qv(t)[1] for t in range(t_steps)]) * f
qv3 = np.array([m.qv(t)[2] for t in range(t_steps)]) * f

ql1 = np.array([m.ql(t)[0] for t in range(t_steps)]) * f
ql2 = np.array([m.ql(t)[1] for t in range(t_steps)]) * f
ql3 = np.array([m.ql(t)[2] for t in range(t_steps)]) * f

plt.plot(days, qv1, label='qv_1')
plt.plot(days, qv2, label='qv_2')
plt.plot(days, qv3, label='qv_3')
plt.plot(days, ql1, label='ql_1')
plt.plot(days, ql2, label='ql_2')
plt.plot(days, ql3, label='ql_3')
plt.title('SLI_ql_qv')
plt.xlabel('[days]')
#plt.xscale('log')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
#plt.ylim(-2.5, 0.5)
plt.show()


################ theta temporal ############################

th1 = np.array([m.theta(t)[0] for t in range(t_steps)])
th2 = np.array([m.theta(t)[1] for t in range(t_steps)])
th3 = np.array([m.theta(t)[2] for t in range(t_steps)])

plt.plot(days, th1, label='th_l1')
plt.plot(days, th2, label='th_l2')
plt.plot(days, th3, label='th_l3')

#plt.plot(days, qev, label='qev')
plt.title('theta')
plt.xlabel('[days]')
#plt.xscale('log')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
#plt.ylim(-2.5, 0.5)
plt.show()


############### temporal  potential #####
pot1 = np.array([m.head(t)[0] for t in range(t_steps)])
pot2 = np.array([m.head(t)[1] for t in range(t_steps)])
pot3 = np.array([m.head(t)[2] for t in range(t_steps)])

plt.plot(days, pot1, label='potential_l1')
plt.plot(days, pot2, label='potential_l2')
plt.plot(days, pot3, label='potential_l3')

#plt.plot(days, qev, label='qev')
plt.title('matric potential')
plt.xlabel('[days]')
#plt.xscale('log')
plt.ylabel('[m]')
plt.grid()
plt.legend()
#plt.ylim(-2.5, 0.5)
plt.show()



