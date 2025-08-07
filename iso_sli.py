'''
Created on 03.06.2024
@author: poudel-b
'''

import copy
# -*- coding: utf-8 -*-

import os
import numpy as np

import Sli
from src import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm


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

    plt.figure(figsize=(6, 10))
    plt.plot(ql[:10], d, label='ql')
    plt.plot(qv[:10], d, label='qv')
    plt.xlabel('mm per day')
    plt.ylabel('depth m')
    plt.title('moisture_flux at 250-day')
    plt.legend()
    plt.show()


def visualize(delta, sli, Isotopologue):
    depth = np.insert(np.cumsum(np.array(sli.dx(0))), 0, 0)
    centers = -(depth[:-1] + depth[1:]) / 2.0
    d = centers.tolist()
    plt.figure(figsize=(6, 8,))

    if Isotopologue == '2H':
        i, j = 2, 'H'
    else:
        i, j = 18, 'O'

    for case in delta.keys():
        plt.plot(delta[case][Isotopologue][-1][:10], d[:10], label=r'testcase_{}'.format(case))
        # plt.plot(delta[2][:10], d, label='testcase_2')

    plt.xlabel(r'$\delta^{{{}}}{}$ ‰ '.format(i, j), fontsize=15)
    plt.ylabel(r'Depth $[m]$', fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=10)
    # plt.title('Initial testcases after 250 days')
    legend = plt.legend(loc='lower center')
    legend.get_frame().set_facecolor('none')  # Set the legend background to be transparent
    legend.get_frame().set_edgecolor('none')  # Remove the legend border if needed

    # Removing grid lines
    plt.grid(False)

    # Removing spines (borders around the plot)
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    # Save the plot with a transparent background
    plt.savefig(r'D:\Isotope transport\Scripts\output\transparent_plot.png', format='png', transparent=True)
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
    elif testcase in [6, 7, 8]:

        ignore['ignorealphai'] = False
        ignore['ignoredli'] = False
        ignore['ignoredvi'] = False
        ignore['ignorealphaik'] = False

        return ignore

    else:
        raise NotImplementedError


def get_sli(testcase=1):
    pth = os.getcwd()
    path = os.path.abspath(os.path.join(pth, "", ".."))

    # imports all the variable files /variables folder: testcase-1, sig=1
    sli = Sli.SlI(path + '\\_sli_\\sli_label3\\iso_variables_{}'.format(1))

    return sli


def _layers(c, sli, testcase):
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
                         extku=1.5)  # 1.5)

    return atm


def update_S(c, dy):

    saturation = np.array([l.S for l in c.layers])
    thre = np.array([(l.theta_sat - l.theta_0) for l in c.layers])
    thetar = np.array([l.theta_0 for l in c.layers])

    updated_saturation = saturation + dy
    updated_theta = updated_saturation * thre + thetar

    return updated_theta


def atm_fluxes(top_node, Ta=30.0):
    Tzero = 273.16000366210938

    return flux_atmosphere(top_layer=top_node, Ta=Ta + Tzero)


def check_timesteps(sli, dt, tsteps):
    t_step = sli.time_step(dt)

    tsteps.nsat_0 = tsteps.nsat  # set the previous nsat values
    tsteps.set_nsat()  # update  current nsat

    tsteps.T_step0 = tsteps.T_step  # set previous time step
    tsteps.T_step = t_step  # update current time step

    if tsteps.T_step == tsteps.T_step0:  # same current time step

        return False
    else:

        return tsteps.repeat()


def water_heat_solve(c, sli, dt, tsteps):

    dy, dT = solve_coupled(c, sli, dt, repeat=False)

    repeat = check_timesteps(sli, dt, tsteps)
    if repeat:
        dy, dT = solve_coupled(c, sli, dt, repeat=True, dy=dy)

    return dy, dT


def solve_coupled(c, sli, dt, repeat=False, dy=None):

    d_t = sli.dt(dt)
    qev = sli.qevap(dt)  # atm.q_evaporation()

    top_node = c.top_layer
    atm = flux_atmosphere(atmosphere=c.atmosphere, top_layer=top_node)
    atm.qevap = qev

    #print(atm.lE0() == atm.ev_pot())

    if repeat:
        atm.repeat = repeat

    yb, Tb = atm.qyb_qTb()
    q, qh = [-qev], [atm.G0()]
    qya, qyb, qTa, qTb = [0], [yb], [0], [Tb]
    qhya, qhyb, qhTa, qhTb = [0], [0], [0], [0]

    for cn in c.connections_v_adv:

        l_node, r_node = cn.left_node, cn.right_node

        hq = heat_flux(left_node=l_node, right_node=r_node, top_layer=top_node)
        if repeat:
            hq.repeat = repeat
            hq.dy = dy
            hq.nodes = c.layers

        q.append(hq.q()), qh.append(hq.qh())

        qya.append(hq.qya()), qyb.append(hq.qyb())
        qTa.append(hq.qTa()), qTb.append(hq.qTb())

        qhya.append(hq.qhya()), qhyb.append(hq.qhyb())
        qhTa.append(hq.qhTa()), qhTb.append(hq.qhTb())

    hs = heat_solve()
    del_cc, del_cch = [], []
    del_ddh, del_ggh = [], []

    for l in c.layers:
        del_cc.append(hs.delta_cc(node=l, dt=d_t))
        del_cch.append(hs.delta_cch(node=l, dt=d_t))
        del_ddh.append(hs.delta_ddh(node=l, dt=d_t))
        del_ggh.append(hs.delta_ggh(node=l, dt=d_t))

    # boundary zero fux
    q.append(0), qh.append(0)

    qya.append(0), qyb.append(0), qhya.append(0), qhyb.append(0)
    qTa.append(0), qTb.append(0), qhTa.append(0), qhTb.append(0)

    aa, aah = np.array(qya[:-1]), np.array(qhya[:-1])
    bb, bbh = np.array(qTa[:-1]), np.array(qhTa[:-1])

    ee, eeh = - np.array(qyb[1:]), - np.array(qhyb[1:])
    ff, ffh = - np.array(qTb[1:]), - np.array(qhTb[1:])

    cc = np.array(qyb[:-1]) - np.array(qya[1:]) - np.array(del_cc)
    cch = np.array(qhyb[:-1]) - np.array(qhya[1:]) + np.array(del_cch)

    dd = np.array(qTb[:-1]) - np.array(qTa[1:])
    ddh = np.array(qhTb[:-1]) - np.array(qhTa[1:]) - np.array(del_ddh)

    gg = - (np.array(q[:-1]) - np.array(q[1:]))
    ggh = - (np.array(qh[:-1]) - np.array(qh[1:])) + np.array(del_ggh)

    x = [aa[1:], bb[1:], cc, dd, ee[:-1], ff[:-1], gg]
    X = [aah[1:], bbh[1:], cch, ddh, eeh[:-1], ffh[:-1], ggh]

    dy, dT = hs.solve_sparse(x, X)

    return dy, dT


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

    #dy, dT = water_heat_solve(c, sli, dt, tsteps)
    #Tsoil = np.array([l.T for l in c.layers]) + dT

    #theta0 = update_S(c, dy)
    # print(sli.qevap(dt))

    #print('dy')
    #print(sli.dy(dt)[1:])
    #print(dy.tolist())

    #print('dT')
    #print(sli.deltaT(dt)[1:])
    #print(dT.tolist())

    # Update layers
    theta = sli.theta(dt)
    T_soil = np.array(sli.T_soil0(dt)) + Tzero_sli
    rH = sli.R_humidity(dt)
    psi = sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)

    # dy, dT = solve_coupled(c, sli, dt)

    # boundary storage
    h0 = sli.pond_height(dt)
    c.update_pond(pond_height=h0)

    # c_aquifer = sli.cali(dt)
    c.update_aquifer(c_iso={"2H": 0.0, "18O": 0.0})


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
    T_surface = 30 + Tzero_sli # sli.Ts(dt) + Tzero_sli
    q_ev = sli.qevapsig(dt)
    ql_surface = sli.ql0(dt)
    qv_surface = sli.qv0(dt)
    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

    # soil fluxes
    ql = sli.qlsig(dt).tolist()     # top flux is not included in ql
    qv = sli.qvsig(dt)[1:].tolist()    # since list includes top vapor flux

    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

    # transpiration
    ql_trans = sli.qex(dt)
    c.update_transpiration(q_trans=ql_trans)

    # boundary storage
    c.update_connection_to_aquifer()


def run_iso(p, sli, **kwargs):

    solutes = ["2H", "18O"]
    c = p.get_cells()[0]  # get current cell of project

    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}
    for dt in range(1, len(sli.get_in_soil()) - 1):  # starting from next time step (n+1)
        # for dt in range(1, 637):

        print(dt)
        print(c.conc_2H_delta)

        c_iso["2H"].append(c.conc_2H), c_iso["18O"].append(c.conc_18O)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        # update storage states and boundaries to current time
        update_boundaries(c, sli, dt), update_storages(c, sli, dt)
        for solute in solutes:
            delta_t = sli.dt(dt)
            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=1e-10, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


def iso_setup(sli, testcase=1, **kwargs):
    # Define boundary storages
    atm = _atm(sli, testcase)  # atmosphere
    pd = iso_storages.iso_pond(pond_height=0)  # define pond
    aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # define aquifer as boundary isotope storage

    p = iso_project()  # create a project
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell

    c = p.get_cells()[0]  # get current cell
    _layers(c, sli, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation(), c.add_transpiration(), c.add_surface_runoff()
    c.add_precipitation(), c.add_pond(pd), c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

    return run_iso(p, sli, **kwargs)


def run_testcases(test_cases):

    delta = {}
    for Testcase in test_cases:
        print('Testcase:', Testcase)
        cases = test_case(testcase=Testcase)

        delta[Testcase] = {}

        sli = get_sli(testcase=Testcase)
        d = iso_setup(sli, testcase=Testcase, **cases)

        # delta at the end of simulation for each test cases
        delta[Testcase]["2H"] = d["2H"]
        delta[Testcase]["18O"] = d["18O"]

    return delta


def slope(testcase, delta):
    d_2H_ali, d_2H_atm = iso_delta.delta_testcases(Isotopologue="2H", testcase=testcase)
    d_18O_ali, d_18O_atm = iso_delta.delta_testcases(Isotopologue="18O", testcase=testcase)

    d_2H_sim, d_18O_sim = delta[testcase]["2H"][0], delta[testcase]["18O"][0]

    slope_initial = (d_2H_ali - d_2H_atm) / (d_18O_ali - d_18O_atm)
    slope_sim = (d_2H_sim - d_2H_atm) / (d_18O_sim - d_18O_atm)

    return slope_initial, slope_sim


def enrichment_max(Isotopologue, testcase, delta):
    delta_soil, delta_vap = iso_delta.delta_testcases(Isotopologue=Isotopologue, testcase=testcase)

    return max(abs(np.array(delta[testcase][Isotopologue]) - delta_soil))


sli = get_sli(testcase=1)
delta = run_testcases([1])

visualize(delta=delta, sli=sli, Isotopologue="2H")
visualize(delta=delta, sli=sli, Isotopologue="18O")

#ignore = test_case(testcase=1)
# d = iso_setup(sli=slI, testcase=2, **ignore)
# moisture(sli)

f = 86400 * 1000  # ms-1 to mm day-1

delta_t = np.array([sli.dt(t) for t in range(len(sli.get_scaler())-1)])
cum_t = np.cumsum(delta_t)
days = cum_t / 86400

plt_days = [50, 100, 150, 200, 250]
plt_steps = [np.argmin(np.abs(days - target)) - 1 for target in plt_days]

depth = - np.cumsum(sli.dx(0))[:10]
qev = [sli.qevapsig(t) * f for t in range(len(sli.get_scaler())-1)]
epot = [sli.E_pot(t) for t in range(len(sli.get_scaler())-1)]

theta = np.array([sli.theta(t)[:10] for t in plt_steps])
ql = np.array([sli.qlsig(t)[:10] for t in plt_steps])
qv = np.array([sli.qvsig(t)[1:11] for t in plt_steps])
T = np.array([sli.T_soil0(t)[:10] for t in plt_steps])
pot = np.array([sli.matric_pot(t)[:10] for t in plt_steps])

#################### delta profile ##########################
delta_2H, delta_18O = delta[1]['2H'], delta[1]['18O']
_delta_18O = [delta_2H[t] for t in plt_steps]
for d, l in zip(_delta_18O, plt_days):
    plt.plot(d[:10], depth, label=str(l) + ' days')

plt.title('delta')
plt.xlabel('delta')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()
#############################################################

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


############### Temporal evaporation potential###############
def lambdav(T):
    return 1.91846e6 * (T / (T - 33.91)) ** 2

q_ev = [sli.qevap(t) for t in range(len(sli.get_scaler())-1)]
Tsoil_l1 = np.array([sli.T_soil0(t)[0] + 273.15 for t in range(len(sli.get_scaler())-1)])
lamda = [lambdav(T) for T in Tsoil_l1]
lE0 = np.array(q_ev) * np.array(lamda) * 1000

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
l1, = ax1.plot(days, qev, '-b', label='qev')
# l2, = ax2.plot(days, np.array(epot) / (np.array(lamda)) * 86400, '-g', label='epot')
l2, = ax2.plot(days, epot, '-g', label='epot')
l3, = ax2.plot(days, lE0, '-y', label='E_evap')

ax1.set_xlabel('[days]')
ax1.set_ylabel('[mm per day]')

labels = [l1, l2, l3]
ax1.legend(labels, [l.get_label() for l in labels])
plt.title('qevap vs E_pot')
plt.grid()
plt.show()

############# temporal surface fluxes #############
q_ev = [sli.qevapsig(t) * f for t in range(len(sli.get_scaler())-1)]
ql0 = [-sli.ql0(t) * f for t in range(len(sli.get_scaler())-1)]
qv0 = [-sli.qv0(t) * f for t in range(len(sli.get_scaler())-1)]
q_tot = np.array(ql0) + np.array(qv0)

plt.plot(days[:1000], q_ev[:1000], label='qev')
plt.plot(days[:1000], ql0[:1000], label='ql0')
plt.plot(days[:1000], qv0[:1000], label='qv0')
plt.plot(days[:1000], q_tot[:1000], label='total_fluxes')

plt.title('surface fluxes')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()

############ql qv Temporal #########################

qv1 = [sli.qvsig(t)[1] * f for t in range(len(sli.get_scaler()))]
qv2 = [sli.qvsig(t)[2] * f for t in range(len(sli.get_scaler()))]
qv3 = [sli.qvsig(t)[3] * f for t in range(len(sli.get_scaler()))]

ql1 = [sli.qlsig(t)[0] * f for t in range(len(sli.get_scaler()))]
ql2 = [sli.qlsig(t)[1] * f for t in range(len(sli.get_scaler()))]
ql3 = [sli.qlsig(t)[2] * f for t in range(len(sli.get_scaler()))]

plt.plot(days, qv1[1:], label='qv_1')
plt.plot(days, qv2[1:], label='qv_2')
plt.plot(days, qv3[1:], label='qv_3')
#plt.plot(days, ql1[1:], label='ql_1')
#plt.plot(days, ql2[1:], label='ql_2')
#plt.plot(days, ql3[1:], label='ql_3')
plt.title('SLI_ql_qv')
plt.xlabel('[days]')
#plt.xscale('log')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
#plt.ylim(-2.5, 0.5)
plt.show()


################ theta temporal ############################

th1 = [sli.theta(t)[0] for t in range(len(sli.get_scaler()))]
th2 = [sli.theta(t)[1] for t in range(len(sli.get_scaler()))]
th3 = [sli.theta(t)[2] for t in range(len(sli.get_scaler()))]

plt.plot(days, th1[1:], label='th_l1')
plt.plot(days, th2[1:], label='th_l2')
plt.plot(days, th3[1:], label='th_l3')

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
pot1 = [sli.matric_pot(t)[0] for t in range(len(sli.get_scaler()))]
pot2 = [sli.matric_pot(t)[1] for t in range(len(sli.get_scaler()))]
pot3 = [sli.matric_pot(t)[2] for t in range(len(sli.get_scaler()))]

plt.plot(days, pot1[1:], label='potential_l1')
plt.plot(days, pot2[1:], label='potential_l2')
plt.plot(days, pot3[1:], label='potential_l3')

#plt.plot(days, qev, label='qev')
plt.title('matric potential')
plt.xlabel('[days]')
#plt.xscale('log')
plt.ylabel('[m]')
plt.grid()
plt.legend()
#plt.ylim(-2.5, 0.5)
plt.show()


################ Temperature temporal ############################

th1 = [sli.T_soil0(t)[0] for t in range(len(sli.get_scaler()))]
th2 = [sli.T_soil0(t)[1] for t in range(len(sli.get_scaler()))]
th3 = [sli.T_soil0(t)[2] for t in range(len(sli.get_scaler()))]

plt.plot(days, th1[1:], label='T_l1')
plt.plot(days, th2[1:], label='T_l2')
plt.plot(days, th3[1:], label='T_l3')

plt.title('Temp')
plt.xlabel('[days]')
plt.ylabel('[C]')
plt.grid()
plt.legend()
plt.show()

###################################################################
cum_qev = np.cumsum(q_ev)

plt.plot(days, q_ev, label='qev')
plt.plot(days, cum_qev, label='cum Evap')

plt.title('Evapotarion')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()

##################################

################### temporal delta ##############################
d1 = [d[0] for d in delta_18O]
d2 = [d[1] for d in delta_18O]
d3 = [d[2] for d in delta_18O]
d4 = [d[3] for d in delta_18O]
d5 = [d[4] for d in delta_18O]

plt.plot(days[1:], np.array(d1), label='delta l1')
plt.plot(days[1:], np.array(d2), label='delta l2')
plt.plot(days[1:], np.array(d3), label='delta l3')
plt.plot(days[1:], np.array(d4), label='delta l4')
plt.plot(days[1:], np.array(d5), label='delta l5')
plt.title('delta')
plt.xlabel('[days]')
plt.ylabel('delta')
plt.grid()
plt.legend()
plt.show()
################### temporal delta ##############################