
import cmf

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cmf_vapor_interface as vp

import matplotlib.cm as cm
import numpy as np

from src import *


def visualize(p, delta, dz, Isotopologue='2H'):

    depth = [-l.center.z for l in p.get_cells()[0].layers]
    plt.figure(figsize=(7, 10))
    for tests in delta.keys():
        plt.plot(delta[tests][Isotopologue][-1][:dz], depth[:dz], label='test_{}'.format(tests))

    plt.xlabel(r'$\delta$ {}'.format(Isotopologue))
    plt.ylabel('depth $[m]$')
    plt.title('Initial testcases after 250 days')
    plt.grid()
    plt.legend(loc='lower right')
    plt.show()


################## cmf ##########################
def cmf_project():
    # Create a project
    p = cmf.project()

    # Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
    c = p.NewCell(0, 0, 0, 1000)
    c.saturated_depth = 0

    return p, c


def rtn_curve():

    # Set up a soil retention curve
    # With parameters as per SISPAT (Braud et al. 2005) Yolo Light clay (Philip, 1957)
    Braud_Ksat = 0.0106272  # 0.0106272  # in m/d
    Braud_phi = 0.35  # porosity
    Braud_alpha = 1/19.3 # 0.005  # 0.002# 0.0040  # 0.005  # 1/19.3  # 0.012 # 0.007 #0.0023  #  air entry potential (Scale value of the water pressure(m))
    Braud_n = 2.22  # 2.22
    Braud_m = 0.099  # 0.1  # 0.4  # 0.25 # 0.27  # 0.45  # 0.099 # 0.20  # 0.055 #0.45  # 0.099
    Braud_theta_r = 0.01
    Braud_eta = 9.14

    vgm = cmf.VGM_BC_RetentionCurve_Windhorst(Ksat=Braud_Ksat, phi=Braud_phi, alpha=Braud_alpha, n=Braud_n, m=Braud_m,
                                              theta_r=Braud_theta_r, eta=Braud_eta)
    vgm.w0 = 0.99999
    vgm.l = 0.67

    # Oversaturation tolerence upto 1% for matrix pot = +1
    # vgm.w0 = 0.99  # See: https://philippkraft.github.io/cmf/cmf_tut_retentioncurve.html Oversaturation
    #vgm.K = vgm.Ksat * pow(vgm.w0, 9.14)

    return vgm


def add_cmf_layers(cell, l_boundaries, r_curve):

    for d in l_boundaries:
        l = cell.add_layer(d, r_curve)
        l.wetness = 1  ##### Check tehts exceeds porosity if  l.wetness = 1.0
        # l.theta = l.porosity
        cell.saturated_depth = 0  # check for wetness = 1, saturated depth = 0


def add_connections(cell, connection=cmf.Richards):
    cell.install_connection(connection)


def cmf_boundary(P):

    cell = P.cells[0]

    summer = cmf.Weather(Tmin=30, Tmax=30, rH=20, wind=2)
    cell.set_weather(summer)

    stress = cmf.ContentStress(theta_d=0.325, theta_w=0.12)
    cell.set_uptakestress(stress)

    # cell.vegetation.Height = 0
    # cell.vegetation.LAI = 0
    # cell.vegetation.CanopyClosure = 0

    # cell.install_connection(cmf.ShuttleworthWallace)
    # cell.install_connection(cmf.PenmanMonteithET)
    # cell.install_connection(cmf.TurcET)
    # cmf.TurcET(cell.layers[0], cell.evaporation)

    ETpot = cmf.timeseries.from_scalar(60)
    cmf.timeseriesETpot(cell.layers[0], cell.evaporation, ETpot)


def cmf_setup():

    # CMF setup
    P, C = cmf_project()  # define project and cell

    L_boundaries = np.array([0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0, 1.25, 1.5,
                           1.9, 2.3, 2.8, 3.5, 4.2, 5.0])

    # L_boundaries = np.cumsum([0.01] * 100)
    add_cmf_layers(C, L_boundaries, rtn_curve())  # define retention curve to all cell layers
    add_connections(C, connection=cmf.Richards)  # apply layer connection

    cmf_boundary(P)

    return P


########################## iso ############################

def _atm(testcase):

    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = 303.17  # sli.Tatm(dt) + Tzero_sli
    Rh_atm = 0.2  # sli.R_humidity_atm(dt)
    wind_speed = 2  # sli.wind_speed(dt)

    #c_iso_2H = iso_delta.delta_to_concentration(0, '2H')  # sli.civa(1) / sli.cva(1)  0.15367056287906838
    #c_iso_18O = iso_delta.delta_to_concentration(0, '18O')

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=0.001,  # canopy height [m] (e.g. 40)
                         d0=0.67 * 0.001,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.1 * 0.001,
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=0,  # leaf area index (e.g. 2.0)
                         extku=0)

    return atm


def _layers(c_iso, C_cmf, testcase):

    # define iso layers as per cmf layer
    layers = C_cmf.layers

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    init_c_iso_2H = iso_delta.delta_to_concentration(ali_2H, '2H')
    init_c_iso_18O = iso_delta.delta_to_concentration(ali_18O, '18O')

    id = 0
    for lr in layers:
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=lr.upper_boundary,
                                                lower_boundary=lr.lower_boundary,
                                                conc_iso_liquid={"2H": init_c_iso_2H, "18O": init_c_iso_18O},
                                                theta=min(lr.theta, lr.porosity),
                                                theta_0=0.01,
                                                theta_sat=lr.porosity,
                                                tortuosity=0.67,
                                                T=303.17,
                                                psi=lr.matrix_potential
                                                )

        id += 1
        c_iso.add_layer(new_layer)

    return c_iso


def heatflux(c_iso, d_t):

    top_node = c_iso.top_layer
    qev = c_iso.q_evap

    q_atm = flux_atmosphere(top_layer=top_node, Ta=c_iso.atmosphere.T)
    q_atm.qevap = c_iso.q_evap

    yb, Tb = q_atm.qyb_qTb()

    q, qh = [-qev], [q_atm.G0()]
    qya, qyb, qTa, qTb = [0], [yb], [0], [Tb]
    qhya, qhyb, qhTa, qhTb = [0], [0], [0], [0]

    for cn in c_iso.connections_l_adv:

        l_node, r_node = cn.left_node, cn.right_node

        hq = heat_flux(left_node=l_node, right_node=r_node, top_layer=top_node)

        dql = hq.qlya() * (l_node.theta - l_node.theta_t0) + hq.qlyb() * (l_node.theta - l_node.theta_t0)
        ql = cn.q_l - dql

        q.append(hq.q_v() + ql), qh.append(hq.qh())

        qya.append(hq.qya()), qyb.append(hq.qyb())
        qTa.append(hq.qTa()), qTb.append(hq.qTb())

        qhya.append(hq.qhya()), qhyb.append(hq.qhyb())
        qhTa.append(hq.qhTa()), qhTb.append(hq.qhTb())

    hs = heat_solve()
    del_cc, del_cch = [], []
    del_ddh, del_ggh = [], []

    for l in c_iso.layers:
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


def h(node):

    try:

        if node.S >= 1:  # saturated

            return node.he

        elif node.S < 1:  # unsaturated

            return node.he * exp(-log(node.S) / node.lamda)
        else:
            raise ValueError

    except ValueError:
        raise NotImplementedError


def _h(S):

    he = -0.193
    # S = (theta - 0.01) / (0.35 - 0.01)

    try:

        if S >= 1:  # saturated

            return he

        elif S < 1:  # unsaturated

            return he * exp(-log(S) / 0.22)
        else:
            raise ValueError

    except ValueError:
        raise NotImplementedError


def _theta(psi):
    alpha = 1/19.3
    n = 2.22
    m = 0.099
    theta_r = 0.00
    theta_sat = 0.35

    S = (1 + (alpha * 100 * psi) ** n) ** -m

    return S * (theta_sat - theta_r) + theta_r


def update_storages(c_iso, c_cmf):

    # atmosphere need not be updated

    # Update layers

    # theta = [min(l.theta, l.porosity) for l in c_cmf.layers]
    T_soil = [303.17] * len(c_cmf.layers)  # np.array([l.T for l in c_iso.layers]) + np.array(dT)
    rH = [None] * len(c_cmf.layers)

    psi = [min(l.matrix_potential, 0) for l in c_cmf.layers]
    theta = [_theta(-p) for p in psi]

    c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c_iso, c_cmf, time):

    """c_iso: iso cell, c_cmf: cmf cell"""
    # Update states and fluxes for
    f = cmf.sec / cmf.day / c_cmf.area  # convert to m per second

    # soil fluxes
    ql = []
    qv = [0] * len(c_iso.layers)
    for lr in c_cmf.layers:
        if lr.upper is not None and lr.lower is None:
            ql.append(0.0 * f)
        else:
            ql.append(lr.flux_to(lr.lower, time) * f)

    # c_iso.update_liquid_fluxes(liquid_fluxes=ql)
    c_iso.update_vapor_fluxes(vapor_fluxes=None)  #vapor_fluxes=None: self compute vapor flux
    qv = np.array(c_iso.vapor_fluxes + [0])
    mask = np.abs(qv) > np.abs(ql)
    ql_liquid = np.where(mask, 0, ql - qv)
    qv_update = np.where(mask, ql, qv)
    c_iso.update_vapor_fluxes(vapor_fluxes=qv_update)

    c_iso.update_liquid_fluxes(liquid_fluxes=ql_liquid)

    f = 1 / c_cmf.area / cmf.day.AsSeconds()  # m3 day-1 to ms-1
    q_evap = c_cmf.evaporation(time) * f

    atm = flux_atmosphere(atmosphere=c_iso.atmosphere, top_layer=c_iso.layers[0])
    f_ql, f_qv = atm.E_liquid() / atm.E_total(), atm.E_vapor() / atm.E_total()

    # Surface variables
    T_surface = 303.17
    ql_surface = - q_evap * f_ql
    qv_surface = - q_evap * f_qv
    c_iso.update_evaporation(q_ev=q_evap, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)


def update_theta_cmf(c_iso, c_cmf, time_step):

    qv = np.array(c_iso.vapor_fluxes)
    dqv = -np.diff(qv)
    delta_qv = np.concatenate(([-qv[0]], dqv, [qv[-1]]))  # net flux into a layer, first layer -ve due to  downward positive flux

    delta_theta = delta_qv * time_step.seconds / c_cmf.layers.thickness

    updated_theta = c_cmf.layers.theta + delta_theta

    for l, theta in zip(c_cmf.layers, updated_theta):

        #t = l.theta + d_theta
        l.wetness = theta / l.porosity

    return


def iso_setup(testcase):

    P = cmf_setup()  # cmf project
    C = P.cells[0]  # current cell cmf

    # setup for isotope (isotope storages)
    atm = _atm(testcase)
    # Iso project
    p = iso_project()
    p.new_cell(atmosphere=atm, area=C.area, x=C.x, y=C.y, z=C.z)

    c = p.get_cells()[0]  # get current iso cell
    _layers(c, C, testcase)  # add cmf layers to the current iso_cell

    #####Install connections#######
    c.install_connections(vapor_diffusion=True, vapor_advection=True)  # install storage connections between the layers
    c.add_evaporation()


    return p, P


def run(p, P, sim_period=50, dt=1, solutes=["2H", "18O"], testcase=1):

    """"
    sim_period: days
    dt: hours
    """

    print('Testcase:', testcase)
    kwargs = iso_delta.test_case_args(testcase=testcase)

    C = P.cells[0]  # current cell cmf_project
    c = p.get_cells()[0]  # current cell of iso_project

    f = cmf.sec / cmf.day / C.area  # convert to m per second  #  m / s to mm / day (m3 / day)

    ### vapor boundary ###
    v = vp.Vaporizer(cell=C)

    # Define solver
    solver = cmf.CVodeBanded(P)
    start = datetime(2024, 1, 1)
    end = start + timedelta(days=sim_period)
    timestep = timedelta(hours=dt)

    delta_t = timestep.seconds  # dt * 3600  # dt.AsSeconds()

    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}
    ev, theta, ql, qv = [], [], [], []
    cv_sat, dv, rh, T, pot, m_pot = [], [], [], [], [], []
    ql_surface, qv_surface = [], []
    for t in solver.run(start, end, timestep):

        print(t)
        print(C.layers.theta.tolist())

        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        update_storages(c_iso=c, c_cmf=C), update_boundaries(c_iso=c, c_cmf=C, time=t)
        for solute in solutes:

            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=None, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

        ev.append(c.q_evap)
        theta.append([l.theta for l in c.layers])
        ql.append(c.liquid_fluxes)
        qv.append(c.vapor_fluxes)
        cv_sat.append([l.cv_sat for l in c.layers])

        fn = flux_node(top_node=c.layers[0])
        dv.append([fn.dv(l) for l in c.layers])
        rh.append([l.relative_humidity(l.psi, l.T) for l in c.layers])
        T.append([l.T for l in c.layers])
        #pot.append([l.psi for l in c.layers])
        m_pot.append([l.matrix_potential for l in C.layers])
        pot.append([l.psi for l in c.layers])
        ql_surface.append(c.ql_surface), qv_surface.append(c.qv_surface)

        # update cmf theta
        # update_theta_cmf(c_iso=c, c_cmf=C, time_step=timestep)
        # delta_theta.append(dth)
        # th_updated.append([l.theta for l in C.layers])

        ##### vapor boundary #####
        # q_v = np.array(c.vapor_fluxes)
        # dqv = - np.diff(q_v)
        # delta_qv = np.concatenate(([-q_v[0]], dqv, [q_v[-1]]))  # net flux into a layer, first layer -ve due to  downward positive flux
        # v.flux = - delta_qv / f  # np.array(q_v) / f

    return c_iso_delta, [ev, theta, ql, qv], [cv_sat, dv, rh, T, pot, m_pot, ql_surface, qv_surface]


def run_testcases(test_cases):

    """"
    simulatio period: days
    dt: hours
    """

    delta = {}
    for Testcase in test_cases:

        delta[Testcase] = {}

        p_iso, p_cmf = iso_setup(testcase=Testcase)
        d, X, Y = run(p_iso, p_cmf, sim_period=250, dt=12, testcase=Testcase)

        # delta at the end of simulation for each test cases
        delta[Testcase]["2H"] = d["2H"]
        delta[Testcase]["18O"] = d["18O"]

    return p_iso, delta, X, Y


p_iso, delta, X, Y = run_testcases(test_cases=[1, 2, 3, 4, 5, 6])

delta_2H_1, delta_18O_1 = delta[1]['2H'], delta[1]['18O']
delta_2H_6, delta_18O_6 = delta[1]['2H'], delta[1]['18O']

dz = 10

# Iso_delta
visualize(p_iso, delta, dz, Isotopologue='2H')
visualize(p_iso, delta, dz, Isotopologue='18O')

f = 86400 * 1000  # factor ms-1 to mm day-1
ev, theta, ql, qv = X
cv, dv, rh, T, pot, m_pot, ql0, qv0 = Y

depth = [-l.center.z for l in p_iso.get_cells()[0].layers][:dz]
days = [50, 100, 150, 200, 250]
t_step = [(d * 2)-1 for d in days]
time = np.arange(0.5, 250.5, 0.5)


###################### Iso ######################
############ delta profile ##################
_delta_18O = [delta_18O_1[t] for t in t_step]
for d, l in zip(_delta_18O, days):
    plt.plot(d[:dz], depth, label=str(l) + ' days')

plt.title('delta')
#plt.xscale('log')
plt.xlabel('delta')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


#######################################################################
############# theta profile ########

_theta = [theta[t] for t in t_step]
#_theta_up = [th_updated[t][:10] for t in t_step]
for th, l, in zip(_theta, days):
    plt.plot(th[:dz], depth, label=str(l) + ' days')
    #plt.plot(thup, depth, color=c, linestyle='--')

plt.title('theta')
plt.xlabel('[m3 / m3 ]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


_ql = [ql[t] for t in t_step]
_qv = [qv[t] for t in t_step]


""""
############ ql profile ##########
for q_l, l in zip(_ql, days):
    plt.plot(np.array(q_l[:dz]) * f, depth, label=str(l) + ' days')

plt.title('ql')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()

############ qv profile ##################
for q_v, l in zip(_qv, days):
    plt.plot(np.array(q_v[:dz]) * f, depth, label=str(l) + ' days')

plt.title('qv')
#plt.xscale('log')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()

"""
############ coupled profile #####################
plt.plot(np.array(_ql[-1][:dz]) * f, depth, label='ql')
plt.plot(np.array(_qv[-1][:dz]) * f, depth, label='qv')
plt.title('ql_qv after 250 days')
#plt.xscale('log')
plt.xlabel('[mm per day]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()

""""
############### potential profile ######################
_potential = [pot[t] for t in t_step]

for pt, l in zip(_potential, days):
    plt.plot(np.array(pt[:dz]), depth, label=str(l) + ' days')


plt.title('matrix potential')
#plt.xscale('log')
plt.xlabel('[m]')
plt.ylabel('depth [m]')
plt.legend()
plt.grid()
plt.show()


##########################################################
############ temporal Evaporation ########################
plt.plot(time, np.array(ev) * f, label='evaporatiom')
plt.title('qevap')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()
"""
############ surface fluxes ########################
plt.plot(time, np.array(ev) * f, label='evaporatiom')
plt.plot(time, np.array(ql0) * f, label='ql0')
plt.plot(time, np.array(qv0) * f, label='qv0')
plt.title('qevap')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()


"""
########## temporatl theta ###################
th1 = [t[0] for t in theta]
th2 = [t[1] for t in theta]
th3 = [t[2] for t in theta]
th4 = [t[3] for t in theta]
th5 = [t[4] for t in theta]

# plt.plot(time, th1, label='theta l1')
plt.plot(time, th2, label='theta l2')
plt.plot(time, th3, label='theta l3')
plt.plot(time, th4, label='theta l4')
plt.plot(time, th5, label='theta l5')

plt.title('theta')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()

############# temporal ql qv ################## 
ql1 = [q[0] * f for q in ql]
ql2 = [q[1] * f for q in ql]
ql3 = [q[2] * f for q in ql]
ql4 = [q[3] * f for q in ql]

qv1 = [q[0] * f for q in qv]
qv2 = [q[1] * f for q in qv]
qv3 = [q[2] * f for q in qv]
qv4 = [q[3] * f for q in qv]

# plt.plot(time[:40], ql1[:40], label='ql_l1')
# plt.plot(time[:40], ql2[:40],  label='ql_l2')
# plt.plot(time[:40], ql3[:40], label='ql_l3')
# plt.plot(time[:40], ql4[:40], label='ql_l4')
plt.plot(time, qv1,  linestyle='--', label='qv_l1')
plt.plot(time, qv2,  linestyle='--', label='qv_l2')
plt.plot(time, qv3,  linestyle='--', label='qv_l3')
plt.plot(time, qv4,  linestyle='--', label='qv_l4')

plt.title('ql_qv')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.grid()
plt.legend()
plt.show()

########## temporal qv components ##################

_cv1 = np.array([c[0] for c in cv])
_cv2 = np.array([c[1] for c in cv])
_cv3 = np.array([c[2] for c in cv])
_cv4 = np.array([c[3] for c in cv])

_dv1 = np.array([c[0] for c in dv])
_dv2 = np.array([c[1] for c in dv])
_dv3 = np.array([c[2] for c in dv])
_dv4 = np.array([c[3] for c in dv])

# plt.plot(time, _dv1, label='dv_l1')
plt.plot(time[:40], _dv2[:40], label='dv_l2')
plt.plot(time[:40], _dv3[:40], label='dv_l3')
plt.plot(time[:40], _dv4[:40], label='dv_l4')
plt.legend()
plt.grid()
plt.show()

_rh1 = np.array([c[0] for c in rh])
_rh2 = np.array([c[1] for c in rh])
_rh3 = np.array([c[2] for c in rh])
_rh4 = np.array([c[3] for c in rh])

#plt.plot(time[:40], _rh1[:40], label='rh_l1')
plt.plot(time[:40], _rh2[:40], label='rh_l2')
plt.plot(time[:40], _rh3[:40], label='rh_l3')
plt.plot(time[:40], _rh4[:40], label='rh_l4')
plt.legend()
plt.grid()
plt.show()

_T1 = np.array([c[0] for c in T])
_T2 = np.array([c[1] for c in T])

_pot1 = np.array([c[0] for c in pot])
_pot2 = np.array([c[1] for c in pot])
_pot3 = np.array([c[2] for c in pot])
_pot4 = np.array([c[3] for c in pot])
_pot5 = np.array([c[4] for c in pot])

plt.plot(time, _pot1, label='pot_l1')
plt.plot(time, _pot2, label='pot_l2')
plt.plot(time, _pot3, label='pot_l3')
plt.plot(time, _pot4, label='pot_l4')
plt.plot(time, _pot5, label='pot_l5')

plt.legend()
plt.grid()
plt.show()

_m_pot1 = np.array([c[0] for c in m_pot])
_m_pot2 = np.array([c[1] for c in m_pot])

plt.plot(time, _T1, label='T_l1')
plt.plot(time, _T2, label='T_l2')
plt.legend()
plt.grid()
plt.show()


############### qv - theta #################

th1 = [t[0] for t in theta]
th2 = [t[1] for t in theta]
th3 = [t[2] for t in theta]
th4 = [t[3] for t in theta]

qv1 = [t[0] * f for t in qv]
qv2 = [t[1] * f for t in qv]
qv3 = [t[2] * f for t in qv]
qv4 = [t[3] * f for t in qv]


ql1 = [t[0] * f for t in ql]
ql2 = [t[1] * f for t in ql]

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

# l1, = ax1.plot(time[:20], th1[:20], color='r', label='theta l1')
l2, = ax1.plot(time[:20], th2[:20], label='theta l2')
l3, = ax1.plot(time[:20], th3[:20],  label='theta l3')
l4, = ax1.plot(time[:20], th4[:20], label='theta l4')

l5, = ax2.plot(time[:20], qv2[:20], label='qv l2')
l6, = ax2.plot(time[:20], qv3[:20], label='qv l3')
l7, = ax2.plot(time[:20], qv4[:20], label='qv l4')


labels = [l2, l3, l4, l5, l6, l7]
ax1.legend(labels, [l.get_label() for l in labels])
ax1.set_xlabel('[days]')
ax1.set_ylabel('[m3 / m3]')
ax2.set_ylabel('[mm per day]')

plt.title('theta')
plt.grid()
plt.show()


################temporal delta ########################
d1 = [d[0] for d in delta_2H_1]
d2 = [d[1] for d in delta_2H_1]
d3 = [d[2] for d in delta_2H_1]
d4 = [d[3] for d in delta_2H_1]
d5 = [d[4] for d in delta_2H_1]

plt.plot(time, np.array(d1), label='delta l1')
plt.plot(time, np.array(d2), label='delta l2')
plt.plot(time, np.array(d3), label='delta l3')
plt.plot(time, np.array(d4), label='delta l4')
plt.plot(time, np.array(d5), label='delta l5')
plt.title('delta')
plt.xlabel('[days]')
plt.ylabel('delta')
plt.grid()
plt.legend()
plt.show()


############ cumulative evaporation ########################
plt.plot(time, np.array(ev) * f, label='evaporatiom')
plt.plot(time, np.cumsum(ev) * f, label='cum evaporatiom')
plt.title('cum Evaporation')
plt.xlabel('[days]')
plt.ylabel('[mm per day]')
plt.yscale('log')
plt.grid()
plt.legend()
plt.show()


"""



