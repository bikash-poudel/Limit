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


def _layers(c, testcase):

    lower_boundaries = np.array([0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0, 1.25, 1.5,
                            1.9, 2.3, 2.8, 3.5, 4.2, 5.0])

    theta = [0.34999] + [0.35] * (len(lower_boundaries) - 1)
    theta_r = 0.01
    theta_sat = 0.35
    tortuosity = 0.67
    T_soil = 303.17

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    init_c_iso_2H = iso_delta.delta_to_concentration(ali_2H, '2H')
    init_c_iso_18O = iso_delta.delta_to_concentration(ali_18O, '18O')
    id = 0
    upper_boundary = 0

    for lower_boundary, th in zip(lower_boundaries, theta):
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=upper_boundary,
                                                lower_boundary=lower_boundary,
                                                conc_iso_liquid={"2H": init_c_iso_2H, "18O": init_c_iso_18O},
                                                theta=th,
                                                theta_0=theta_r,
                                                theta_sat=theta_sat,
                                                tortuosity=tortuosity,
                                                T=T_soil
                                                )
        id += 1
        upper_boundary = lower_boundary
        c.add_layer(new_layer)

    return c


def _atm(testcase):

    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = 303.17
    Rh_atm = 0.2
    wind_speed = 2

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=10,  # 10,  # canopy height [m] (e.g. 40)
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


def h(S):

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


def q_surface(atm):

        qev = atm.q_evaporation()

        if atm.E_total() >= atm.ev_pot():
            qls = - qev
            qvs = 0
        elif atm.E_total() < atm.ev_pot():
            qls = - atm.ql_evap()
            qvs = - atm.qv_evap()
        else:
            raise NotImplementedError

        return qev, qls, qvs


def solve_coupled(c, dt):

    atm = flux_atmosphere(c.atmosphere, c.top_layer)
    qev, qls, qvs = q_surface(atm)

    yb, Tb = atm.qyb_qTb()

    q, qv, qh = [-qev], [qvs], [atm.G0()]
    qya, qyb, qTa, qTb = [0], [yb], [0], [Tb]
    qvya, qvyb, qlyb = [0], [atm.qvyb()], [atm.qlyb()]
    qhya, qhyb, qhTa, qhTb = [0], [0], [0], [0]

    for cn in c.connections_v_adv:

        l_node, r_node = cn.left_node, cn.right_node
        # hv = heat_flux(left_node=l_node, right_node=r_node, top_layer=c.top_layer)
        vp = vapor_flux(left_node=l_node, right_node=r_node, top_layer=c.top_layer)

        q.append(vp.q()), qv.append(vp.q_v()), # qh.append(vp.qh())

        qya.append(vp.qya()), qyb.append(vp.qyb())
        qTa.append(vp.qTa()), qTb.append(vp.qTb())

        # qhya.append(vp.qhya()), qhyb.append(vp.qhyb())
        # qhTa.append(vp.qhTa()), qhTb.append(vp.qhTb())

        qvya.append(vp.qvya()), qvyb.append(vp.qvyb())


    # boundary zero fux
    q.append(0), qv.append(0), qh.append(0)
    qya.append(0), qyb.append(0), qTa.append(0), qTb.append(0)

    qhya.append(0), qhyb.append(0), qhTa.append(0), qhTb.append(0)
    qvya.append(0), qvyb.append(0)

    vs = vapor_solve()
    del_cc, del_cch = [], []
    del_ddh, del_ggh = [], []

    for l in c.layers:
        del_cc.append(vs.delta_cc(node=l, dt=dt))
        #del_cch.append(hs.delta_cch(node=l, dt=dt))
        # del_ddh.append(hs.delta_ddh(node=l, dt=dt))
        # del_ggh.append(hs.delta_ggh(node=l, dt=dt))

    aa, bb, ee, ff = np.array(qya[1:-1]), np.array(qTa[1:-1]), - np.array(qyb[1:-1]), - np.array(qTb[1:-1])

    cc = np.array(qyb[:-1]) - np.array(qya[1:]) - np.array(del_cc)
    dd = np.array(qTb[:-1]) - np.array(qTa[1:])
    gg = - (np.array(q[:-1]) - np.array(q[1:]))

    X = [aa, bb, cc, dd, ee, ff, gg]
    dy = vs.solve_sparse(X)

    return dy, [qev, qls, qvs], [q, qv], [qya, qyb, qTa, qTb, qvya, qvyb, qlyb]


def update_storages(c, dy):

    theta = update_S(c, dy)

    psi = [h(l.S) for l in c.layers]
    T = [303.17] * len(c.layers)
    rH = None

    c.update_layers(theta=theta, T=T, rH=rH, psi=psi)


def update_boundaries(c, q_surface, q_layers, d_q, dy):

    qev, qls, qvs = q_surface
    q, qv = q_layers
    qya, qyb, qTa, qTb, qvya, qvyb, qlyb = d_q

    # Surface variables
    T_surface = 303.17
    q_ev = qev - qyb[0] * dy[0]
    ql_surface = qls + qlyb[0] * dy[0]
    qv_surface = qvs + qvyb[0] * dy[0]

    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

    # soil fluxes
    q_tot = np.array(q[1:-1]) + np.array(qya[1:-1]) * dy[:-1] + np.array(qyb[1:-1]) * dy[1:]
    qv = np.array(qv[1:-1]) + np.array(qvya[1:-1]) * dy[:-1] + np.array(qvyb[1:-1]) * dy[1:]

    qn = np.array(q[-1]) + np.array(qya[-1]) * dy[-1]
    q_tot, qv = np.append(q_tot, qn), np.append(qv, 0)

    ql = q_tot - qv

    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)


def run_iso(p, sim_period=250, **kwargs):

    dt = 60
    T_last = sim_period * 86400 + dt

    c = p.get_cells()[0]
    for dt in range(dt, T_last, dt):

        print(dt / 86400, ' days')
        dy, qsurface, q_layers, dq = solve_coupled(c, dt)

        update_storages(c, dy)
        update_boundaries(c, qsurface, q_layers, dq, dy)

        print([l.theta for l in c.layers])


def iso_setup(testcase=1, **kwargs):

    # Define boundary storages
    atm = _atm(testcase)  # atmosphere

    p = iso_project()  # create a project
    p.new_cell(atmosphere=atm, area=1, x=0, y=0, z=0)  # add new cell

    c = p.get_cells()[0]  # get current cell
    _layers(c, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation()

    return run_iso(p, **kwargs)


run = iso_setup()


