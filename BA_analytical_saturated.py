

import numpy as np
import iso_cmf as iso
import iso_cmf_Barnes_Alison as BA
from src import iso_fluxes as isof
import matplotlib.pyplot as plt


def delta_surface(alpha_i, alpha_iK, h_a, delta_alim, delta_v):
    """
    Compute surface isotopic composition δi_s (Eq. 41, BA83).

    Parameters
    ----------
    alpha_i : float
        Equilibrium fractionation factor (dimensionless)
    alpha_iK : float
        Kinetic fractionation factor (dimensionless)
    h_a : float
        Relative humidity (0–1)
    delta_alim : float
        Isotopic ratio of alimentation water [permil or chosen units]
    delta_v : float
        Isotopic ratio of atmospheric vapor [permil or chosen units]
    """
    num = (1 - h_a) * alpha_iK * (1 + delta_alim/1000) + h_a * (1 + delta_v / 1000)
    return 1000 * ((num / alpha_i) - 1)


def isotope_profile(z, delta_alim, delta_surf, q_ev, Di_l, rho_w=1000.0):
    """
    Analytical isotopic profile for saturated soil (Braud et al. 2005, Eq. 40).

    Parameters
    ----------
    z : array_like
        Depth in soil [m]
    delta_alim : float
        Isotopic ratio of alimentation water [permil or chosen units]
    delta_surf : float
        Isotopic ratio at soil surface [permil or chosen units]
    E_ms : float
        Evaporation flux [m/s]
    Di_l : array
        Diffusion coefficient of isotope in water [m^2/s]
    rho_w : float
        Water density [kg/m^3], default=1000
    """
    # Convert E from m/s to kg m-2 s-1
    E = q_ev
    return delta_alim + (delta_surf - delta_alim) * np.exp(- (E / Di_l) * z)


def solve_analytical(cell, Isotopologues=["2H", "18O"], delta_ali={"2H": 0, "18O": 0}, **kwargs):

    v_diff = isof.vapor_diffusion_base_class()
    l_diff = isof.liquid_diffusion_base_class()

    ev = cell.connection_evap
    atm = cell.atmosphere
    l = cell.layers[0]  # top layer

    dz = np.array([l.upper_boundary for l in cell.layers])

    ########### surface concentration ###########################
    ha = cell.atmosphere.Rh
    qi = cell.q_evap

    print(qi)

    delta_atm = {"2H": -100, "18O": -14}
    #ds = {"2H": -67, "18O": -31.9}
    d_s, delta = {}, {}

    for solute in Isotopologues:

        d_s[solute], delta[solute] = [], []

        d_v = v_diff.dv_free_air(T=atm.T)
        dv_i = v_diff.dv_i(dv=d_v, Isotopologue=solute, **kwargs)

        alpha_ik = ev.alpha_i_k(dv=d_v, dv_i=dv_i, formulation="BarnesAllsion", **kwargs)

        ds = delta_surface(alpha_i=l.alpha_i(Isotopologue=solute, **kwargs),
                           alpha_iK=alpha_ik,
                           h_a=ha,
                           delta_alim=delta_ali[solute],
                           delta_v=delta_atm[solute])

        dil = np.array([l_diff.dl_i(T=l.T, Isotopologue=solute, theta=l.theta, tortuosity=l.tortuosity)
                        for l in cell.layers])

        dlt = isotope_profile(z=dz,
                              delta_alim=delta_ali[solute],
                              delta_surf=ds,
                              q_ev=qi,
                              Di_l=dil
                              )
        conc = np.array([iso.iso_delta.delta_to_concentration(d, solute)
                         for d in dlt])

        cell.update_c_layers(conc_iso=conc, Isotopologue=solute)

        d_s[solute].append(ds)
        delta[solute].append(dlt)

    return delta, dz


def setup_run(dt=1, sim_period=250, solutes=["2H", "18O"], delta_ali={"2H": 0, "18O": 0}):

    p, P = BA.iso_setup()

    """
    P: cmf project
    p: iso project
    simulation period: days
    dt: hours
    """

    kwargs = iso.iso_delta.test_case_args(testcase=6, BA=True)

    C = P.cells[0]  # current cell cmf_project
    c = p.get_cells()[0]  # current cell of iso_project

    ###################### Solver #######################################################
    solver = iso.cmf.CVodeBanded(P)
    start = iso.datetime(2024, 1, 1)
    end = start + iso.timedelta(days=sim_period)
    timestep = iso.timedelta(hours=dt)

    ################# output variables  #################################################
    c_iso_delta = {'2H': [], '18O': []}
    ######################################################################################

    delta = {"2H": [], "18O": []}
    for t in solver.run(start, end, timestep):
        print(t)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        iso.update_storages(c_iso=c, c_cmf=C)
        BA.update_boundaries(c_iso=c, c_cmf=C, time=t, dt=dt)

        #############################################################################################
        ############################# run solver ####################################################

        dlt, dz = solve_analytical(cell=c, Isotopologues=solutes, delta_ali=delta_ali, **kwargs)

        delta["2H"].append(dlt['2H'][0])
        delta["18O"].append(dlt['18O'][0])

        ############################################################################################
        ############################################################################################

    return delta, dz


d, dz = setup_run(dt=1, sim_period=250)  # hour, days

d_2H, d_18O = d["2H"][-1], d["18O"][-1]

plt.plot(d_2H, -dz, label='2H')
plt.xlabel('delta')
plt.ylabel('depth m')
plt.show()

plt.plot(d_18O, -dz, label='2H')
plt.xlabel('delta')
plt.ylabel('depth m')
plt.show()

