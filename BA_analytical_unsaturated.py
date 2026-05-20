

import numpy as np
from math import exp, log

import iso_cmf as iso
import iso_cmf_Barnes_Alison as BA
from src import iso_fluxes as isof

from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt


def get_cv_sat(T):

    """
    Returns the saturated water vapor volumetric mass (m**3_H20/m**3) .

    T_atmosphere = Temperature in Kelvin

    pv_sat = saturated water vapor volumetric mass (m**3_H20/m**3)

    TODO: Check for which temperature range it is valid
    TODO: recheck the result with SLI
    Control status: Checked on 03.06.2024 --> Results are the same as for SLI
    """

    # SLI: sli_utils.f90::L1416
    Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
    R = 8.3142995834350586  # universal gas constant (j / mol / k)

    cv_sat = p_sat(T) * Mw / R / T / 1000  # m3/m3
    # cv_sat = 0.002166 * self.p_sat / self.T * 1000  # m3/m3

    return cv_sat


def p_sat(T):

    """
    Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch and valid over the range 0 to 60 Degree celsius

    T_atmosphere = Temperature in Kelvin

    p_sat = Saturation vapour pressure, ps, in pascal

    <--checked and own implementation is okay for temp. above 273.15 Kelvin
    """
    # SLI: cable_sli_utils.f90::L2195-2195
    Tzero = 273.16000366210938
    p_sat = 610.59999465942383 * exp(
        17.270000457763672 * (T - Tzero) /
        ((T - Tzero) + 237.30000305175781))  # Saturation vapour pressure, ps, in pascal

    if T < 273.16:  # below 0 Degree celsius
        p_sat = -4.86 + 0.855 * p_sat + 0.000244 * p_sat ** 2  # ps for ice

    return p_sat


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


def semi_analytical(z, T_z, h_z, D_l, D_v, q_evap, alpha_i, alpha_ik, delta_surface, delta_alim):

    # Constants
    rho_w = 1000  # Water density [kg/m³]
    g = 9.81  # Gravity [m/s²]
    R = 461.5  # Gas constant [J/kg/K]

    E = q_evap  # Evaporation rate [m/s]

    h_u = np.exp(g * h_z / (R * T_z))

    rho_vsat = np.exp(31.3716 - 6014.79 / T_z - 0.00792495 * T_z) / T_z * 1e-3

    z_1 = D_l / E
    z_v = D_v * rho_vsat / (rho_w * E)

    denom = z_1 + h_u * z_v

    P_z = 1.0 / denom

    log_term = np.log(h_u * rho_vsat * (alpha_ik - alpha_i))

    d_log_term_dz = np.gradient(log_term, z)

    Q_z = (h_u * z_v * (alpha_ik - alpha_i) / denom) * d_log_term_dz + P_z * delta_alim

    # Initialize solution array

    delta_analytical = np.zeros_like(z)
    delta_analytical[0] = delta_surface  # Your input surface value

    # Forward Euler integration
    for i in range(len(z) - 1):
        dz = z[i + 1] - z[i]
        delta_analytical[i + 1] = delta_analytical[i] + dz * (-P_z[i] * delta_analytical[i] + Q_z[i])

    return z, delta_analytical


def solve_analytical(cell, Isotopologues=["2H", "18O"], delta_ali={"2H": 0, "18O": 0}, **kwargs):

    v_diff = isof.vapor_diffusion_base_class()
    l_diff = isof.liquid_diffusion_base_class()

    ev = cell.connection_evap
    atm = cell.atmosphere

    layers = cell.layers
    z = np.array([l.upper_boundary for l in layers])
    Temp = np.array([l.T for l in layers])
    head = np.array([l.psi for l in layers])

    # dz = np.array([l.thickness for l in layers])
    dz = np.array([l.thickness for l in layers])

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

        div = np.array([v_diff.dv_soil_air(T=l.T, theta=l.theta, theta_sat=l.theta_sat, tortuosity=l.tortuosity)
                        for l in cell.layers])
        dil = np.array([l_diff.dl_i(T=l.T, Isotopologue=solute, theta=l.theta, tortuosity=l.tortuosity)
                        for l in cell.layers])
        alpha_i = np.array([l.alpha_i(Isotopologue=solute, **kwargs) for l in layers])

        dz, dlt = semi_analytical(z, T_z=Temp, h_z=head,
                                  D_l=dil, D_v=div, q_evap=cell.q_evap,
                                  alpha_i=alpha_i,
                                  alpha_ik=alpha_ik,
                                  delta_surface=delta_atm[solute],
                                  delta_alim=delta_ali[solute])

        delta[solute].append(dlt)

    return dz, d_s, delta


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

        BA.update_boundaries(c_iso=c, c_cmf=C, time=t, dt=dt)
        BA.update_storages(c_iso=c, c_cmf=C, dt=dt)

    #############################################################################################
    ############################# run solver ####################################################

    dz, ds, dlt = solve_analytical(cell=c, Isotopologues=solutes, delta_ali=delta_ali, **kwargs)

    delta["2H"].append(dlt['2H'][0])
    delta["18O"].append(dlt['18O'][0])

    ############################################################################################
    ############################################################################################

    return c, dz, ds, delta


# cell, dz, ds, d = setup_run(dt=1, sim_period=250)  # hour, days


# d_2H, d_18O = d["2H"][0], d["18O"][0]

"""
plt.plot(d_18O, -dz, label='18O')
plt.xlabel('delta')
plt.ylabel('depth m')
plt.legend()
plt.show()

plt.plot(d_2H, -dz, label='2H')
plt.xlabel('delta')
plt.ylabel('depth m')
plt.legend()
plt.show()

"""

