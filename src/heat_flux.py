
from math import exp
import numpy as np


def vapor_flux_to_mass_flux(q_vap, T, RH):
    """
    Convert vapor flux from [m/s] to [kg/m2/s].

    Parameters
    ----------
    q_vap : float
        Vapor flux [m/s]
    T_C : float
        Temperature in K
    RH : float
        Relative humidity [0-1]

    Returns
    -------
    J_vap : float
        Mass vapor flux [kg/m2/s]
    rho_v : float
        Vapor density [kg/m3]
    """
    # Constants
    Mw = 0.018015  # kg/mol
    R = 8.314  # J/mol/K

    # Convert temperature to Celcius
    T_C = T - 273.15

    # Saturation vapor pressure (Pa) – Tetens formula
    p_sat = 610.78 * exp((17.27 * T_C) / (T_C + 237.3))

    # Partial vapor pressure (Pa)
    p_v = RH * p_sat

    # Vapor density [kg/m3]
    rho_v = (p_v * Mw) / (R * T)

    # Convert flux [m/s] → [kg/m2/s]
    J_vap = rho_v * q_vap

    return J_vap


def compute_surface_temperature(T1, dz1, lam1, evap=0.0, T_atm=25.0, h_c=10.0, L_v=2.45e6):
    """
    Compute soil surface temperature from a simple surface energy balance.

    Parameters
    ----------
    T1 : float
        Temperature of first soil node [°C]
    dz1 : float
        Thickness of first soil layer [m]
    lam1 : float
        Thermal conductivity of first layer [W/m/K]
    evap : float
        Evaporation flux [kg/m2/s]
    T_atm : float
        Atmospheric temperature [°C]
    h_c : float
        Convective heat transfer coefficient [W/m2/K]
    L_v : float
        Latent heat of vaporization [J/kg]

    Returns
    -------
    T_s : float
        Computed surface temperature [°C]
    """
    # Heat conduction into soil (top half-layer)
    # G = lam * (T_s - T1) / (dz1/2)
    # Sensible heat flux to atmosphere: H = h_c * (T_s - T_atm)
    # Latent heat flux: LE = L_v * evap

    # Energy balance: H + LE + G = 0
    # Solve for T_s:
    T_s = (lam1 * T1 / (dz1/2) + h_c * T_atm - L_v * evap) / (lam1 / (dz1/2) + h_c)

    return T_s


def compute_temperature(T, theta, RH, q_liq, q_vap, soil, dt, dz,
                               T_atm=25.0, q_evap=0.0, h_c=10.0, S=None):
    """
    Update soil temperature in a 1-D layered column with calculated surface temperature.
    Safe indexing is used for advection terms to avoid out-of-bounds errors.
    """

    N = len(T)
    rho_s, c_s = soil['rho_s'], soil['c_s']
    rho_w, c_w = soil['rho_w'], soil['c_w']
    lambda_dry, lambda_sat = soil['lambda_dry'], soil['lambda_sat']

    #
    # vapor flux [m/s] to [kg/m2/s]
    J_vap = np.array([vapor_flux_to_mass_flux(qv, t, rh)
                      for qv, t, rh in zip(q_vap, T, RH)])
    J_vap = np.insert(J_vap, 0, q_vap[0])

    # evap flux from [m/s] to [kg/m2/s]
    evap = np.array(-q_evap) * 1000

    # Volumetric heat capacity and thermal conductivity
    C = (1 - theta) * rho_s * c_s + theta * rho_w * c_w
    lam = lambda_dry + (lambda_sat - lambda_dry) * (theta**1.5)
    L_v = 2.45e6  # latent heat

    # Compute surface temperature
    T_surf = compute_surface_temperature(T[0], dz[0], lam[0], evap, T_atm, h_c, L_v)

    # Build tridiagonal system
    A = np.zeros((N, N))
    RHS = np.zeros(N)

    for i in range(N):

        dz_i = dz[i]

        # Conduction coefficients
        if i > 0:
            lam_im = 0.5 * (lam[i] + lam[i-1])
            dz_im = 0.5 * (dz[i] + dz[i-1])
            A[i, i-1] = -lam_im / (dz_im * dz_i)
        if i < N-1:
            lam_ip = 0.5 * (lam[i] + lam[i+1])
            dz_ip = 0.5 * (dz[i] + dz[i+1])
            A[i, i+1] = -lam_ip / (dz_ip * dz_i)
        A[i, i] = C[i]/dt - (A[i, i-1] if i > 0 else 0) - (A[i, i+1] if i < N-1 else 0)

        # RHS: old temperature + advection + latent + sources
        RHS[i] = C[i]*T[i]/dt

        # --- Advection terms with safe upwind indexing ---
        # Liquid
        T_up_in  = T[i-1] if i > 0 else T[i]       # bottom interface
        T_up_out = T[i+1] if i < N-1 else T[i]    # top interface
        adv_liq = rho_w * c_w * (q_liq[i+1]*T_up_out - q_liq[i]*T_up_in)

        # Vapor
        adv_vap = 1850 * (J_vap[i+1]*T_up_out - J_vap[i]*T_up_in)

        # Latent heat from vapor flux
        latent = L_v * (J_vap[i+1] - J_vap[i])

        RHS[i] += -(adv_liq + adv_vap)/dz_i + latent/dz_i

        if S is not None:
            RHS[i] += S[i]

    # Top boundary: surface temperature as Dirichlet BC
    A[0, :] = 0.0
    A[0, 0] = 1.0
    RHS[0] = T_surf

    # Bottom boundary: fixed temperature (can be replaced by flux BC later)
    A[-1, :] = 0.0
    A[-1, -1] = 1.0
    RHS[-1] = T[-1]

    # Solve linear system
    T_new = np.linalg.solve(A, RHS)
    return T_new, T_surf

