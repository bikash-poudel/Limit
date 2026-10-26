
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


def upstream_temperature(face_flux, T, T_top, T_bottom):

    face_flux = np.asarray(face_flux, dtype=float)
    T = np.asarray(T, dtype=float)

    if face_flux.shape != (len(T) + 1,):
        raise ValueError("Expected N+1 interface fluxes for N layers.")

    donor = np.empty(len(T) + 1)

    donor[0] = T_top if face_flux[0] > 0 else T[0]
    donor[-1] = T_bottom if face_flux[-1] < 0 else T[-1]

    donor[1:-1] = np.where(
        face_flux[1:-1] >= 0,
        T[:-1],
        T[1:],
    )
    return donor


def compute_temperature(T, theta, theta_old, theta_sat, RH, q_liq, q_vap, soil, dt, dz,
                               T_atm=303.15, q_evap=0.0, h_c=10.0, S=None):
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
    # J_vap = np.array([vapor_flux_to_mass_flux(qv, t, rh)
      #                for qv, t, rh in zip(q_vap, T, RH)])
    # J_vap = np.insert(J_vap, 0, q_vap[0])

    # N layers have N+1 interfaces.
    q_vap = np.asarray(q_vap, dtype=float)
    T = np.asarray(T, dtype=float)
    RH = np.asarray(RH, dtype=float)

    # Water-equivalent vapor flux [m/s] -> mass flux [kg/m²/s].
    q_vap = np.asarray(q_vap, dtype=float)

    if q_vap.shape != (N + 1,):
        raise ValueError("q_vap must contain N+1 interface fluxes.")

    J_vap = rho_w * q_vap

    # Surface latent cooling is included by the surface energy balance.
    J_latent = J_vap.copy()
    J_latent[0] = 0.0

    # evap flux from [m/s] to [kg/m2/s]
    # evap = np.array(-q_evap) * 100
    evap = rho_w * q_evap

    theta = np.asarray(theta, dtype=float)
    theta_old = np.asarray(theta_old, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)

    # Volumetric heat capacity, neglecting pore-air heat capacity [J/m³/K].
    C = (1.0 - theta_sat) * rho_s * c_s + theta * rho_w * c_w
    # C = (1 - theta) * rho_s * c_s + theta * rho_w * c_w

    C_old =  (1.0 - theta_sat) * rho_s * c_s + theta_old * rho_w * c_w

    dz = np.asarray(dz, dtype=float)
    q_liq = np.asarray(q_liq, dtype=float)

    # Outgoing sensible-heat transport coefficients for each layer.
    outgoing = (
            rho_w * c_w * (
            np.maximum(q_liq[1:], 0.0)
            + np.maximum(-q_liq[:-1], 0.0)
    )
            + 1850.0 * (
                    np.maximum(J_vap[1:], 0.0)
                    + np.maximum(-J_vap[:-1], 0.0)
            )
    )

    active = outgoing > 0.0

    if np.any(active):
        dt_limit = np.min(C_old[active] * dz[active] / outgoing[active])

        if dt > 0.8 * dt_limit:
            raise ValueError(
                f"Explicit heat advection: dt={dt:.1f} s exceeds "
                f"the conservative limit {0.8 * dt_limit:.1f} s."
            )

    lam = lambda_dry + (lambda_sat - lambda_dry) * (theta**1.5)
    L_v = 2.45e6  # latent heat

    T_reference = 273.15

    # Compute surface temperature
    # T_surf = compute_surface_temperature(T[0], dz[0], lam[0], evap, T_atm, h_c, L_v)

    # Build tridiagonal system
    A = np.zeros((N, N))
    RHS = np.zeros(N)

    q_liq = np.asarray(q_liq, dtype=float)

    if q_liq.shape != (N + 1,):
        raise ValueError("q_liq must contain N+1 interface fluxes.")

    T_liq_face = upstream_temperature(q_liq, T, T_atm, T[-1])
    T_vap_face = upstream_temperature(J_vap, T, T_atm, T[-1])

    T_reference = 273.15

    heat_liq = rho_w * c_w * q_liq * (T_liq_face - T_reference)
    heat_vap = 1850.0 * J_vap * (T_vap_face - T_reference)

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
        # RHS[i] = C[i]*T[i]/dt

        RHS[i] = ( C_old[i] * (T[i] - T_reference) + C[i] * T_reference) / dt

        # --- Advection terms with safe upwind indexing ---
        # Liquid
        # T_up_in  = T[i-1] if i > 0 else T[i]       # bottom interface
        # T_up_out = T[i+1] if i < N-1 else T[i]    # top interface
        # adv_liq = rho_w * c_w * (q_liq[i+1]*T_up_out - q_liq[i]*T_up_in)

        # Vapor
        # adv_vap = 1850 * (J_vap[i+1]*T_up_out - J_vap[i]*T_up_in)

        adv_liq = heat_liq[i + 1] - heat_liq[i]
        adv_vap = heat_vap[i + 1] - heat_vap[i]

        # Latent heat from vapor flux
        # latent = L_v * (J_vap[i+1] - J_vap[i])
        # RHS[i] += -(adv_liq + adv_vap)/dz_i + latent/dz_i

        # Divergence of latent heat transport [W/m²].
        # latent = L_v * (J_vap[i + 1] - J_vap[i])

        latent = L_v * (J_latent[i + 1] - J_latent[i])

        RHS[i] -= (adv_liq + adv_vap + latent) / dz_i

        if S is not None:
            RHS[i] += S[i]

    # Surface energy balance, coupled implicitly to the first layer.
    # Temperatures are in Kelvin; evaporation is positive outward.
    g_surface = lam[0] / (0.5 * dz[0])  # W m^-2 K^-1
    denominator = g_surface + h_c

    # Downward conductive heat flux:
    # G = g_surface * (T_surface - T_new[0])
    #   = conductance * (T_atm - T_new[0]) - latent_sink
    conductance = g_surface * h_c / denominator
    latent_sink = g_surface * L_v * evap / denominator

    A[0, 0] += conductance / dz[0]
    RHS[0] += (conductance * T_atm - latent_sink) / dz[0]

    # Bottom boundary: fixed temperature (can be replaced by flux BC later)
    A[-1, :] = 0.0
    A[-1, -1] = 1.0
    RHS[-1] = T[-1]

    T_new = np.linalg.solve(A, RHS)

    # Surface temperature consistent with the updated first-layer temperature.
    T_surf = (
                     g_surface * T_new[0] + h_c * T_atm - L_v * evap
             ) / denominator

    return T_new, T_surf



