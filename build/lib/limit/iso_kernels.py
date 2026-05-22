# iso_kernels.py
# Single source of truth for all storage computations.
# Both scalar methods in iso_storage and the batch kernel call these functions.

import numpy as np
from numba import njit, prange
from math import exp

# ═══════════════════════════════════════════════════════════════════
# ------------------------------------------------------------------#
# ═══════════════════════════════════════════════════════════════════

@njit(cache=True)
def _calc_p_sat(T):
    """Saturation vapour pressure in pascal"""
    Tzero = 273.16000366210938
    return 610.59999465942383 * exp(
        17.270000457763672 * (T - Tzero) /
        ((T - Tzero) + 237.30000305175781)
    )


@njit(cache=True)
def _calc_cv_sat(T):
    """
    Saturated vapor volumetric mass .

    """
    Mw = 0.018015999346971512
    R  = 8.3142995834350586
    return _calc_p_sat(T) * Mw / R / T / 1000.0


@njit(cache=True)
def _calc_relative_humidity(psi, T, S):
    """
    Soil relative humidity.
    """
    Mw      = 0.018015999346971512
    gravity = 9.8000001907348633
    R       = 8.3142995834350586
    rhmin   = 5.0000000745058060E-002

    if S < 1.0:
        return max(exp(Mw * gravity * psi / (R * T)), rhmin)
    else:
        return 1.0


@njit(cache=True)
def _calc_cv(psi, T, T0, S_current):
    """
    Water vapor volumetric mass .
    """
    rh = _calc_relative_humidity(psi, T, S_current)
    return rh * _calc_cv_sat(T)


@njit(cache=True)
def _calc_S(theta, theta_0, theta_sat):
    """Effective saturation — calc_S / get_saturation in iso_storages."""
    return (theta - theta_0) / (theta_sat - theta_0)


@njit(cache=True)
def _calc_liq_saturation(S, cv):
    """Liquid saturation — get_liq_saturation in iso_storages."""
    return (S - cv) / (1.0 - cv)


@njit(cache=True)
def _calc_Sl(s_liq, pond_height, thickness, theta_sat, theta_0, has_pond):
    """
    Liquid saturation with pond — Sl property in iso_storages.
    """
    if has_pond:
        return s_liq + pond_height / thickness / (theta_sat - theta_0)
    return s_liq


@njit(cache=True)
def _calc_alpha_i(T, iso_idx):
    """
    Liquid-vapour isotopic fractionation — alpha_i in iso_storages.
    iso_idx: 0 for 2H, 1 for 18O
    """
    if iso_idx == 0:  # 2H
        return exp(-(24844.0 / T ** 2 + (-76.248) / T + 0.052612))
    else:             # 18O
        return exp(-(1137.0 / T ** 2 + (-0.4156) / T - 0.0020667))


@njit(cache=True)
def _calc_d_alpha(T, iso_idx):
    """
    Derivative of alpha_i w.r.t T — _d_alpha in iso_storages.
    iso_idx: 0 for 2H, 1 for 18O
    """
    if iso_idx == 0:  # 2H
        return (2.0 * 24844.0 / T ** 3 + (-76.248) / T ** 2) / exp(
            24844.0 / T ** 2 + (-76.248) / T + 0.052612)
    else:             # 18O
        return (2.0 * 1137.0 / T ** 3 + (-0.4156) / T ** 2) / exp(
            1137.0 / T ** 2 + (-0.4156) / T - 0.0020667)


@njit(cache=True)
def _calc_beta(alpha_i, d_alpha, dT):
    """beta = alpha_i + d_alpha * dT — beta in iso_storages."""
    return alpha_i + d_alpha * dT


@njit(cache=True)
def _calc_del_eff_S(
    theta_sat, theta_0, theta_t0,
    psi_0, T0,
    S_current, T_current,
    Sl, cv, s_liq, dT,
    alpha_i, d_alpha,
    pond_height_t0, thickness, has_pond
):
    """
    del_eff_saturation in iso_storages.
    S_current: current self.S — used in relative_humidity check for cv_0.
    """
    # S at previous time
    S0 = _calc_S(theta_t0, theta_0, theta_sat)

    # cv at previous time — note: uses current S_current for rh check (matches original)
    cv_0 = _calc_cv(psi_0, T0, T0, S_current)

    # liquid saturation at previous time
    liq_S_0 = _calc_liq_saturation(S0, cv_0)

    # Sl at previous time with pond
    Sl_0 = _calc_Sl(liq_S_0, pond_height_t0, thickness,
                    theta_sat, theta_0, has_pond)

    # deltas
    d_sliq = Sl - Sl_0
    d_cv   = cv - cv_0

    # beta components
    beta  = _calc_beta(alpha_i, d_alpha, dT)
    dbeta = d_alpha * dT

    if s_liq < 1.0:
        return (d_sliq + beta * d_cv + cv * dbeta
                - (s_liq * beta * d_cv
                   + cv * beta * d_sliq
                   + cv * s_liq * dbeta))
    else:
        return d_sliq


@njit(cache=True)
def _calc_eff_S(Sl, cv, s_liq, theta_0, theta_sat, alpha_i, d_alpha, dT):
    """eff_saturation in iso_storages."""
    if s_liq < 1.0:
        beta = _calc_beta(alpha_i, d_alpha, dT)
        return (Sl + cv * beta - cv * s_liq * beta
                + theta_0 / theta_sat)
    else:
        return Sl + theta_0 / theta_sat


@njit(cache=True)
def _calc_avail_vol(theta_sat, thickness, area):
    """get_available_liquid_volume in iso_storages."""
    return theta_sat * thickness * area


@njit(cache=True)
def _calc_eff_liq_vol(eff_S, avail_vol):
    """get_eff_liquid_volume in iso_storages."""
    return eff_S * avail_vol


@njit(cache=True)
def _calc_storage_i(conc_iso, del_eff_S, avail_vol):
    """get_storage_i in iso_storages."""
    return conc_iso * del_eff_S * avail_vol


# ═══════════════════════════════════════════════════════════════════
# ------------------------------------------------------------------#
# ═══════════════════════════════════════════════════════════════════

@njit(cache=True, parallel=True)
def storage_batch(
    # soil state — shape (n,)
    theta_sat, theta_0, theta_t0,
    psi_0, T0, T, dT,
    S_current,      # current self.S for each storage
    Sl, cv, s_liq,
    iso_idx,        # scalar: 0 for 2H, 1 for 18O
    conc_iso,
    area, thickness,
    pond_height_t0,
    pond_height,    # current pond height
    has_pond,
    # outputs — shape (n,)
    out_del_eff_S,
    out_eff_S,
    out_eff_liq_vol,
    out_storage_i
):
    n = len(theta_sat)

    for i in prange(n):

        alpha_i = _calc_alpha_i(T[i], iso_idx)
        d_alpha = _calc_d_alpha(T[i], iso_idx)

        out_del_eff_S[i] = _calc_del_eff_S(
            theta_sat[i], theta_0[i], theta_t0[i],
            psi_0[i], T0[i],
            S_current[i], T[i],
            Sl[i], cv[i], s_liq[i], dT[i],
            alpha_i, d_alpha,
            pond_height_t0[i], thickness[i], has_pond[i]
        )

        out_eff_S[i] = _calc_eff_S(
            Sl[i], cv[i], s_liq[i],
            theta_0[i], theta_sat[i],
            alpha_i, d_alpha, dT[i]
        )

        avail_vol = _calc_avail_vol(theta_sat[i], thickness[i], area[i])

        out_eff_liq_vol[i] = _calc_eff_liq_vol(out_eff_S[i], avail_vol)

        out_storage_i[i] = _calc_storage_i(
            conc_iso[i], out_del_eff_S[i], avail_vol
        )

