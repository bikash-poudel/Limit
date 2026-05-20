

import numpy as np


# ----------------------------------------------------------------------
# 1. CONSTANT AND ANALYTIC PROFILE PARAMETERS (From Table 5)
# ----------------------------------------------------------------------

# Parameters for H2(18)O
ALPHA_IK = 1.0285  # alpha_iK for H2(18)O, calculated as Dv/Dvi  (Merlivat 1978a value)
DELTA_ALIM = 0.0  # delta_i_alim^l (per mille)
DELTA_A_V = -14.0  # delta_i_a^v (per mille)

# General Constants (from paper text)
G = 9.81  # Acceleration of gravity (m/s^2)  (Assumed standard value)
R_VAPOR = 461.5  # Perfect gas mass constant for water vapour (J/kg/K)
RHO_W = 1000.0  # Density of liquid water (kg/m^3) (Assumed standard value)
EVAP_RATE_E = 0.982e-5  # Steady state evaporation rate E (kg/m^2/s)
D_V_REF = 2.82e-5  # Ordinary water vapor diffusivity (m^2/s) (Approx. value at 303 K)


# Temperature Profile from Table 5
def T_profile_K(z: float) -> float:
    """T(z) = 20(1 + exp(-20z)) in Celsius, converted to Kelvin."""
    T_C = 20.0 * (1.0 + np.exp(-20.0 * z))
    return T_C + 273.15  #


# ----------------------------------------------------------------------
# 2. THERMODYNAMIC FUNCTIONS (Depth-dependent, based on T(z))
# ----------------------------------------------------------------------

def alpha_i_star(T_K: float) -> float:
    """Equilibrium fractionation factor alpha_i* for H2(18)O (Majoube, 1971)"""
    T_2_pow = T_K * T_K
    T_1_pow = T_K
    # 1000 * ln(alpha_i*) = 1.137e6/T^2 - 0.4156e3/T - 2.0667 (Formula from Appendix C)
    ln_alpha_1000 = (1.137e6 / T_2_pow) - (0.4156e3 / T_1_pow) - 2.0667
    return np.exp(ln_alpha_1000 / 1000.0)  # [cite: 558]


def rho_v_sat(T_K: float) -> float:
    """Saturated water vapor volumetric mass rho_v_sat (kg/m^3)"""
    # Using the Tetens formula for saturated vapor pressure P_sat (Pa) as approximation
    # P_sat = 610.78 * exp(17.27 * T_C / (T_C + 237.3))
    T_C = T_K - 273.15
    P_sat = 610.78 * np.exp(17.27 * T_C / (T_C + 237.3))
    # rho_v_sat = P_sat / (R_VAPOR * T_K)
    return P_sat / (R_VAPOR * T_K)  # (Used for R_VAPOR)


def z_v_factor(T_K: float) -> float:
    """Vapor transport factor z_v (Eq. 45)"""
    # D_i^v = D_v / alpha_iK
    # D_v is also a function of T, using D_v ~ T^1.75. We approximate it for simplicity.
    # z_v = [ (Dv / alpha_iK) * rho_v_sat ] / (rho_w * E)
    D_i_v = D_V_REF / ALPHA_IK
    rho_v_sat_val = rho_v_sat(T_K)
    return (D_i_v * rho_v_sat_val) / (RHO_W * EVAP_RATE_E)  # (Used for zv formula)


# ----------------------------------------------------------------------
# 3. PLACEHOLDER/MISSING SOIL FUNCTIONS (Need full SiSPAT solution)
# ----------------------------------------------------------------------
# NOTE: These functions require solving the heat/water equations from Appendix B
# and using the Yolo Light Clay hydraulic parameters from Table 1 (not provided here).

def liquid_water_flux_C(z: float) -> float:
    """Placeholder for Liquid Water Flux C(z) = q_l (m/s)"""
    # C(z) must be solved from Richards' equation and is highly variable with depth.
    # Below the evaporation front, it typically tends toward the drainage flux (if any).
    # We must use C(z) from a solved profile for a meaningful result.
    return 0.0  # Placeholder: Assuming C is near-zero below a shallow surface layer


def water_pressure_head_h(z: float) -> float:
    """Placeholder for Water Pressure Head h(z) (m)"""
    # h(z) must be solved from Richards' equation.
    # It is used to calculate h_u(z) and D_i_l(z)
    return -10.0 * np.exp(-10.0 * z)  # Placeholder: Example non-saturated profile (in meters)


def total_liquid_diffusivity_D_i_l(z: float) -> float:
    """Placeholder for Total Liquid Diffusivity D_i_l (m^2/s)"""
    # D_i_l = D_i^lo * tau * theta_l + Lambda*|q_l|. Requires theta_l(h) and q_l(z)
    return 1.0e-9  # Placeholder: A typical value for liquid diffusion


# ----------------------------------------------------------------------
# 4. MAIN DIFFERENTIAL EQUATION (Equation 43)
# ----------------------------------------------------------------------

def d_delta_dz_equation_43(z: float, delta_i_l: float) -> float:
    """
    Calculates the derivative d(delta_i_l)/dz for the non-saturated case (Eq. 43).
    """
    # 1. Calculate Required Depth-Dependent Values
    T_K = T_profile_K(z)
    alpha_star = alpha_i_star(T_K)
    D_i_l = total_liquid_diffusivity_D_i_l(z)

    # Relative humidity from Kelvin Law (Eq. 46): h_u = exp(g * h / (R * T))
    h_val = water_pressure_head_h(z)
    h_u = np.exp(G * h_val / (R_VAPOR * T_K))  #

    # Vapor transport factor
    z_v = z_v_factor(T_K)

    # Liquid Water Flux
    C = liquid_water_flux_C(z)  # C is q_l

    if D_i_l <= 0:
        return np.nan  # Avoid division by zero in dry layers

    # 2. Equation Terms

    # Term 1: [ C * (delta_i^l - delta_i_alim^l) ] / D_i^l
    term1 = (C * (delta_i_l - DELTA_ALIM)) / D_i_l

    # Term 2: Vapor/Fractionation Term
    # Square brackets content: [ A * (delta_i^l - delta_i_alim^l) - B ]
    term2_bracket_A = (ALPHA_IK / alpha_star) * (delta_i_l - DELTA_ALIM)
    term2_bracket_B = (DELTA_ALIM - DELTA_A_V)
    term2_numerator = h_u * z_v * (term2_bracket_A - term2_bracket_B)

    # Term 2 denominator: D_i_l * (1 + h_u * z_v)
    term2_denominator = D_i_l * (1.0 + h_u * z_v)
    term2 = term2_numerator / term2_denominator

    # 3. Final Derivative: d(delta)/dz = Term 1 - Term 2
    d_delta_dz = term1 - term2

    return d_delta_dz


# ----------------------------------------------------------------------
# 5. DEMONSTRATION (Single point calculation using placeholders)
# ----------------------------------------------------------------------

# We must assume an isotopic value for delta_i_l at a specific depth 'z'
# Let's test at z=0.01m, assuming a high delta_i_l (19.92 per mille, the maximum)

z_test = 0.01  # m
# Maximum delta_18O concentration from the theoretical result
delta_i_l_test = 19.92

derivative_value = d_delta_dz_equation_43(
    z=z_test,
    delta_i_l=delta_i_l_test
)

# ----------------------------------------------------------------------
# 6. RESULTS
# ----------------------------------------------------------------------

# Print out the key calculated thermodynamic values for context
T_K_test = T_profile_K(z_test)
alpha_star_test = alpha_i_star(T_K_test)
rho_v_sat_test = rho_v_sat(T_K_test)
z_v_test = z_v_factor(T_K_test)
h_u_test = np.exp(G * water_pressure_head_h(z_test) / (R_VAPOR * T_K_test))

print(f"--- Thermodynamic and Transport Parameters at z={z_test} m ---")
print(f"Temperature T(z): {T_K_test:.2f} K ({T_K_test - 273.15:.2f} °C)")
print(f"Equilibrium Factor alpha_i*(T): {alpha_star_test:.4f}")
print(f"Relative Humidity h_u(z) (Placeholder used): {h_u_test:.4f}")
print(f"Vapor Transport Factor z_v(T): {z_v_test:.4f}")
print(f"Liquid Water Flux C(z) (Placeholder used): {liquid_water_flux_C(z_test):.2e} m/s")
print("\n--- Equation (43) Derivative Result ---")
print(f"d(delta_18O)/dz at z={z_test}m and delta_l={delta_i_l_test}‰: {derivative_value:.2f} ‰/m")