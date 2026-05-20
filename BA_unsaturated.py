import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import BA_analytical_unsaturated as BA

# ------------------------------------------------------------------
# 1. CONSTANTS
# ------------------------------------------------------------------
ISOTOPE = '2H'  # Running for HDO as requested
E_s = 1.043e-5
delta_inf = 0.0
n_k = 1.0
rho_w = 1000.0
tau_w = 0.67
tau_g = 0.67
n_soil = 0.35
T_a = 313.15  # Air temperature (K)
h_a = 0.2  # Air relative humidity
g = 9.81
M_w = 0.018015
R_u = 8.314

# Isotope-specific constants for HDO ('2H')
ratio_D_v = 1.0166
sigma_i = 0.98331
delta_ia_v = -100.0  # Atmospheric isotope ratio (per mil)
R_std = 155.76e-6


# ------------------------------------------------------------------
# 2. HELPER FUNCTIONS (PHYSICAL EQUATIONS)
# ------------------------------------------------------------------

def calc_alpha_star(T):
    if ISOTOPE == '18O':
        return np.exp(2.0667e-3 + 0.4156 / T - 1.137e3 / (T ** 2))
    else:
        return np.exp(-52.612e-3 + 76.248 / T - 24.844e3 / (T ** 2))


def calc_alpha_k():
    return ratio_D_v ** n_k


def calc_rho_sat_v(T):
    return 1e-3 * np.exp(31.3716 - 6014.79 / T - 7.92495e-3 * T) / T


def calc_Hr(h, T):
    return np.exp(M_w * g * h / (R_u * T))


def calc_D_li(T):
    return sigma_i * 1e-9 * np.exp(-535400 / (T ** 2) + 1393.3 / T + 2.1876)


def calc_D_vi(T):
    return calc_D_v(T) / ratio_D_v


def calc_D_v(T):
    return 2.12e-5 * (T / 273.16) ** 2


def calc_D_l_star(D_li, theta_l):
    return D_li * tau_w * theta_l


def calc_D_v_star(D_vi, theta_l):
    theta_g = n_soil - theta_l
    return theta_g * tau_g * D_vi


def calc_z_l(D_l_star):
    return rho_w * D_l_star / E_s


def calc_z_v(D_v_star, rho_sat_v):
    return D_v_star * rho_sat_v / E_s


def calculate_surface_delta(T_s, H_rs):
    alpha_star_s = calc_alpha_star(T_s)
    alpha_k_s = calc_alpha_k()
    rho_sat_v_s = calc_rho_sat_v(T_s)
    rho_sat_v_a = calc_rho_sat_v(T_a)
    h_a_prime = h_a * (rho_sat_v_a / rho_sat_v_s)
    R_Rstd_inf = (delta_inf / 1000.0) + 1.0
    R_Rstd_a = (delta_ia_v / 1000.0) + 1.0

    numerator = (R_Rstd_inf * (H_rs - h_a_prime) * alpha_k_s) + (h_a_prime * R_Rstd_a)
    denominator = H_rs * alpha_star_s

    R_Rstd_L = numerator / denominator
    delta_0 = (R_Rstd_L - 1.0) * 1000.0

    return delta_0


# ------------------------------------------------------------------
# 3. ODE MODEL DEFINITION (Eq. 64)
# ------------------------------------------------------------------

def ode_model(z, delta_l, interp_funcs):
    """
    Defines the differential equation: d(delta_l)/dz = f(z, delta_l).
    Uses interpolated functions for all depth-dependent coefficients.
    """
    delta_l = delta_l[0]

    # --- Fetch coefficients at current depth z ---
    alpha_star = interp_funcs['alpha_star'](z)
    alpha_k = interp_funcs['alpha_k'](z)
    Hr = interp_funcs['Hr'](z)
    z_l = interp_funcs['z_l'](z)
    z_v = interp_funcs['z_v'](z)
    # The derivative term is now a function of z from the interpolation
    d_log_term_dz = interp_funcs['d_log_term_dz'](z)

    # --- Calculate d(delta_l)/dz ---
    denom = (z_l + Hr * z_v)

    # Check for near-zero or non-physical denominator
    if denom <= 1e-12:
        return [0.0]  # If flow is negligible, gradient is zero

    # Term A (Driving term)
    term_A = (Hr * z_v / denom) * (alpha_k - alpha_star) * d_log_term_dz

    # Term B (Restoring term)
    term_B = (1.0 / denom) * (delta_l - delta_inf)

    d_delta_dz = term_A - term_B

    return [d_delta_dz]


# ------------------------------------------------------------------
# 4. MAIN SOLVER FUNCTION
# ------------------------------------------------------------------

def solve_continuous_isotope_profile(z_nodes, T_profile, h_profile, theta_l_profile):
    # --- 1. Pre-calculate all depth-dependent coefficients ---
    alpha_star_profile = calc_alpha_star(T_profile)
    alpha_k_profile = calc_alpha_k()
    rho_sat_v_profile = calc_rho_sat_v(T_profile)
    Hr_profile = np.array([calc_Hr(h, T) for h, T in zip(h_profile, T_profile)])

    D_li_profile = calc_D_li(T_profile)
    D_vi_profile = calc_D_vi(T_profile)
    D_l_star_profile = calc_D_l_star(D_li_profile, theta_l_profile)
    D_v_star_profile = calc_D_v_star(D_vi_profile, theta_l_profile)
    z_l_profile = calc_z_l(D_l_star_profile)
    z_v_profile = calc_z_v(D_v_star_profile, rho_sat_v_profile)

    # Calculate the derivative term (CRITICAL STEP)
    # Note: np.gradient is only used here on the input profile nodes, then interpolated.
    log_term_profile = np.log(Hr_profile * rho_sat_v_profile * (alpha_k_profile - alpha_star_profile))
    d_log_term_dz_profile = np.gradient(log_term_profile, z_nodes)

    # --- 2. Create interpolation functions ---
    interp_funcs = {
        'alpha_star': interp1d(z_nodes, alpha_star_profile, fill_value="extrapolate"),
        'alpha_k': interp1d(z_nodes, np.full_like(z_nodes, alpha_k_profile), fill_value="extrapolate"),
        'Hr': interp1d(z_nodes, Hr_profile, fill_value="extrapolate"),
        'z_l': interp1d(z_nodes, z_l_profile, fill_value="extrapolate"),
        'z_v': interp1d(z_nodes, z_v_profile, fill_value="extrapolate"),
        'd_log_term_dz': interp1d(z_nodes, d_log_term_dz_profile, fill_value="extrapolate")
    }

    # --- 3. Calculate the surface boundary condition (Initial Value) ---
    T_s = T_profile[0]
    H_rs = Hr_profile[0]
    delta_0_initial = calculate_surface_delta(T_s, H_rs)

    print(f"Calculated surface boundary condition (delta_0): {delta_0_initial:.2f} ‰")

    # --- 4. Solve the ODE using the stable solver ---
    z_span = [z_nodes.min(), z_nodes.max()]

    # Use method='Radau' or 'BDF' for very stiff problems, but RK45 often works for this.
    sol = solve_ivp(
        ode_model,
        z_span,
        [delta_0_initial],
        args=(interp_funcs,),
        dense_output=True,
        method='RK45'
    )

    # --- 5. Return the solution ---
    z_sol = np.linspace(z_span[0], z_span[1], 200)
    delta_l_sol = sol.sol(z_sol)[0]

    return z_sol, delta_l_sol


# ------------------------------------------------------------------
# 5. EXAMPLE USAGE (using the same dry zone input profiles)
# ------------------------------------------------------------------

if __name__ == "__main__":

    print(f"Running CONTINUOUS simulation for {ISOTOPE} (Deuterium).")

    # --- Input Data (Use your known good model inputs here) ---
    cell, dz, ds, d = BA.setup_run(dt=1, sim_period=250)
    layers = cell.layers
    z_nodes = np.array([l.center.z for l in layers])  # np.linspace(0, 1.0, 101)

    z_nodes = np.linspace(0, 1.0, 101)  # 1m profile, 101 equidistant nodes
    T_profile = 20.0 * (1.0 + np.exp(-20.0 * z_nodes)) + 273.15

    # Using the highly variable, dry profiles necessary for the Figure 4 shape
    h_profile = -10.0 * np.exp(-50 * z_nodes) - 1.0
    h_profile[z_nodes > 0.05] = -100.0
    h_profile[z_nodes > 0.2] = -5.0
    h_profile[z_nodes > 0.5] = -0.01
    h_profile[0] = -10.0

    theta_l_profile = 0.05 * np.ones_like(z_nodes)
    theta_l_profile[z_nodes > 0.02] = 0.05 + 0.30 * (1 - np.exp(-50 * (z_nodes[z_nodes > 0.02] - 0.02)))
    theta_l_profile[theta_l_profile > n_soil] = n_soil
    theta_l_profile[0] = 0.05
    # ----------------------------------------------------------

    # --- Call the continuous solver function ---
    try:
        z_solution, delta_l_solution = solve_continuous_isotope_profile(
            z_nodes,
            T_profile,
            h_profile,
            theta_l_profile
        )

        print(f"Simulation complete.")
        print(f"Isotope value at bottom (z={z_solution[-1]:.2f}m): {delta_l_solution[-1]:.2f} ‰")
        
        # --- Plot the results ---
        plt.figure(figsize=(6, 8))
        plt.plot(delta_l_solution, z_solution, 'b-', label='Continuous Solution (RK45)')
        plt.gca().invert_yaxis()
        plt.xlabel(f'Soil Water $\delta^{ISOTOPE}$ (‰)')
        plt.ylabel('Depth (m)')
        plt.title(f'Robust Solution of Semi-Analytical Isotope Equation ($\delta^{ISOTOPE}$)')
        plt.legend()
        plt.grid(True)
        plt.show()

    except Exception as e:
        print(f"\nAn error occurred during continuous simulation: {e}")
        print("If the error persists, check your input profiles for non-monotonicity or large jumps.")

