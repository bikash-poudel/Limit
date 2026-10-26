
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss

import iso_cmf as iso
from limit import iso_fluxes as isof


def solve_unsaturated(cell, n_reference=1000):
    """Integrate Braud Eq. (43) with an exchange-derived liquid surface value.

The boundary value applies at z=0, not at the first-layer centre.
Compare numerical layer concentrations with reference layer means.
Reads the completed numerical state without changing concentrations.
"""
    layers = cell.layers
    ev = cell.connection_evap
    kwargs = iso.iso_delta.test_case_args(testcase=6, BA=True)

    upper_raw = np.array([l.upper_boundary for l in layers])
    lower_raw = np.array([l.lower_boundary for l in layers])

    direction = np.sign(lower_raw[-1] - upper_raw[0])
    upper = direction * (upper_raw - upper_raw[0])
    lower = direction * (lower_raw - upper_raw[0])
    z = 0.5 * (upper + lower)

    if direction == 0 or np.any(lower <= upper):
        raise ValueError("Invalid layer geometry.")

    if not isinstance(n_reference, (int, np.integer)) or n_reference < 2:
        raise ValueError("n_reference must be an integer of at least 2.")

    reference_grid = np.linspace(
        0.0, lower[-1], n_reference + 1
    )
    reference_step = lower[-1] / n_reference

    E = float(cell.q_evap) / float(cell.area)
    if E <= 0:
        raise ValueError("This reference requires positive evaporation.")

    liquid = isof.liquid_diffusion_base_class()
    vapor = isof.vapor_diffusion_base_class()

    # Liquid-water-equivalent vapor content [m³/m³].
    cv = np.array([l.cv for l in layers], dtype=float)

    output = Path("analytical_unsaturated")
    output.mkdir(parents=True, exist_ok=True)

    results = {}
    fig, axes = plt.subplots(
        1, 2, figsize=(8, 4.5), sharey=True,
        constrained_layout=True,
    )

    for ax, isotope in zip(axes, ["2H", "18O"]):
        Dl = np.array([
            liquid.dl_i(
                T=l.T,
                Isotopologue=isotope,
                theta=l.theta,
                tortuosity=l.tortuosity,
                formulation="Melayah",  # Braud Appendix C, Eq. (C.1)
                **kwargs,
            )
            for l in layers
        ])

        Dv = np.array([
            vapor.dv_soil_air(
                T=l.T,
                theta=l.theta,
                theta_sat=l.theta_sat,
                tortuosity=l.tortuosity,
            )
            for l in layers
        ])

        Dvi = np.array([
            vapor.dv_i(dv=d, Isotopologue=isotope, **kwargs)
            for d in Dv
        ])

        alpha = np.array([
            l.alpha_i(Isotopologue=isotope, **kwargs)
            for l in layers
        ])

        # BA84 pure molecular kinetic factor, nK = 1.
        dv_air = vapor.dv_free_air(
            T=ev.T_surface, Pa=cell.atmosphere.Pa
        )
        alpha_k = dv_air / vapor.dv_i(
            dv=dv_air, Isotopologue=isotope, **kwargs
        )

        h_exchange = float(cell.atmosphere.cv_a / ev.cv_surface)

        alpha_surface = layers[0].alpha_i(
            Isotopologue=isotope,
            T=ev.T_surface,
            **kwargs,
        )

        delta_feed = 0.0
        delta_vapor = {"2H": -100.0, "18O": -14.0}[isotope]

        delta_surface = 1000.0 * (
                (
                        (1.0 - h_exchange)
                        * alpha_k
                        * (1.0 + delta_feed / 1000.0)
                        + h_exchange * (1.0 + delta_vapor / 1000.0)
                ) / alpha_surface - 1.0
        )

        """Braud Eq. (43), with an exchange-derived liquid interface value.

        Reads the completed numerical state without changing concentrations.
        Depth is positive downward; diffusivities are in m²/s.
        """

        # Extend coefficient profiles to the physical boundaries.
        numerical = np.asarray(
            cell.conc_2H_delta if isotope == "2H"
            else cell.conc_18O_delta,
            dtype=float,
        )

        """
        # Surface-conditioned reference.    
        # Extend vapor content with zero gradient across the upper half-layer.
        zz = np.r_[0.0, z, lower[-1]]

        # cv_nodes = np.r_[ev.cv_surface, cv, cv[-1]]
        # Zero surface vapor flux: zero vapor-content gradient
        cv_nodes = np.r_[cv[0], cv, cv[-1]]
        
        alpha_nodes = np.r_[alpha_surface, alpha, alpha[-1]]
        """

        """
        # Extend coefficient profiles to the surface and bottom boundaries.
        zz = np.r_[0.0, z, lower[-1]]
        # cv_nodes = np.r_[cv[0], cv, cv[-1]]
        # cv_nodes = np.r_[ev.cv_surface, cv, cv[-1]]
        cv_nodes = np.r_[ev.cv_surface, cv, cv[-1]]
        alpha_nodes = np.r_[alpha_surface, alpha, alpha[-1]]

        b = PchipInterpolator(
            zz, cv_nodes * (alpha_k - alpha_nodes)
        )
        cv_interp = PchipInterpolator(zz, cv_nodes)

        def rhs(depth, value):
            dl = np.interp(depth, z, Dl)
            dvi = np.interp(depth, z, Dvi)
            denominator = dl + dvi * cv_interp(depth)

            if denominator <= 0:
                raise ValueError("Nonpositive reference transport coefficient.")

            P = E / denominator
            Q = 1000.0 * dvi * b.derivative()(depth) / denominator

            # Feed delta = 0‰ for the current numerical configuration.
            return [-P * (value[0] - delta_feed) + Q]
        
        """

        # Exact prescribed temperature profile.
        def temperature_at(depth):
            return 293.15 + 20.0 * np.exp(-20.0 * depth)

        temperature_nodes = np.r_[
            float(ev.T_surface),
            np.array([layer.T for layer in layers], dtype=float),
            float(layers[-1].T),
        ]
        temperature_depths = np.r_[0.0, z, lower[-1]]

        temperature_interp = PchipInterpolator(
            temperature_depths, temperature_nodes
        )
        temperature_gradient = temperature_interp.derivative()

        # Saturated liquid-water-equivalent vapor content and its T derivative.
        # Uses the saturation-pressure expression supplied earlier.
        def saturated_vapor_at(T):
            Tzero = 273.16000366210938
            a = 17.270000457763672
            b = 237.30000305175781
            Mw = 0.018015999346971512
            Rgas = 8.3142995834350586

            offset = T - Tzero
            pressure = 610.59999465942383 * np.exp(
                a * offset / (offset + b)
            )

            cv_sat = pressure * Mw / (Rgas * T * 1000.0)

            dlog_cv_sat_dT = (
                    a * b / (offset + b) ** 2 - 1.0 / T
            )

            return cv_sat, dlog_cv_sat_dT

        # Use the same physical surface temperature as the prescribed profile.
        # T_surface_reference = temperature_at(0.0)
        T_surface_reference = float(ev.T_surface)


        # Infer a surface head consistent with the supplied interface vapor content.
        g = 9.81
        Rv = 461.5

        cv_sat_surface, _ = saturated_vapor_at(T_surface_reference)
        humidity_surface = float(ev.cv_surface) / cv_sat_surface

        if not 0.0 < humidity_surface <= 1.0:
            raise ValueError(
                "Interface vapor content implies relative humidity outside (0, 1]."
            )

        head_surface = (
                Rv * T_surface_reference / g * np.log(humidity_surface)
        )

        head = np.array([layer.psi for layer in layers], dtype=float)
        theta = np.array([layer.theta for layer in layers], dtype=float)
        porosity = np.array(
            [layer.theta_sat for layer in layers], dtype=float
        )
        tortuosity = np.array(
            [layer.tortuosity for layer in layers], dtype=float
        )

        zz = np.r_[0.0, z, lower[-1]]

        # Reconstruction assumptions:
        # - Surface head reproduces ev.cv_surface through Kelvin's law.
        # - Surface water content uses the first-layer representative value.
        # - Bottom values extend the last-layer values.
        head_at = PchipInterpolator(
            zz, np.r_[head_surface, head, head[-1]]
        )
        head_gradient_at = head_at.derivative()

        theta_at = PchipInterpolator(
            zz, np.r_[theta[0], theta, theta[-1]]
        )

        def rhs(depth, value):
            # Exact prescribed temperature and spatial derivative.
            # T_here = temperature_at(depth)
            # dT_dz = -400.0 * np.exp(-20.0 * depth)

            T_here = float(temperature_interp(depth))
            dT_dz = float(temperature_gradient(depth))

            head_here = float(head_at(depth))
            dh_dz = float(head_gradient_at(depth))

            theta_here = float(theta_at(depth))
            porosity_here = float(np.interp(depth, z, porosity))
            tortuosity_here = float(np.interp(depth, z, tortuosity))

            # Kelvin humidity and vapor content evaluated at this depth.
            cv_sat_here, dlog_cv_sat_dT = saturated_vapor_at(T_here)

            humidity_here = np.exp(
                g * head_here / (Rv * T_here)
            )
            cv_here = humidity_here * cv_sat_here

            # Chain-rule derivative of vapor content.
            dlog_cv_dz = (
                    g / Rv * (
                    dh_dz / T_here
                    - head_here * dT_dz / T_here ** 2
            )
                    + dlog_cv_sat_dT * dT_dz
            )
            dcv_dz = cv_here * dlog_cv_dz

            # Evaluate fractionation using LIMIT's temperature formulation.
            alpha_here = layers[0].alpha_i(
                Isotopologue=isotope,
                T=T_here,
                **kwargs,
            )

            # Small centred temperature difference for d(alpha)/dT.
            # This differentiates the fractionation function, not layer profiles.
            epsilon_T = 1e-3
            dalpha_dT = (
                                layers[0].alpha_i(
                                    Isotopologue=isotope,
                                    T=T_here + epsilon_T,
                                    **kwargs,
                                )
                                - layers[0].alpha_i(
                            Isotopologue=isotope,
                            T=T_here - epsilon_T,
                            **kwargs,
                        )
                        ) / (2.0 * epsilon_T)

            dalpha_dz = dalpha_dT * dT_dz

            # Evaluate diffusivities at the requested depth.
            dl_here = liquid.dl_i(
                T=T_here,
                Isotopologue=isotope,
                theta=theta_here,
                tortuosity=tortuosity_here,
                formulation="Melayah",
                **kwargs,
            )

            dv_here = vapor.dv_soil_air(
                T=T_here,
                theta=theta_here,
                theta_sat=porosity_here,
                tortuosity=tortuosity_here,
            )

            dvi_here = vapor.dv_i(
                dv=dv_here,
                Isotopologue=isotope,
                **kwargs,
            )

            denominator = dl_here + dvi_here * cv_here
            if not np.isfinite(denominator) or denominator <= 0:
                raise ValueError("Invalid reference transport coefficient.")

            # d[cv * (alpha_k - alpha)] / dz
            product_gradient = (
                    (alpha_k - alpha_here) * dcv_dz
                    - cv_here * dalpha_dz
            )

            P = E / denominator
            Q = 1000.0 * dvi_here * product_gradient / denominator

            return [-P * (value[0] - delta_feed) + Q]

        solution = solve_ivp(
            rhs,
            (0.0, lower[-1]),
            [delta_surface],
            method="Radau",
            rtol=1e-8,
            atol=1e-9,
            dense_output=True,
            t_eval=reference_grid,
            max_step=reference_step,
        )

        if not solution.success:
            raise RuntimeError(solution.message)

        nodes, weights = leggauss(8)

        points = (
                z[:, None]
                + 0.5 * (lower - upper)[:, None] * nodes
        )

        reference_mean = (
                0.5
                * solution.sol(points.ravel())[0].reshape(points.shape)
                @ weights
        )

        # dense_z = np.linspace(0.0, lower[-1], 3001)
        dense_z = reference_grid
        centers_z = dense_z - np.diff(dense_z, prepend=0) / 2
        # centers_z = np.array([(a + b) / 2 for a, b in zip([0] + dense_z[:-1], dense_z)])

        analytical = solution.sol(centers_z)[0][1:]

        difference = numerical - reference_mean
        print(
            f"{isotope}: exchange-conditioned surface reference; "
            f"h={h_exchange:.6g}; "
            f"surface delta={delta_surface:.6g}‰; "
            f"RMSE={np.sqrt(np.mean(difference ** 2)):.6g}‰; "
            f"maximum difference={np.max(np.abs(difference)):.6g}‰"
        )

        # Analytical solution: blue dashed line.
        ax.plot(
            #solution.sol(dense_z)[0][:-1],
            # -dense_z,
            analytical,
            -z,
            linestyle="--",
            color="tab:blue",
            linewidth=1.4,
            label="Analytical",
        )

        # LIMIT numerical solution: red solid line.
        ax.plot(
            numerical,
            -z,
            linestyle="-",
            color="tab:red",
            linewidth=1.4,
            label="Numerical",
        )

        ax.set_xlabel(
            r"$\delta^{2}\mathrm{H}$ (‰)"
            if isotope == "2H"
            else r"$\delta^{18}\mathrm{O}$ (‰)",
            fontsize=10,
        )

        # Show up to the first 0.5 m, or the full column if shallower.
        ax.set_ylim(-min(0.5, lower[-1]), 0.0)

        ax.tick_params(
            axis="both",
            labelsize=9,
            direction="out",
            length=3,
            width=0.8,
            top=False,
            right=False,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        # Place panel labels inside the upper-left corner.
        panel_label = "(a)" if isotope == "2H" else "(b)"
        ax.text(
            -0.16, 1.02,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            clip_on=False,
        )

        # Consistent legend position.
        ax.legend(
            loc="lower right",
            frameon=False,
            fontsize=8,
            handlelength=2.5,
            handletextpad=0.6,
        )

        results[isotope] = {
            "depth": z.copy(),
            "depth_": centers_z.copy(),
            "numerical": numerical.copy(),
            "analytical": analytical.copy(),
            "reference_mean": reference_mean,
            "difference": difference,
            "surface_delta": delta_surface,
        }

    axes[0].set_ylabel("Depth [m]", fontsize=10)

    fig.savefig(
        output / "unsaturated_comparison.pdf",
        bbox_inches="tight",
    )

    fig.savefig(
        output / "unsaturated_comparison.png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.show()
    plt.close(fig)

    return results


def semi_analytical(
    z, T_z, h_z, D_l, D_v, q_evap,
    alpha_i, alpha_ik, delta_surface, delta_alim,
    subdivisions=32,
):
    """Braud Eq. (43) on a refined reference grid.

    D_v is effective isotope vapor diffusivity [m²/s].
    Evaporation is water-equivalent [m/s]; delta is in per mil.
    Physical input profiles are reconstructed linearly between nodes.
    """
    z = np.asarray(z, dtype=float)
    profiles = [
        np.asarray(a, dtype=float)
        for a in (T_z, h_z, D_l, D_v, alpha_i)
    ]
    T_z, h_z, D_l, D_v, alpha_i = profiles

    if (
        z.ndim != 1
        or len(z) < 2
        or np.any(np.diff(z) <= 0)
        or any(a.shape != z.shape for a in profiles)
        or not np.all(np.isfinite(np.concatenate([z] + profiles)))
    ):
        raise ValueError("Need finite, matching profiles at increasing depths.")

    if q_evap <= 0 or subdivisions < 1:
        raise ValueError("Positive evaporation and subdivisions required.")

    rho_w = 1000.0
    g = 9.81
    R = 461.5

    output_z = [float(z[0])]
    output_delta = [float(delta_surface)]

    for i, length in enumerate(np.diff(z)):
        # Fixed slopes of the reconstructed physical profiles.
        dT_dz = (T_z[i + 1] - T_z[i]) / length
        dh_dz = (h_z[i + 1] - h_z[i]) / length
        dalpha_dz = (alpha_i[i + 1] - alpha_i[i]) / length

        edges = np.linspace(z[i], z[i + 1], subdivisions + 1)

        for left, right in zip(edges[:-1], edges[1:]):
            midpoint = 0.5 * (left + right)
            fraction = (midpoint - z[i]) / length
            step = right - left

            T = T_z[i] + fraction * (T_z[i + 1] - T_z[i])
            head = h_z[i] + fraction * (h_z[i + 1] - h_z[i])
            dl = D_l[i] + fraction * (D_l[i + 1] - D_l[i])
            dv = D_v[i] + fraction * (D_v[i + 1] - D_v[i])
            alpha = (
                alpha_i[i]
                + fraction * (alpha_i[i + 1] - alpha_i[i])
            )

            humidity = np.exp(g * head / (R * T))
            rho_vsat = (
                np.exp(31.3716 - 6014.79 / T - 0.00792495 * T)
                / T * 1e-3
            )
            cv = humidity * rho_vsat / rho_w

            # Exact chain-rule derivative for this reconstruction.
            dlog_rho_sat_dT = (
                6014.79 / T**2 - 0.00792495 - 1.0 / T
            )
            dlog_cv_dz = (
                g / R * (dh_dz / T - head * dT_dz / T**2)
                + dlog_rho_sat_dT * dT_dz
            )
            dcv_dz = cv * dlog_cv_dz

            product_gradient = (
                (alpha_ik - alpha) * dcv_dz
                - cv * dalpha_dz
            )

            transport = dl + dv * cv
            if not np.isfinite(transport) or transport <= 0:
                raise ValueError("Invalid reference transport coefficient.")

            P = q_evap / transport
            Q = (
                1000.0 * dv * product_gradient / transport
                + P * delta_alim
            )

            # Stable update with coefficients frozen at the midpoint.
            decay = np.exp(-P * step)
            value = (
                decay * output_delta[-1]
                + (-np.expm1(-P * step)) * Q / P
            )

            output_z.append(float(right))
            output_delta.append(float(value))

    return np.asarray(output_z), np.asarray(output_delta)


def solve_unsaturated_1(cell):
    """Compare LIMIT with the grid-based approximation to Braud Eq. (43).

    Calls semi_analytical(); does not modify the numerical state.
    Delta values are in per mil; depth is positive downward.
    """
    layers = cell.layers
    ev = cell.connection_evap
    kwargs = iso.iso_delta.test_case_args(testcase=6, BA=True)

    upper_raw = np.array(
        [layer.upper_boundary for layer in layers], dtype=float
    )
    lower_raw = np.array(
        [layer.lower_boundary for layer in layers], dtype=float
    )

    surface = upper_raw[0]
    direction = np.sign(lower_raw[-1] - surface)
    upper = direction * (upper_raw - surface)
    lower = direction * (lower_raw - surface)
    centres = 0.5 * (upper + lower)

    if direction == 0 or np.any(lower <= upper):
        raise ValueError("Invalid layer geometry.")

    area = float(cell.area)
    if area <= 0:
        raise ValueError("Cell area must be positive.")

    E = float(cell.q_evap) / area
    if E <= 0:
        raise ValueError("Positive evaporation is required.")

    liquid = isof.liquid_diffusion_base_class()
    vapor = isof.vapor_diffusion_base_class()

    T_surface = float(ev.T_surface)
    cv_surface = float(ev.cv_surface)

    if cv_surface <= 0:
        raise ValueError("Positive surface vapor content is required.")

    h_exchange = 0.2 # float(cell.atmosphere.cv_a) / cv_surface

    # Surface vapor density consistent with semi_analytical().
    rho_vsat_surface = (
        np.exp(
            31.3716
            - 6014.79 / T_surface
            - 0.00792495 * T_surface
        )
        / T_surface * 1e-3
    )

    humidity_surface = cv_surface * 1000.0 / rho_vsat_surface

    if not 0.0 < humidity_surface <= 1.0:
        raise ValueError(
            "Surface vapor content implies humidity outside (0, 1]."
        )

    # Infer surface head through Kelvin's equation.
    head_surface = (
        461.5 * T_surface / 9.81 * np.log(humidity_surface)
    )

    # Grid: physical surface, layer centres, physical bottom.
    z_nodes = np.r_[0.0, centres, lower[-1]]

    T_layers = np.array([layer.T for layer in layers], dtype=float)
    head_layers = np.array(
        [layer.psi for layer in layers], dtype=float
    )

    T_nodes = np.r_[T_surface, T_layers, T_layers[-1]]
    head_nodes = np.r_[
        head_surface, head_layers, head_layers[-1]
    ]

    # Surface theta and tortuosity use first-layer representative values.
    theta_layers = np.array(
        [layer.theta for layer in layers], dtype=float
    )
    porosity_layers = np.array(
        [layer.theta_sat for layer in layers], dtype=float
    )
    tortuosity_layers = np.array(
        [layer.tortuosity for layer in layers], dtype=float
    )

    theta_nodes = np.r_[
        theta_layers[0], theta_layers, theta_layers[-1]
    ]
    porosity_nodes = np.r_[
        porosity_layers[0], porosity_layers, porosity_layers[-1]
    ]
    tortuosity_nodes = np.r_[
        tortuosity_layers[0],
        tortuosity_layers,
        tortuosity_layers[-1],
    ]

    fig, axes = plt.subplots(
        1, 2, figsize=(8, 4.5),
        sharey=True, constrained_layout=True,
    )
    results = {}

    for ax, isotope in zip(axes, ["2H", "18O"]):
        Dl_nodes = np.array([
            liquid.dl_i(
                T=T,
                Isotopologue=isotope,
                theta=theta,
                tortuosity=tau,
                formulation="Melayah",
                **kwargs,
            )
            for T, theta, tau in zip(
                T_nodes, theta_nodes, tortuosity_nodes
            )
        ])

        Dv_nodes = np.array([
            vapor.dv_soil_air(
                T=T,
                theta=theta,
                theta_sat=porosity,
                tortuosity=tau,
            )
            for T, theta, porosity, tau in zip(
                T_nodes, theta_nodes,
                porosity_nodes, tortuosity_nodes,
            )
        ])

        Dvi_nodes = np.array([
            vapor.dv_i(
                dv=Dv, Isotopologue=isotope, **kwargs
            )
            for Dv in Dv_nodes
        ])

        alpha_nodes = np.array([
            layers[0].alpha_i(
                Isotopologue=isotope, T=T, **kwargs
            )
            for T in T_nodes
        ])

        dv_air = vapor.dv_free_air(
            T=T_surface, Pa=cell.atmosphere.Pa
        )
        alpha_k = dv_air / vapor.dv_i(
            dv=dv_air, Isotopologue=isotope, **kwargs
        )

        delta_feed = 0.0
        delta_vapor = {"2H": -100.0, "18O": -14.0}[isotope]

        delta_surface =(
            (
                (1.0 - h_exchange)
                * alpha_k
                * (1.0 + delta_feed)
                + h_exchange * (1.0 + delta_vapor)
            )
            / alpha_nodes[0] - 1.0
        )

        previous_z = None
        previous_delta = None

        for subdivisions in (1, 1):
            reference_z, reference_delta = semi_analytical(
                z=z_nodes,
                T_z=T_nodes,
                h_z=head_nodes,
                D_l=Dl_nodes,
                D_v=Dvi_nodes,
                q_evap=E,
                alpha_i=alpha_nodes,
                alpha_ik=alpha_k,
                delta_surface=delta_surface,
                delta_alim=delta_feed,
                subdivisions=subdivisions,
            )

            if previous_z is not None:
                previous_on_current = np.interp(
                    reference_z, previous_z, previous_delta
                )
                change = np.max(
                    np.abs(reference_delta - previous_on_current)
                )
                print(
                    f"{isotope}: {subdivisions} substeps/interval; "
                    f"maximum reference change={change:.6g}‰"
                )

            previous_z = reference_z
            previous_delta = reference_delta

        # Layer means of the piecewise-linear reference output.
        # This interpolation is only for comparison, not RHS evaluation.
        reference_mean = np.empty(len(layers), dtype=float)

        for j, (top, bottom) in enumerate(zip(upper, lower)):
            inside = (
                (reference_z > top) & (reference_z < bottom)
            )
            sample_z = np.r_[top, reference_z[inside], bottom]
            sample_delta = np.interp(
                sample_z, reference_z, reference_delta
            )

            integral = np.sum(
                0.5
                * (sample_delta[:-1] + sample_delta[1:])
                * np.diff(sample_z)
            )
            reference_mean[j] = integral / (bottom - top)

        numerical = np.asarray(
            cell.conc_2H_delta if isotope == "2H"
            else cell.conc_18O_delta,
            dtype=float,
        )

        difference = numerical - reference_mean

        print(
            f"{isotope}: grid-based Eq. (43); "
            f"surface delta={delta_surface:.6g}‰; "
            f"RMSE={np.sqrt(np.mean(difference**2)):.6g}‰; "
            f"maximum difference="
            f"{np.max(np.abs(difference)):.6g}‰"
        )

        # Analytical reference: dashed line with dots.
        ax.plot(
            reference_delta[1:-1], -centres,
            color="tab:blue",
            linestyle="--",
            marker=".",
            markersize=3,
            linewidth=1.2,
            label="Eq. (43), soil profile",
        )

        # Numerical solution: continuous line.
        ax.plot(
            numerical, -centres,
            color="tab:red",
            linestyle="-",
            linewidth=1.4,
            label="LIMIT",
        )

        ax.set_xlabel(
            r"$\delta^{2}$H (‰)" if isotope == "2H"
            else r"$\delta^{18}$O (‰)"
        )
        ax.set_ylim(-min(0.3, lower[-1]), 0.0)
        ax.legend(fontsize=7, frameon=False)

        results[isotope] = {
            "depth": centres.copy(),
            "numerical": numerical.copy(),
            "reference_mean": reference_mean,
            "difference": difference,
            "surface_delta": delta_surface,
            "reference_depth": reference_z.copy(),
            "reference_delta": reference_delta.copy(),
        }

    axes[0].set_ylabel("Depth [m]")
    plt.show(block=True)
    return results


def compare_assigned_flux(cell):
    """Local exponential reference using LIMIT's final face coefficients.

    Conditioned on the numerical first-layer value.
    This is an internal-transport diagnostic, not Braud Eq. (43).
    """
    import numpy as np
    import matplotlib.pyplot as plt

    kwargs = iso.iso_delta.test_case_args(testcase=6, BA=True)
    layers = cell.layers
    n = len(layers)

    direction = np.sign(
        layers[-1].lower_boundary - layers[0].upper_boundary
    )
    surface = layers[0].upper_boundary

    z = np.array([
        direction * (
            0.5 * (l.upper_boundary + l.lower_boundary) - surface
        )
        for l in layers
    ])

    groups = [
        cell.connections_l_adv,
        cell.connections_v_adv,
        cell.connections_l_diff,
        cell.connections_v_diff,
    ]

    if any(len(group) != n - 1 for group in groups):
        raise ValueError("Expected N-1 connections in each internal group.")

    fig, axes = plt.subplots(
        1, 2, figsize=(8, 4.5),
        sharey=True, constrained_layout=True,
    )
    results = {}

    for ax, isotope in zip(axes, ["2H", "18O"]):
        c_standard = iso.iso_delta.delta_to_concentration(0.0, isotope)

        numerical_ratio = np.array([
            l.get_conc_iso_liquid(isotope) / c_standard
            for l in layers
        ])

        reference_ratio = np.empty(n)
        reference_ratio[0] = numerical_ratio[0]

        # Signed downward isotope flux divided by standard concentration.
        # Feed water is 0‰; steady bottom supply balances evaporation.
        J = -float(cell.q_evap)

        curve_z = []
        curve_ratio = []

        for j, (la, va, ld, vd) in enumerate(zip(*groups)):
            if la.climate != "arid" or va.climate != "arid":
                raise ValueError("This diagnostic assumes arid-mode connections.")

            if any(
                connection.left_node is not layers[j]
                or connection.right_node is not layers[j + 1]
                for connection in (la, va, ld, vd)
            ):
                raise ValueError("Connections must follow surface-to-bottom order.")

            A = (
                la.get_flux()
                + va.calc_flux_liquid(Isotopologue=isotope, **kwargs)
            )

            G = (
                ld.calc_flux(Isotopologue=isotope, **kwargs)
                + vd.calc_flux_liquid(Isotopologue=isotope, **kwargs)
            )

            if not np.isfinite(G) or G <= 0:
                raise ValueError("Positive finite diffusion conductance required.")

            # Solve J = A*r - K*dr/dz over this centre-to-centre segment.
            # With K = G*distance, the segment Peclet number is A/G.
            fraction = np.linspace(0.0, 1.0, 21)
            pe = A / G
            start = reference_ratio[j]

            if abs(pe) < 1e-8:
                segment = start - (J / G) * fraction
            else:
                segment = (
                    start
                    + np.expm1(pe * fraction) * (start - J / A)
                )

            reference_ratio[j + 1] = segment[-1]
            curve_z.extend(z[j] + fraction * (z[j + 1] - z[j]))
            curve_ratio.extend(segment)

        numerical_delta = 1000.0 * (numerical_ratio - 1.0)
        reference_delta = 1000.0 * (reference_ratio - 1.0)
        curve_delta = 1000.0 * (np.asarray(curve_ratio) - 1.0)

        ax.plot(
            curve_delta, -np.asarray(curve_z), "b-",
            label="Assigned-flux exponential diagnostic",
        )
        ax.plot(
            numerical_delta, -z, "r.", markersize=4,
            label="LIMIT",
        )
        ax.set_xlabel(
            r"$\delta^{2}$H (‰)" if isotope == "2H"
            else r"$\delta^{18}$O (‰)"
        )
        ax.set_ylim(-min(0.1, z[-1]), 0.0)
        ax.legend(fontsize=7, frameon=False)

        print(
            isotope,
            "maximum centre-value difference below first layer [‰]:",
            np.max(np.abs(
                numerical_delta[1:] - reference_delta[1:]
            )),
        )

        results[isotope] = {
            "depth": z.copy(),
            "numerical": numerical_delta,
            "reference_at_centres": reference_delta,
        }

    axes[0].set_ylabel("Depth [m]")
    plt.show(block=True)
    return results


