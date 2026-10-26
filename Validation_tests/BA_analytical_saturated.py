
"""S

aturated isothermal reference: Braud et al. (2005a), Eqs. 40–41.

Pure molecular kinetic fractionation is Braud's BA83 reference (nK=1).
If the numerical boundary uses another kinetic factor, pass the SAME factor
explicitly through alpha_k; the result then generalizes Eq.41 rather than
reproducing Table 3's pure-molecular setting.

"""

import warnings
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from limit import iso_fluxes as isof


def surface_delta(alpha_v_over_l, alpha_k, humidity, delta_feed, delta_vapor):

    if not 0 <= humidity < 1:
        raise ValueError('Atmospheric relative humidity must be a fraction in [0,1).')
    if alpha_v_over_l <= 0 or alpha_k <= 0:
        raise ValueError('Fractionation factors must be positive.')
    return 1000 * (((1-humidity)*alpha_k*(1+delta_feed/1000)
                    + humidity*(1+delta_vapor/1000))/alpha_v_over_l - 1)


def profile(z, evaporation_ms, diffusivity_m2s, delta_feed, delta_surface):

    z = np.asarray(z, dtype=float)
    if np.any(z < 0) or evaporation_ms <= 0 or diffusivity_m2s <= 0:
        raise ValueError('Need nonnegative depth, positive evaporation and diffusivity.')
    return delta_feed + (delta_surface-delta_feed)*np.exp(
        -evaporation_ms*z/diffusivity_m2s)


def solve_saturated(cell, *, delta_feed=None, delta_vapor=None,
                      alpha_k=None, output_dir='analytical_saturated',
                      saturation_tolerance=0.0055, temperature_tolerance=0.05,
                      show=True):

    """Extract final LIMIT fields, plot reference and save comparison values.

    delta_feed / delta_vapor: isotope dictionaries in per mil. Defaults match
    the supplied benchmark: feed 0/0, vapor -100/-14. Pass explicitly when
    forcing differs. alpha_k: optional isotope dict of kinetic factors.

    """

    delta_feed = {'2H': 0., '18O': 0.} if delta_feed is None else delta_feed
    delta_vapor = {'2H': -100., '18O': -14.} if delta_vapor is None else delta_vapor
    layers = cell.layers

    upper_raw = np.array(
        [l.upper_boundary for l in layers], dtype=float
    )
    lower_raw = np.array(
        [l.lower_boundary for l in layers], dtype=float
    )

    # Convert either coordinate orientation to depth below the surface.
    surface = upper_raw[0]
    direction = np.sign(lower_raw[-1] - surface)

    if direction == 0:
        raise ValueError("The soil column has zero total depth.")

    upper = direction * (upper_raw - surface)
    lower = direction * (lower_raw - surface)

    thickness = lower - upper
    z = 0.5 * (upper + lower)

    if np.any(upper < -1e-12) or np.any(thickness <= 0):
        raise ValueError(
            "Layers must be ordered from the soil surface downward "
            "with positive thickness."
        )

    T = np.array([l.T for l in layers])
    theta = np.array([l.theta for l in layers])
    theta_sat = np.array([l.theta_sat for l in layers])
    if not np.all(np.isfinite(np.r_[T, theta, theta_sat])):
        raise ValueError('Non-finite final soil state.')

    max_deficit = np.max(np.abs(theta-theta_sat))
    if max_deficit > saturation_tolerance:
        raise ValueError(f'Saturated reference not applicable: maximum water-content '
                         f'departure = {max_deficit:.6g} m3/m3. Run the saturated configuration.')

    if np.ptp(T) > temperature_tolerance or abs(cell.atmosphere.T-T.mean()) > temperature_tolerance:
        raise ValueError('Isothermal reference requires uniform soil and matching air '
                         'temperature. Disable heat evolution for this benchmark.')

    area = float(cell.area)
    if area <= 0:
        raise ValueError('Cell area must be positive.')
    E = float(cell.q_evap)/area

    h = float(cell.atmosphere.Rh)

    liquid = isof.liquid_diffusion_base_class()
    vapor = isof.vapor_diffusion_base_class()
    dense_z = np.linspace(0, lower[-1], 2001)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5), sharey=True,
                             constrained_layout=True)
    results = {}
    for ax, isotope in zip(axes, ['2H', '18O']):
        D = np.array([liquid.dl_i(T=l.T, Isotopologue=isotope,
                                  theta=l.theta, tortuosity=l.tortuosity)
                      for l in layers])

        # Constant effective diffusivity for the saturated reference.
        Dsat = np.array([
            liquid.dl_i(
                T=float(T.mean()),
                Isotopologue=isotope,
                theta=l.theta_sat,
                tortuosity=l.tortuosity,
            )
            for l in layers
        ])

        if not np.all(np.isfinite(Dsat)) or np.any(Dsat <= 0):
            raise ValueError("Invalid saturated liquid diffusivity.")

        if not np.allclose(Dsat, Dsat[0], rtol=1e-6, atol=0):
            raise ValueError(
                "The analytical solution requires homogeneous diffusivity."
            )

        Di = float(Dsat[0])

        alpha = float(layers[0].alpha_i(
            Isotopologue=isotope,
            T=float(T.mean()),
        ))

        if alpha_k is None:
            dv = vapor.dv_free_air(T=float(T.mean()))
            dvi = vapor.dv_i(dv=dv, Isotopologue=isotope)
            kinetic = float(dv / dvi)
        else:
            kinetic = float(alpha_k[isotope])

        ds = surface_delta(
            alpha,
            kinetic,
            h,
            delta_feed[isotope],
            delta_vapor[isotope],
        )

        curve = profile(
            dense_z, E, Di, delta_feed[isotope], ds
        )

        # Exact analytical average over each numerical layer.
        k = E / Di
        average = (
                delta_feed[isotope]
                + (ds - delta_feed[isotope])
                * np.exp(-k * upper)
                * (-np.expm1(-k * thickness))
                / (k * thickness)
        )

        numerical = np.asarray(
            getattr(cell, f"conc_{isotope}_delta"),
            dtype=float,
        ).copy()

        if numerical.shape != z.shape:
            raise ValueError("Numerical profile and depth sizes differ.")

        if not np.all(np.isfinite(numerical)):
            raise ValueError("Numerical isotope profile contains invalid values.")

        residual = numerical - average
        rmse = float(np.sqrt(np.mean(residual ** 2)))

        print(
            f"{isotope}: "
            f"surface analytical delta={ds:.6g} per mil; "
            f"decay={k:.6g} /m; "
            f"RMSE={rmse:.6g} per mil"
        )

        # Analytical solution: blue dashed line.
        ax.plot(
            curve,
            -dense_z,
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

        # Panel labels match the unsaturated figure.
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

        np.savetxt(
            output / f"saturated_{isotope}.csv",
            np.column_stack([z, numerical, average, residual]),
            delimiter=",",
            header=(
                "depth_m,numerical_permil,"
                "analytical_layer_mean_permil,residual_permil"
            ),
            comments="",
        )

        results[isotope] = {
            "delta_surface": ds,
            "decay_per_m": k,
            "rmse_permil": rmse,
        }

    # Outside the isotope loop, but still inside compare_saturated().
    # axes[0].set_ylabel("Depth [m]")
    axes[0].set_ylabel("Depth [m]", fontsize=10)

    fig.savefig(
        output / "saturated_comparison.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        output / "saturated_comparison.png",
        dpi=600,
        bbox_inches="tight",
    )

    if show:
        plt.show()
    else:

        plt.close(fig)

    return results


    