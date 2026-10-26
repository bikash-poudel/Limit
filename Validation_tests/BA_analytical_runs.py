
import sys

sys.path.insert(0,r"D:\Isotope transport\Scripts\Limit\Validation_tests")

import importlib
import numpy as np
import Barnes_Alison_saturated_isothermal as Ba1
import Barnes_Alison_unsaturated_nonisothermal as Ba2
importlib.reload(Ba1)
importlib.reload(Ba2)

import BA_analytical_saturated as saturated_analytical
import BA_analytical_unsaturated as unsaturated_analytical
importlib.reload(saturated_analytical)
importlib.reload(unsaturated_analytical)

import matplotlib.pyplot as plt

BA = "unsaturated"
dt, sim = 1, 250
if BA == "saturated":
    print("Starting saturated case")

    p_iso, delta, X, Y = Ba1.run_testcases(solutes=['2H', '18O'],
                                          bottom_potential= 1.3, #1.30,
                                          dt=dt,  # hours
                                          sim_period=sim,  # days
                                          BA=True
                                          )

    importlib.reload(saturated_analytical)
    print("Starting saturated analytical comparison")

    comparison = saturated_analytical.solve_saturated(
        p_iso.get_cells()[0],
        delta_feed={"2H": 0.0, "18O": 0.0},
        delta_vapor={"2H": -100.0, "18O": -14.0},
    )

    print("Saturated case finished")

else:
    print("Starting unsaturated case")
    # ----------------------------------------------------
    p_iso, delta, X, Y = Ba2.run_testcases(
        solutes=["2H", "18O"],
        bottom_potential=-0.4,
        dt=dt,          # Hours:
        sim_period=sim,  # Days:
        BA=True,
    )

    c = p_iso.get_cells()[0]
    temperatures = np.asarray(X[4], dtype=float)

    print("Temperatures finite:", np.isfinite(temperatures).all())
    print("Temperature range [K]:",
          temperatures.min(), temperatures.max())
    print("Final surface temperature [K]:", c.connection_evap.T_surface)
    print("Final theta range:",
          min(l.theta for l in c.layers),
          max(l.theta for l in c.layers))

    for isotope in ["2H", "18O"]:
        values = np.asarray(delta[6][isotope], dtype=float)
        print(isotope, "finite:", np.isfinite(values).all(),
              "range [‰]:", (values.min(), values.max()))

    comparison = unsaturated_analytical.solve_unsaturated(p_iso.get_cells()[0],
                                                         n_reference = 100)
    diff = comparison['2H']['analytical'] - comparison['2H']['numerical']
    print('max diff:', min(diff))


