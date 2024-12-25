
import numpy as np
from math import sin, pi

from project import project
import storage


def _atm():
    atm = storage.atmosphere(T=283.15,  # temperature in the atmosphere in Kelvin
                             Rh_atmosphere=0.20,  # relative humidity of the atmosphere (-)
                             Pa_atmosphere=10 ** 5,  # atmospheric pressure (Pa)
                             initial_wind_speed=2.0,  # wind speed at top of canopy in m/s
                             R_net=0.0,  # net radiation absorbed
                             wind_speed=2.0,  # wind speed at top of canopy in m/s
                             hc=10,  # canopy height [m] (e.g. 40)
                             d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                             z0m=0.1 * 10,  # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                             LAI=1.1,  # leaf area index (e.g. 2.0)
                             extku=1.5)  # extinction coeff't for windspeed (e.g. 0.5))

    return atm


def _layers(c):

    lower_boundaries = np.arange(2.5, 54, 2.5)

    theta_init = 0.2
    tmp_ini = 25

    id = 0
    upper_boundary = 0
    for lower_boundary in lower_boundaries:
        new_layer = storage.soil_layer(upper_boundary,
                                       lower_boundary,
                                       theta=0.2,
                                       theta_sat=0.547,
                                       tortuosity=0.67,
                                       rH=0,
                                       head=1.0,
                                       T=25,
                                       T0=25,
                                       K_s=3.8e-3,  # cm /s
                                       phi_e=-13.0,  # cm
                                       b=6.53,
                                       clay=0.249,
                                       sand=0.022,
                                       silt=0.729,
                                       som=0.044,
                                       rho_b=1.2,
                                       rho_l=1.0
                                       )

        new_layer.head = new_layer.calc_head(theta_init)
        new_layer.T = tmp_ini

        id += 1
        upper_boundary = lower_boundary
        c.add_layer(new_layer)


def run(p, initial_dt=1, simulation_time=86400*250):

    c = p.cells[0]

    T = 0
    count = 0
    while T < simulation_time:
        count += 1
        dt = initial_dt

        print('count :', count)
        print('dt :', dt, T / 86400, ': days')
        print(c.layers[0].theta, c.layers[-1].theta)

        q_ev = 1e-7 + 2e-8 * sin(2 * pi * (T + dt) / 86400)  # cm / s
        c.q_evap = q_ev

        initial_dt = p.run(dt=dt, dt_max=60, epsilon=1e-4)

        T += initial_dt


def soil_setup():

    atm = _atm()

    p = project()
    p.new_cell(atm)

    c = p.cells[0]
    _layers(c)

    c.install_connections()
    c.add_evaporation(top_layer=c.layers[-1])

    c.hleft, c.hright, c.Tleft, c.Tright = -9.274162934176326e+03, 8.570390159801457e+05, 25, 30

    simulation_time = 86400 * 10
    initial_dt = 1

    run(p, initial_dt, simulation_time)

    return p


prj = soil_setup()




