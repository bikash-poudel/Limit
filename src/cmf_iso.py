# loads libraries

import cmf
import numpy as np
from datetime import datetime, timedelta
from math import exp
import matplotlib.pyplot as plt

from Boundary_conditions import BoundaryCondition, iso_atmosphere
from iso_layers import iso_layer
from iso_time import iso_time
from Visualize import Visualize

import solve_iso_transport as iso


## Check for pv_sat,  0.002166 fron SLI Cvsat??? gas constant ???


def main():
    # Initialize

    # Likelihood Tests (Mathieu and Bariac (1996).)
    ignorealphai = True
    ignorealphaik = True
    ignoredl = True
    ignoredvi = True

    ignore = [ignorealphai, ignorealphaik, ignoredl, ignoredvi]

    # Initial isotope concentrations
    c0_2H_s = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
    c0_18O_s = iso.delta_to_concentration(delta_i=-8, solute_i="18O")
    c0_2H_atm = iso.delta_to_concentration(delta_i=-65, solute_i="2H")
    c0_18O_atm = iso.delta_to_concentration(delta_i=-8, solute_i="18O")

    # CMF setup
    P, C = iso_project()  # define project and cell

    # Top layer thickness of e10-5 m as per SISPAT_iso ???
    L_boundaries = np.arange(0.01, 1.01, 0.05)  # define layer thicknesses
    add_cmf_layers(C, L_boundaries, rtn_curve())  # define retention curve to all cell layers
    # TODO: is Richards connection to 1st layer default ???
    add_connections(C, connection=cmf.Richards)  # apply layer connection
    et_pot = -1.005e-5  # kg/(m2 * s) OR mm / s
    et_cmf = True  # True: obtain evaporation flux from cmf, False: Provide neuman evaporation flux

    # Simulation Period
    start = cmf.Time(1, 1, 2020)
    end = start + timedelta(days=200)  # run cmf for 300 days
    dt = cmf.h  # time step

    # Define solver
    solver = cmf.CVodeBanded(P, 1e-6)
    solver.t = start

    # setup for isotope (isotope storages)
    atm = iso_atm(c0_2H_atm, c0_18O_atm)  # define atmosphere with initial isotope concentrations (isotope storage)
    i_layers = layers_iso(C, c0_2H_s, c0_18O_s)  # isotope layers (isotope storage)
    my_time = iso_time(final_time=(end - start) / cmf.h, delta_time=1, time_units='hours')  # Time function for isotope computation to be replaced with timedelta
    X = [atm, i_layers, my_time, et_pot, et_cmf]

    # Run cmf with iso
    C_iso, fluxes, theta, hs, evap, pot = run(C, solver, end, dt, X, ignore)

    C2H = np.array(C_iso['2H']).T  # Extracting isotope concentration time series for '2H'
    C2H_delta = to_delta(C2H)  # convert concentration to delta notation

    ql_up = np.array(fluxes[0]).T
    qv_up = np.array(fluxes[2]).T
    th = np.array(theta).T

    # Visualize
    plot = Visualize(i_layers, my_time)

    plot.iso_profile(C2H_delta, solute_i='2H', print_time_steps=24 * 5)
    #plot.iso_breakthrough(C2H_delta, solute_i='2H', print_steps=5)
    plot.profiles([ql_up[:, -1], qv_up[:, -1]], label=['ql_up', 'qv_up'], title='profiles')
    # plot.profile(th[:, -1], label='theta a the end')
    plot.breakthrough(th[0], label='top layer theta')
    plot.breakthrough(evap, label='evaporation')
    # plot.breakthrough(pot, y_label='m', label=None, title='matrix potential')
    plot.breakthrough(hs, y_label='[-]', label='surface relative humidity')


def iso_project():
    # Create a project
    p = cmf.project()

    # Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
    c = p.NewCell(0, 0, 0, 1000)

    return p, c


def rtn_curve():
    # Set up a soil retention curve
    # With parameters as per SISPAT (Braud et al. 2005) Yolo Light clay (Philip, 1957)
    Braud_Ksat = 0.0106272  # in m/d
    Braud_phi = 0.35  # porosity
    Braud_alpha = 0.193  # Scale value of the water pressure(m)
    Braud_n = 2.22
    Braud_m = 0.099
    Braud_theta_r = 0.01
    Braud_eta = 9.14

    vgm = cmf.VanGenuchtenMualem(Ksat=Braud_Ksat,
                                 phi=Braud_phi,
                                 alpha=Braud_alpha,
                                 n=Braud_n,
                                 m=Braud_m,
                                 theta_r=Braud_theta_r)

    # Oversaturation tolerence upto 1% for matrix pot = +1
    vgm.w0 = 0.9995  # See: https://philippkraft.github.io/cmf/cmf_tut_retentioncurve.html Oversaturation

    return vgm


def add_cmf_layers(cell, l_boundaries, r_curve):
    for d in l_boundaries:
        l = cell.add_layer(d, r_curve)
        l.wetness = 1  ##### Check tehts exceeds porosity if  l.wetness = 1.0
        # c.saturated_depth = 0  # check for wetness = 1, saturated depth = 0,


def add_connections(cell, connection=cmf.Richards):
    cell.install_connection(connection)


def et_boundary(cell, E_pot=-1.005e-5):
    potET = abs(E_pot) * 86400  # mm/day
    ETpot = cmf.timeseries.from_scalar(potET)

    for layer in cell.layers:
        cmf.timeseriesETpot(layer, cell.evaporation, ETpot)  # connection to top only or to all layers ???
        # et_pot_connection = cmf.timeseriesETpot(layer, c.evaporation, ETpot)


def add_boundary(cell):
    # Create a Neumann Boundary condition connected to top layer to account for flows due to evaporation
    E_act = cmf.NeumannBoundary.create(cell.layers[0])
    E_act.flux = cmf.timeseries.from_scalar(0.0)

    return E_act


def get_fluxes(cell, i_layers, time):

    ql_up, ql_down = get_ql(cell, time)
    qv_up, qv_down = get_qv(i_layers)

    return [ql_up, ql_down, qv_up, qv_down]


def get_ql(cell, time):
    # get fluxes of the next time step
    ql_up = []  # flux in m3/d during the last time step
    ql_down = []

    f = cmf.sec / cmf.day / cell.area  # convert to m per second
    qqq = 0.0

    for lr in cell.layers:

        if lr.upper is None and lr.lower is not None:
            ql_up.append(0.0 * f + qqq)
            ql_down.append(lr.flux_to(lr.lower, time) * f + qqq)

        elif lr.upper is not None and lr.lower is not None:
            ql_up.append(lr.flux_to(lr.upper, time) * f + qqq)
            ql_down.append(lr.flux_to(lr.lower, time) * f + qqq)
        elif lr.upper is not None and lr.lower is None:
            ql_up.append(lr.flux_to(lr.upper, time) * f + qqq)
            ql_down.append(0.0 * f)

    return ql_up, ql_down


def get_qv(i_layers):

    qv_up = []
    qv_down = []

    for l in range(len(i_layers)):
        if l == 0:
            qv_up.append(0.0)
            qv_down.append(vapour_flux(i_layers[l+1], i_layers[l]))

        elif l != 0 and l != len(i_layers)-1:
            qv_up.append(vapour_flux(i_layers[l-1], i_layers[l]))
            qv_down.append(vapour_flux(i_layers[l], i_layers[l+1]))

        elif l == len(i_layers)-1 :
            qv_up.append(vapour_flux(i_layers[l-1], i_layers[l]))
            qv_down.append(0.0)

    return qv_up, qv_down


def vapour_flux(l_up, l_down):

    dz = l_down.center - l_up.center

    dv_air_up = iso.dv(T=l_up.T)
    dv_soil_up = iso.dv_soil_air(dv_air_up, l_up.theta, l_up.theta_sat, l_up.tortuosity)
    dv_air_down = iso.dv(T=l_down.T)
    dv_soil_down = iso.dv_soil_air(dv_air_down, l_down.theta, l_down.theta_sat, l_down.tortuosity)

    hr_up = iso_layer.soil_relative_humidity(l_up.psi, l_up.T)
    hr_down = iso_layer.soil_relative_humidity(l_down.psi, l_down.T)

    mean_dv = (dv_soil_up * l_up.thickness + dv_soil_down * l_down.thickness) / dz
    mean_hr = (hr_up * l_up.thickness + hr_down * l_down.thickness) / dz
    mean_pv_sat = (iso.pv_sat(l_up.T) * l_up.thickness + iso.pv_sat(l_down.T) * l_down.thickness) / dz

    q_vh = mean_dv * mean_pv_sat * (hr_up-hr_down) / dz
    q_vt = mean_dv * mean_hr * (iso.pv_sat(l_up.T) - iso.pv_sat(l_down.T)) / dz

    qv = q_vh + q_vt  # m/s

    return qv


def iso_atm(initial_c_2H_atmosphere, initial_c_18O_atmosphere, T_atmosphere=303.0, rH_atmosphere=0.20):
    atm_iso = iso_atmosphere(
        initial_c_atmosphere={"2H": initial_c_2H_atmosphere, "18O": initial_c_18O_atmosphere},
        # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!)
        # (currently supported "2H" and/or "18O")
        initial_T_atmosphere=T_atmosphere,  # temperature in the atmosphere in Kelvin
        initial_Rh_atmosphere=rH_atmosphere)  # relative humidity of the atmosphere (-)

    return atm_iso


def layers_iso(cell, initial_c_2H_soil, initial_c_18O_soil):
    # Define layer properties for isotopes as per the cell layers

    layers = cell.layers

    my_layers = []
    for lr in layers:
        l_iso = iso_layer(upper_boundary=lr.upper_boundary,
                          lower_boundary=lr.lower_boundary,
                          initial_c_solutes={"2H": initial_c_2H_soil, "18O": initial_c_18O_soil},
                          # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!)
                          # (currently supported "2H" and/or "18O")
                          initial_theta=lr.theta,  # volumetric water content m3/m3
                          theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                          theta_sat=lr.porosity,  # volumetric water content m3/m3 at saturation
                          initial_T=303,  # soil temperature in Kelvin
                          porosity=lr.porosity,  # porosity of the soil m3/m3
                          tortuosity=0.67,  # tortuosity of the soil m/m
                          psi=lr.matrix_potential
                          )  # soil maric potential [m]

        my_layers.append(l_iso)

    return my_layers


def iso_BC(ET_pot=-1.005e-5):
    # ET_pot=-1.005e-5 # kg/(m2 * s) or[ mm s-1]

    # Boundary conditions
    U_boundary_conc = iso.delta_to_concentration(delta_i=-70, solute_i='2H')
    L_boundary_conc = iso.delta_to_concentration(delta_i=-65, solute_i='2H')

    BC = BoundaryCondition()
    BC.upper_boundary('atmosphere', ET_pot)  # neuman[flux] / dirichlet[constant conc] / atmosphere[potential evaporation]
    BC.lower_boundary('neuman', 0)  # TODO: need to check lower boundary

    return BC


def soil_relative_humidity(h, T, R=461.5, g=9.81):
    """
    Calculates the relative humidity at given soil depth for given pressure head and temperature  certain
    soil depth according to Braud et al. 2005 Eq. 46

    Description
    == == == == == =

     @param h: pressure head [m]

     @param R: perfect gas mass constant for water vapour [J kg-1]

     @param T: Temperaturtre in kelvin [k]

     @param g: acceleration due to gravity [ms-2]
    """

    hu = exp(g * h / (R * T))

    return hu


def run(cell, solver, end, dt, X, ignore):

    atm, layers, time, pot_et, et_cmf = X
    ign_alpha, ign_alphak, ign_dli, ign_dvi = ignore

    if not et_cmf:
        e_act = add_boundary(cell)  # add neuman flux at boundary
        evap = [e_act.fluxes(solver.t)[0][0]]  # evaporation flux
    elif et_cmf:  # use CMF to calculate the ETact based on a fixed ETpot
        et_boundary(cell, pot_et)
        evap = [cell.evaporation(solver.t)]  # evaporation flux

    q1, q2, q3, q4 = get_fluxes(cell, layers, solver.t)
    ql_up, ql_down, qv_up, qv_down = [q1], [q2], [q3], [q4]

    theta = [cell.layers.theta]  # water content
    hs = [soil_relative_humidity(cell.layers[0].matrix_potential, 303)]  # surface relative humidity
    pot = [cell.layers.potential]  # soil matrix potential

    ciso = {'2H': [[l.c_solutes['2H'] for l in layers]],  # initial isotopic concentration
            '18O': [[l.c_solutes['18O'] for l in layers]]}

    # The run time loop. Iterates over the outer timestep of the model
    for t in solver.run(solver.t, end, dt):
        # run cmf for one time step
        # solver.reset() # - Resets the solver history. Only needed if you really change something in the system

        theta.append(cell.layers.theta)
        Q = get_fluxes(cell, layers, solver.t)
        fluxes = Q + [[theta[-1], theta[-1]]]
        hs_surface = iso_layer.soil_relative_humidity(layers[0].psi, layers[0].T)

        # Run iso_top simulation
        Ci_t = iso.run_1D_model(atm, layers, fluxes, iso_BC(pot_et), time, hs_surface, solutes=['2H'],
                                ignore_alpha_i=ign_alpha,
                                ignore_alpha_i_k=ign_alphak,
                                ignore_dl_i=ign_dli,
                                ignore_dv_i=ign_dvi)

        update_layers = iso.update_c_i(Ci_t['2H'], '2H', layers)
        layers = update_layers

        ciso['2H'].append(Ci_t['2H'])

        # calculate actual Evaporation for the next time step
        l0 = layers[0]
        if not et_cmf:
            ha = iso_atmosphere.soil_surface_Rh(atm.Rh_atmosphere, l0.T, atm.T_atmosphere)  #
            ev = pot_et / 1000 * cell.area * 86400 * (hs_surface - ha) / (1 - ha)   # m3 per day
            e_act.flux = cmf.timeseries.from_scalar(ev)

        elif et_cmf:  # use CMF to calculate the ETact based on a fixed ETpot
            ev = cell.evaporation(solver.t)

        qlup, qldown, qvup, qvdown = Q
        ql_up.append(qlup),  ql_down.append(qldown), qv_up.append(qvup), qv_down.append(qv_down)

        hs.append(hs_surface)
        evap.append(ev)
        pot.append(cell.layers.potential)

    return ciso, [ql_up, ql_down, qv_up, qv_down], theta, hs, evap, pot


def to_delta(c_iso):
    # Convert Result into delta notation
    my_vectorized_function = np.vectorize(iso.concentration_to_delta, excluded=['solute_i'])
    c_delta = my_vectorized_function(c_iso, solute_i='2H')
    # C18O_delta = my_vectorized_function(C18O, solute_i='18O')

    return c_delta


if __name__ == '__main__':
    main()
