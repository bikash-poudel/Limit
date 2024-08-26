
import cmf

from datetime import datetime, timedelta
import numpy as np

import matplotlib.pyplot as plt

from src import *


def visualize(p, delta, Isotopologue='2H'):

    depth = [-l.center.z for l in p.get_cells()[0].layers]

    for tests in delta.keys():
        plt.plot(delta[tests][Isotopologue], depth, label='test_{}'.format(tests))

    plt.xlabel('delta_{}'.format(Isotopologue))
    plt.ylabel('depth [m]')
    plt.title('initial testcases after 250 days')
    plt.legend()
    plt.gca().set_aspect(aspect=100)
    plt.show()


###### cmf ######
def cmf_project():
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


def cmf_boundary(P):

    cell = P.cells[0]

    #evap = P.NewNeumannBoundary('evap', cell.layers[0])
    #evap.flux = - 0.866592  # cum day-1

    ETpot = cmf.timeseries.from_scalar(0.2)
    cmf.timeseriesETpot(cell.layers[0], cell.evaporation, ETpot)

    # Create the boundary condition
    #gw = P.NewOutlet('groundwater', x=0, y=0, z=-1.01)
    q_hot = P.NewNeumannBoundary('q_bot', cell.layers[-1])
    q_hot.flux = 0.008

    # Set the potential
    # gw.potential = 1.01
    # Connect the lowest layer to the groundwater using Richards percolation
    # cmf.Richards(cell.layers[-1], gw)


    # summer = cmf.Weather(Tmin=16, Tmax=28, rH=50, wind=1.0, sunshine=0.9, daylength=16, Rs=26)
    # cell.set_weather(summer)
    # cmf.TurcET(cell.layers[0], cell.evaporation)


def cmf_setup():
    # CMF setup
    P, C = cmf_project()  # define project and cell

    L_boundaries = np.arange(0.01, 1.01, 0.05)  # define layer thicknesses
    add_cmf_layers(C, L_boundaries, rtn_curve())  # define retention curve to all cell layers
    add_connections(C, connection=cmf.Richards)  # apply layer connection

    cmf_boundary(P)

    return P


###### iso ######
def _layers(c_iso, C_cmf):

    # define iso layers as per cmf layer
    layers = C_cmf.layers

    init_c_iso_2H = iso_delta.delta_to_concentration(0, '2H')  # 0.15367056287906838
    init_c_iso_18O = iso_delta.delta_to_concentration(0, '18O')

    id = 0
    for lr in layers:
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=lr.upper_boundary,
                                                lower_boundary=lr.lower_boundary,
                                                conc_iso_liquid={"2H": init_c_iso_2H, "18O": init_c_iso_18O},
                                                theta=lr.theta,
                                                theta_0=0.01,
                                                theta_sat=lr.porosity,
                                                tortuosity=0.67,
                                                T=303.17,
                                                psi=lr.matrix_potential
                                                )
        id += 1
        c_iso.add_layer(new_layer)
    return c_iso


def _atm():
    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = 303.17  # sli.Tatm(dt) + Tzero_sli
    Rh_atm = 0.2  # sli.R_humidity_atm(dt)
    wind_speed = 0.2  # sli.wind_speed(dt)

    c_iso_2H = iso_delta.delta_to_concentration(67, '2H')  # sli.civa(1) / sli.cva(1)  0.15367056287906838
    c_iso_18O = iso_delta.delta_to_concentration(31.9, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=Tatm, Rh_atmosphere=Rh_atm,
                         Pa_atmosphere=patm,
                         wind_speed=wind_speed,
                         hc=10,  # canopy height [m] (e.g. 40)
                         d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.1 * 10,
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=1.1,  # leaf area index (e.g. 2.0)
                         extku=1.5)
    return atm


def update_storages(c_iso, c_cmf):

    # atmosphere need not be updated

    # Update layers
    theta = [l.theta for l in c_cmf.layers]
    T_soil = [303.17] * len(c_cmf.layers)
    rH = [None] * len(c_cmf.layers)
    psi = [l.matrix_potential for l in c_cmf.layers]

    c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c_iso, c_cmf, time):

    """c_iso: iso cell, c_cmf: cmf cell"""
    # Update states and fluxes for
    f = cmf.sec / cmf.day / c_cmf.area  # convert to m per second

    # soil fluxes
    ql = []
    qv = [0] * len(c_iso.layers)
    for lr in c_cmf.layers:
        if lr.upper is not None and lr.lower is None:
            ql.append(0.0 * f)
        else:
            #ql.append(lr.flux_to(lr.lower, time) * f)
            ql.append(lr.flux_to(lr.lower, time) * f)

    c_iso.update_liquid_fluxes(liquid_fluxes=ql)
    #c_iso.update_vapor_fluxes(vapor_fluxes=qv)

    # boundary storage
    c_iso_2H = iso_delta.delta_to_concentration(0, '2H')  # sli.civa(1) / sli.cva(1)  0.15367056287906838
    c_iso_18O = iso_delta.delta_to_concentration(0, '18O')

    f = 1 / c_cmf.area / 86400  # m3 day-1 to ms-1
    q_evap = c_cmf.layers[0].fluxes(time)[-1][0] * f  # c_cmf.evporation.fluxes(time)[0][0] * f  # c_cmf.layers[0].fluxes(time)[-1][0] * f

    c_iso.update_evaporation(q_ev=q_evap, T_surface=303.17, ql_surface=0.0)
    # c_iso.update_dirichlet_boundary(c_dirichlet={"2H": c_iso_2H, "18O": c_iso_18O})


def iso_setup():

    P = cmf_setup()  # cmf project
    C = P.cells[0]  # current cell cmf

    # setup for isotope (isotope storages)
    atm = _atm()
    # Iso project
    p = iso_project()
    p.new_cell(atmosphere=atm, area=C.area, x=C.x, y=C.y, z=C.z)

    c = p.get_cells()[0]  # get current iso cell
    _layers(c, C)  # add cmf layers to the current iso_cell

    #####Install connections#######
    c.install_connections(vapor_advection=False, vapor_diffusion=False)  # install storage connections between the layers

    c.add_evaporation()
    # c.add_dirichlet_boundary(soil_layer=c.layers[-1])

    return p, P


def run(p, P, sim_period=50, dt=1, **kwargs):

    """"
    sim_period: days
    dt: hours
    """

    C = P.cells[0]  # current cell cmf_project
    c = p.get_cells()[0]  # current cell of iso_project

    # Define solver
    solver = cmf.CVodeBanded(P)
    start = datetime(2024, 1, 1)
    end = start + timedelta(days=sim_period)
    timestep = timedelta(hours=dt)

    solutes = ["2H", "18O"]
    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}
    for t in solver.run(start, end, timestep):

        print(t)
        c_iso["2H"].append(c.conc_2H), c_iso["18O"].append(c.conc_18O)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        print([l.theta for l in C.layers])

        update_storages(c_iso=c, c_cmf=C), update_boundaries(c_iso=c, c_cmf=C, time=t)
        for solute in solutes:

            delta_t = dt * 3600  # dt.AsSeconds()
            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=None, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


def run_testcases(test_cases):

    delta = {}
    for Testcase in test_cases:
        print('Testcase:', Testcase)
        cases = iso_delta.test_case_args(testcase=Testcase)

        delta[Testcase] = {}

        p_iso, p_cmf = iso_setup()
        d = run(p_iso, p_cmf, sim_period=50, dt=1, **cases)

        # delta at the end of simulation for each test cases
        delta[Testcase]["2H"] = d["2H"][-1]
        delta[Testcase]["18O"] = d["18O"][-1]

    return p_iso, delta


p_iso, delta = run_testcases(test_cases=[4, 5, 6])
visualize(p_iso, delta, Isotopologue='2H')
visualize(p_iso, delta, Isotopologue='18O')


