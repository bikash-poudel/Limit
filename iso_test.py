import cmf

from datetime import timedelta
import numpy as np

from src import *


def test_case(testcase=1):
    ignore = {'ignoredvi': True, 'ignoredli': True, 'ignorealphai': True,
              'ignorealphaik': True}  # Testcases: Mathieu and Bariac (1996)

    if testcase == 1 or testcase == 2:
        return ignore
    elif testcase == 3 or testcase == 4:
        ignore['ignorealphai'] = False
        return ignore
    elif testcase == 5:
        ignore['ignorealphai'] = False
        ignore['ignoredli'] = False

        return ignore
    elif testcase in [6, 7, 8]:

        ignore['ignorealphai'] = False
        ignore['ignoredli'] = False
        ignore['ignoredvi'] = False
        ignore['ignorealphaik'] = False

        return ignore

    else:
        raise NotImplementedError


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


def cmf_setup():
    # CMF setup
    P, C = cmf_project()  # define project and cell

    L_boundaries = np.arange(0.01, 1.01, 0.05)  # define layer thicknesses
    add_cmf_layers(C, L_boundaries, rtn_curve())  # define retention curve to all cell layers
    add_connections(C, connection=cmf.Richards)  # apply layer connection

    return P


###### iso ######
def _layers(c_iso, C_cmf):

    # define iso layers as per cmf layer
    layers = C_cmf.layers

    init_c_iso_2H = iso_delta.delta_to_concentration(-65, '2H')  # 0.15367056287906838
    init_c_iso_18O = iso_delta.delta_to_concentration(-8, '18O')

    id = 0
    for lr in layers:
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=lr.upper_boundary,
                                                lower_boundary=lr.lower_boundary,
                                                conc_iso_liquid={"2H": init_c_iso_2H, "18O": init_c_iso_18O},
                                                theta=lr.theta,
                                                theta_0=0.0,
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

    c_iso_2H = iso_delta.delta_to_concentration(-65, '2H')  # sli.civa(1) / sli.cva(1)  0.15367056287906838
    c_iso_18O = iso_delta.delta_to_concentration(-8, '18O')

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
    # Update layers
    theta = [l.theta for l in c_cmf.layers]
    T_soil = [303.17] * len(c_cmf.layers)
    rH = [None] * len(c_cmf.layers)
    psi = [l.matrix_potential for l in c_cmf.layers]

    c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c_iso, c_cmf, time):
    # Update states and fluxes for
    f = cmf.sec / cmf.day / c_cmf.area  # convert to m per second

    # soil fluxes
    ql = []
    for lr in c_cmf.layers:
        if lr.upper is not None and lr.lower is None:
            ql.append(0.0 * f)
        else:
            ql.append(lr.flux_to(lr.lower, time) * f)

    qv = [0] * len(c_cmf.layers)

    c_iso.update_liquid_fluxes(liquid_fluxes=ql)
    c_iso.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

    # boundary storage
    c_iso.update_neuman_boundary(q_neuman=0.0, c_neuman={"2H": 1.0, "18O": 1.0})


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
    c.install_connections()  # install storage connections between the layers
    c.add_neuman_boundary(soil_layer=c.layers[0])

    return p, P


def run(p, P, **kwargs):

    C = P.cells[0]  # current cell cmf_project
    c = p.get_cells()[0]  # current cell of iso_project

    # Simulation Period
    start = cmf.Time(1, 1, 2020)
    end = start + timedelta(days=50)  # run cmf for 300 days
    dt = cmf.h  # time step

    # Define solver
    solver = cmf.CVodeBanded(P, 1e-6)
    solver.t = start

    solutes = ["2H", "18O"]
    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}
    while solver.t < end:

        print(solver.t)

        solver(dt)

        c_iso["2H"].append(c.conc_2H), c_iso["18O"].append(c.conc_18O)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        update_storages(c_iso=c, c_cmf=C), update_boundaries(c_iso=c, c_cmf=C, time=solver.t)
        for solute in solutes:
            delta_t = 3600  # dt.AsSeconds()
            dc = p.run(Isotopologue=solute, delta_time=delta_t, error_tol=1e-11, **kwargs)

            c_t = list(np.array(c.get_conc_layers(Isotopologue=solute)) + np.array(dc))
            c.update_c_layers(conc_iso=c_t, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


ignore = iso_delta.test_case_args(testcase=1)
p_iso, p_cmf = iso_setup()
delta = run(p_iso, p_cmf, **ignore)

