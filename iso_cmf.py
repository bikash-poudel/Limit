
import cmf

from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from src import *


def visualize(p, delta, dz, Isotopologue='2H'):
    depth = [-l.center.z for l in p.get_cells()[0].layers]

    fig, ax = plt.subplots(figsize=(7, 10))
    # plt.figure(figsize=(7, 10))
    for tests in delta.keys():
        ax.plot(delta[tests][Isotopologue][-1][:dz], depth[:dz], label='test_{}'.format(tests))

    ax.set_xlabel(r'$\delta$ {}'.format(Isotopologue))
    ax.set_ylabel('depth $[m]$')
    ax.set_title('Initial testcases after 250 days')
    ax.grid()
    ax.legend(loc='lower right')
    plt.show()

    return fig


def _theta(psi):
    alpha = 1 / 19.3
    n = 2.22
    m = 0.099

    theta_sat = 0.35
    theta_r = 0.01

    S = (1 + (alpha * 100 * psi) ** n) ** -m

    return S * (theta_sat - theta_r) + theta_r


###########################################################
############################ cmf ##########################
###########################################################

def cmf_project():
    # Create a project
    p = cmf.project()

    # Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
    c = p.NewCell(0, 0, 0, 1000)

    return p, c


def rtn_curve():
    
    # Set up a soil retention curve
    # With parameters as per SISPAT (Braud et al. 2005) Yolo Light clay (Philip, 1957)

    vgm = cmf.VGM_BC_RetentionCurve_Windhorst(Ksat=0.0106272,  # m/d
                                              phi=0.35,  # porosity
                                              alpha=1 / 19.3,
                                              # inverse of air entry potential (Scale value of the water pressure(m))
                                              n=2.22,
                                              m=0.099,
                                              theta_r=0.01,
                                              eta=9.14,
                                              )

    vgm.l = 0.67
    # vgm.fit_w0(w1=1.01, Psi_p=1.0)    # Oversaturation tolerence upto 1% for matrix pot = +1
    # print(vgm.w0)
    vgm.w0 = 0.99999

    return vgm


def add_cmf_layers(cell, l_boundaries, r_curve):

    for d in l_boundaries:
        cell.add_layer(d, r_curve)
    cell.saturated_depth = 0


def add_connections(cell, connection=cmf.Richards):
    cell.install_connection(connection)


def cmf_boundary(P):

    cell = P.cells[0]

    summer = cmf.Weather(Tmin=30, Tmax=30, rH=20, wind=2.0)
    cell.set_weather(summer)

    stress = cmf.ContentStress(theta_d=0.18,
                               theta_w=0.07)  # mpot = −153 for WP, −3.3 for FC, Standard plant wilting threshold, 	Gravity drainage ends

    cell.set_uptakestress(stress)
    cell.vegetation.RootDepth = cell.layers[0].lower_boundary
    # cell.vegetation.CanopyPARExtinction = 1.5
    # cell.vegetation.Height = 10
    # cell.vegetation.CanopyClosure = 0

    cmf.PenmanMonteithET(cell.layers[0], cell.evaporation)

    # cell.install_connection(cmf.ShuttleworthWallace)

    # ETpot = cmf.timeseries.from_scalar(2)
    # cmf.timeseriesETpot(cell.layers[0], cell.evaporation, ETpot)


def cmf_setup():

    # CMF setup

    P, C = cmf_project()  # define project and cell

    L_boundaries = np.array([0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0, 1.25, 1.5,
                             1.9, 2.3, 2.8, 3.5, 4.2, 5.0])

    # L_boundaries = np.cumsum([0.01] * 100)

    add_cmf_layers(C, L_boundaries, rtn_curve())  # define retention curve to all cell layers
    add_connections(C, connection=cmf.Richards)  # apply layer connection
    cmf_boundary(P)

    return P


###########################################################
########################## iso ############################
###########################################################

def _atm(testcase):

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                         conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                         T=303.17,  # K
                         Rh_atmosphere=0.2,
                         Pa_atmosphere=1,  # as in sli
                         wind_speed=2,  # m/s
                         hc=10,  # canopy height [m] (e.g. 40)
                         d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                         z0m=0.1 * 10,
                         # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                         LAI=1.1,  # leaf area index (e.g. 2.0)
                         extku=1.5)

    return atm


def _layers(c_iso, C_cmf, testcase):

    """ c_iso: iso cell, c_cmf: cmf cell """

    ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=testcase)
    ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=testcase)
    init_c_iso_2H = iso_delta.delta_to_concentration(ali_2H, '2H')
    init_c_iso_18O = iso_delta.delta_to_concentration(ali_18O, '18O')

    # define iso layers according to cmf-layers
    id = 0
    for lr in C_cmf.layers:

        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=lr.upper_boundary,
                                                lower_boundary=lr.lower_boundary,
                                                conc_iso_liquid={"2H": init_c_iso_2H, "18O": init_c_iso_18O},
                                                theta=min(lr.theta, lr.porosity),
                                                theta_0=0.01,
                                                theta_sat=lr.porosity,
                                                tortuosity=0.67,
                                                T=303.17,
                                                psi=lr.matrix_potential
                                                )

        id += 1
        c_iso.add_layer(new_layer)

    return c_iso


def update_storages(c_iso, c_cmf):

    """ c_iso: iso cell, c_cmf: cmf cell """

    # atmosphere need not be updated #

    # Update layers
    theta = [min(l.theta, l.porosity) for l in c_cmf.layers]
    T_soil = [303.17] * len(c_cmf.layers)
    rH = [None] * len(c_cmf.layers)
    psi = [min(l.matrix_potential, 0) for l in c_cmf.layers]

    c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c_iso, c_cmf, time):

    """ c_iso: iso cell, c_cmf: cmf cell """

    f = 1 / c_cmf.area / cmf.day.AsSeconds()  # factor [mm / day] to [m / s]

    ###################### Soil fluxes ########################
    ql = []
    for lr in c_cmf.layers:
        if lr.upper is not None and lr.lower is None:
            ql.append(0.0 * f)
        else:
            ql.append(lr.flux_to(lr.lower, time) * f)

    c_iso.update_vapor_fluxes(vapor_fluxes=None)  # vapor_fluxes = None: self compute vapor flux

    ####### disntegrate fluxes into liquid and vapor fluxes ##########
    qv = np.array(c_iso.vapor_fluxes + [0])
    mask = np.abs(qv) > np.abs(ql)
    ql_liquid = np.where(mask, 0, ql - qv)
    qv_update = np.where(mask, ql, qv)
    ##################################################################

    c_iso.update_vapor_fluxes(vapor_fluxes=qv_update)
    c_iso.update_liquid_fluxes(liquid_fluxes=ql_liquid)

    ##################### Surface fluxes ##############################
    atm = flux_atmosphere(atmosphere=c_iso.atmosphere, top_layer=c_iso.layers[0])
    fql, fqv = atm.E_liquid() / atm.E_total(), atm.E_vapor() / atm.E_total()  # factors determining surface fluxes

    q_evap = c_cmf.evaporation(time) * f
    ql_surface = - q_evap * fql
    qv_surface = - q_evap * fqv
    c_iso.update_evaporation(q_ev=q_evap, ql_surface=ql_surface, qv_surface=qv_surface, T_surface=303.17)


def iso_setup(testcase):
    ################## cmf ################
    P = cmf_setup()  # cmf project
    C = P.cells[0]  # current cell cmf

    ################# iso ##################
    atm = _atm(testcase)

    p = iso_project()
    p.new_cell(atmosphere=atm, area=C.area, x=C.x, y=C.y, z=C.z)

    c = p.get_cells()[0]  # get current iso cell
    _layers(c, C, testcase)  # add cmf layers to the current iso_cell

    ########### Install connections ###########
    c.install_connections(liquid_diffusion=True, vapor_diffusion=True,
                          liquid_advection=True, vapor_advection=True)  # install storage connections between the layers
    c.add_evaporation()

    return p, P


def run(P, p, sim_period=50, dt=1, solutes=["2H", "18O"], testcase=1):

    """
    P: cmf project
    p: iso project
    simulation period: days
    dt: hours
    """

    print('Testcase:', testcase)
    kwargs = iso_delta.test_case_args(testcase=testcase)

    C = P.cells[0]  # current cell cmf_project
    c = p.get_cells()[0]  # current cell of iso_project

    ###################### Solver ###################################################################
    solver = cmf.CVodeBanded(P)
    start = datetime(2024, 1, 1)
    end = start + timedelta(days=sim_period)
    timestep = timedelta(hours=dt)

    ################# output variables  #############################################################
    c_iso_delta = {'2H': [], '18O': []}
    m_pot, theta, ql, qv = [], [], [], []
    ev, ql_surface, qv_surface = [], [], []
    error = {'2H': [], '18O': []}
    #################################################################################################

    delta_t = timestep.seconds
    for t in solver.run(start, end, timestep):

        print(t)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        update_storages(c_iso=c, c_cmf=C), update_boundaries(c_iso=c, c_cmf=C, time=t)

        ############################################################################################
        ############################# run solver ####################################################
        for solute in solutes:

            dc, err = p.run(Isotopologue=solute, delta_time=delta_t, error_tolerance=1e-09, **kwargs)
            error[solute].append(err)

        ############################################################################################
        ############################################################################################
        theta.append([l.theta for l in c.layers])
        m_pot.append([l.matrix_potential for l in C.layers])
        ql.append(c.liquid_fluxes), qv.append(c.vapor_fluxes)

        ev.append(c.q_evap), ql_surface.append(c.ql_surface), qv_surface.append(c.qv_surface)
        ############################################################################################
        ############################################################################################

    return c_iso_delta, [theta, m_pot, ql, qv], [error, ev, ql_surface, qv_surface]


def run_testcases(test_cases=[1], solutes=["2H", "18O"], dt=12, sim_period=250):
    """"
    test_cases: Braud test cases
    sim_period: days
    dt: hours
    """

    delta = {}
    for Testcase in test_cases:
        delta[Testcase] = {}

        P_iso, P_cmf = iso_setup(testcase=Testcase)

        d, X, Y = run(P_cmf,
                      P_iso,
                      sim_period=sim_period,
                      dt=dt,
                      solutes=solutes,
                      testcase=Testcase)

        delta[Testcase]["2H"] = d["2H"]
        delta[Testcase]["18O"] = d["18O"]

    return P_iso, delta, X, Y


# p_iso, delta, X, Y = run_testcases(test_cases=[1, 2, 3, 4, 5, 6])


