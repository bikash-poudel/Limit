
import cmf
import iso_cmf as iso

from datetime import datetime, timedelta
import numpy as np
from math import exp

###########################################################
############################ CMF ##########################
###########################################################


def cmf_boundary(P):

    cell = P.cells[0]

    #################### Atmospheric Boundary ###################################
    summer = cmf.Weather(Tmin=30, Tmax=30, rH=20, wind=2.0)

    stress = cmf.ContentStress(theta_d=0.08, # 0.18,
                               theta_w=0.01) # 0.07)  # mpot = −153 for WP, −3.3 for FC, Standard plant wilting threshold, 	Gravity drainage ends

    cell.set_weather(summer)
    cell.set_uptakestress(stress)
    cell.vegetation.RootDepth = cell.layers[0].lower_boundary

    # cmf.PenmanMonteithET(cell.layers[0], cell.evaporation)

    ETpot = cmf.timeseries.from_scalar(7)
    cmf.timeseriesETpot(cell.layers[0], cell.evaporation, ETpot)

    ##################### lower boundary ########################################
    # Create a  outlet (Dirichlet) boundary
    # TODO: Is surface water created (surface water flux ultimately is non zero)
    outlet = P.NewOutlet('Lower boundary')
    outlet.potential = - 0.7  # 1.31 (18O), 1.33 (2H) for steadystate, isothermal, -0.6(2H), -0.8(18O)  for non isothermal
    outlet.is_source = True
    cmf.Richards(cell.layers[-1], outlet)

    """

    ####################### free drain to GW #####################################
    # Create the boundary condition
    gw = P.NewOutlet('groundwater', x=0, y=0, z=-6)
    # # Connect the lowest layer to the groundwater using Richards percolation
    gw.potential = 0
    cmf.Richards(cell.layers[-1], gw)
    # cmf.FreeDrainagePercolation(cell.layers[-1], gw) # for non steady state
    
    """


def cmf_setup():

    # CMF setup
    P, C = iso.cmf_project()  # define project and cell

    # L_boundaries = np.array([0.008, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 1.0,
      #                       1.25, 1.5, 1.9, 2.3, 2.8, 3.5, 4.2, 5.0])

    # L_boundaries = np.array([0.008, 0.016, 0.025, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1,
      #                        0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.6, 0.8, 0.95, 1.0])

    L_boundaries = np.cumsum([0.01] * 100)

    iso.add_cmf_layers(C, L_boundaries, iso.rtn_curve())  # define retention curve to all cell layers
    iso.add_connections(C, connection=cmf.Richards)  # apply layer connection
    cmf_boundary(P)

    return P


###########################################################
############################ ISO ##########################
###########################################################


def _atm():

    atm_2H, atm_18O = -100, -14
    c_iso_2H = iso.iso_delta.delta_to_concentration(atm_2H, '2H')
    c_iso_18O = iso.iso_delta.delta_to_concentration(atm_18O, '18O')

    atm = iso.iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                             conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                             T=313.17,  # K
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


def aquifer():
    atm_2H, atm_18O = 0, 0
    c_2H = iso.iso_delta.delta_to_concentration(atm_2H, '2H')
    c_18O = iso.iso_delta.delta_to_concentration(atm_18O, '18O')

    aq = iso.iso_aquifer(conc_iso_liquid={"2H": c_2H, "18O": c_18O},
                         T=303.17
                         )

    return aq


def _layers(c_iso, C_cmf):

    """ c_iso: iso cell, c_cmf: cmf cell """

    ali_2H, ali_18O = 0, 0
    init_c_iso_2H = iso.iso_delta.delta_to_concentration(ali_2H, '2H')
    init_c_iso_18O = iso.iso_delta.delta_to_concentration(ali_18O, '18O')

    # define iso layers according to cmf-layers
    id = 0
    for lr in C_cmf.layers:

        new_layer = iso.iso_storages.iso_soil_layer(ID=id,
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
        """
        # For BA unsaturate non-isothermal case
        for l in c_iso.layers:
            z = -l.center.z
            l.T = 20 * (1 + exp(-20 * z)) + 273.15
        """

    return c_iso


def update_heat(c, dt):

    layers = c.layers
    Tatm = c.atmosphere.T
    qevap = c.q_evap / c.area

    soil = {'rho_s': 2650, 'c_s': 840, 'rho_w': 1000, 'c_w': 4186, 'lambda_dry': 0.3, 'lambda_sat': 2.0}
    hc = 10

    thickness = np.array([l.thickness for l in layers])
    T = np.array([l.T for l in layers])
    theta = np.array([l.theta for l in layers])
    rH = np.array([l.relative_humidity(l.psi, l.T) for l in layers])

    qliq = np.insert(c.liquid_fluxes, 0, c.q_evap) / c.area
    qvap = np.insert(c.vapor_fluxes, 0, 0) / c.area

    T, Ts = iso.compute_temperature(T=T,
                                    theta=theta,
                                    RH=rH,
                                    q_liq=qliq,
                                    q_vap=qvap,
                                    q_evap=qevap,
                                    soil=soil,
                                    dz=thickness,
                                    dt=dt,
                                    T_atm=Tatm,
                                    h_c=hc
                                    )

    return T, Ts


def update_storages(c_iso, c_cmf, dt):

    """ c_iso: iso cell, c_cmf: cmf cell """

    # atmosphere need not be updated #

    T_soil, Ts = update_heat(c=c_iso, dt=dt)

    # Update layers
    theta = [min(l.theta, l.porosity) for l in c_cmf.layers]
    # T_soil = [303.17] * len(c_cmf.layers)
    rH = [None] * len(c_cmf.layers)
    psi = [min(l.matrix_potential, 0) for l in c_cmf.layers]

    c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)


def update_boundaries(c_iso, c_cmf, time, dt):
    """ c_iso: iso cell, c_cmf: cmf cell """

    f = 1 / cmf.day.AsSeconds()  # factor [mm / day] to [m / s]

    ###################### Soil fluxes #####################
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
    atm = iso.flux_atmosphere(atmosphere=c_iso.atmosphere, top_layer=c_iso.layers[0])
    fql, fqv = atm.E_liquid() / atm.E_total(), atm.E_vapor() / atm.E_total()  # factors determining surface fluxes

    T_soil, Ts = update_heat(c=c_iso, dt=dt)

    q_evap = c_cmf.evaporation(time) * f
    ql_surface = - q_evap  # * fql
    qv_surface =  0 # - q_evap * fqv
    c_iso.update_evaporation(q_ev=q_evap, ql_surface=ql_surface, qv_surface=qv_surface, T_surface=Ts)

    ##################### Lower boundary fluxes ##############################
    l = c_cmf.layers[-1]  # bottom layer
    q = l.flux_to(l.connected_nodes[-1], time) * f  # flux from the dirichlet boundary

    c_iso_2H = iso.iso_delta.delta_to_concentration(0, '2H')
    c_iso_18O = iso.iso_delta.delta_to_concentration(0, '18O')
    c_iso.update_aquifer(q_aq=q, c_iso={"2H": c_iso_2H, "18O": c_iso_18O})


def iso_setup():

    ################## cmf ################
    P = cmf_setup()  # cmf project
    C = P.cells[0]  # current cell cmf

    ################# iso ##################
    atm = _atm()

    p = iso.iso_project()
    p.new_cell(atmosphere=atm, area=C.area, x=C.x, y=C.y, z=C.z)

    c = p.get_cells()[0]  # get current iso cell
    _layers(c, C)  # add cmf layers to the current iso_cell

    ########### Install connections ###########
    c.install_connections(liquid_diffusion=True, vapor_diffusion=True,
                          liquid_advection=True, vapor_advection=True)  # install storage connections between the layers
    c.add_evaporation()
    c.add_aquifer(aquifer=aquifer(), soil_layer=c.layers[-1])

    p.build_mapping()

    return p, P


def run(P, p, sim_period=50, dt=1, solutes=["2H", "18O"], testcase=1, BA=True):

    """
    P: cmf project
    p: iso project
    simulation period: days
    dt: hours
    """

    print('Testcase:', testcase)
    kwargs = iso.iso_delta.test_case_args(testcase=testcase, BA=BA)

    C = P.cells[0]  # current cell cmf_project
    c = p.get_cells()[0]  # current cell of iso_project

    ###################### Solver #######################################################
    solver = cmf.CVodeBanded(P)
    start = datetime(2024, 1, 1)
    end = start + timedelta(days=sim_period)
    timestep = timedelta(hours=dt)

    ################# output variables  #################################################
    c_iso_delta = {'2H': [], '18O': []}
    m_pot, theta, ql, qv, T = [], [], [], [], []
    ev, ql_surface, qv_surface = [], [], []
    error = {'2H': [], '18O': []}
    ######################################################################################

    delta_t = timestep.seconds
    for t in solver.run(start, end, timestep):

        print(t)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        update_boundaries(c_iso=c, c_cmf=C, time=t, dt=dt)
        update_storages(c_iso=c, c_cmf=C, dt=dt)

        #############################################################################################
        ############################# run solver ####################################################
        for solute in solutes:
            dc, err = p.run(Isotopologue=solute, delta_time=delta_t, error_tolerance=1e-12, **kwargs)
            error[solute].append(err)

        ############################################################################################
        ############################################################################################
        theta.append([l.theta for l in c.layers])
        m_pot.append([l.matrix_potential for l in C.layers])
        ql.append(c.liquid_fluxes), qv.append(c.vapor_fluxes)
        T.append([l.T for l in c.layers])

        ev.append(c.q_evap), ql_surface.append(c.ql_surface), qv_surface.append(c.qv_surface)
        ############################################################################################
        ############################################################################################

    return c_iso_delta, [theta, m_pot, ql, qv, T], [error, ev, ql_surface, qv_surface]


def run_testcases(solutes=["2H", "18O"], dt=12, sim_period=250, BA=True):
    """"
    test_cases: Braud test cases - [6] for BA test
    sim_period: days
    dt: hours
    """

    delta = {6: {}}
    P_iso, P_cmf = iso_setup()

    d, X, Y = run(P_cmf,
                  P_iso,
                  sim_period=sim_period,
                  dt=dt,
                  solutes=solutes,
                  testcase=6,
                  BA=BA
                  )

    delta[6]["2H"] = d["2H"]
    delta[6]["18O"] = d["18O"]

    return P_iso, delta, X, Y

# delta, X, Y = run_testcases()
