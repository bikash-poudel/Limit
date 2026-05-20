import cmf
import numpy as np

import vapor_model as vp
from src import *

from datetime import datetime, timedelta

# Plot results
from pylab import *


###### cmf ######
class cmf_iso(object):

    def __init__(self, r_curve='VGM'):

        self.retention = r_curve

    @property
    def project(self):
        # Create a project
        p = cmf.project()
        # Create a cell at position (0,0,0) with 1000m2 size (making conversion from m3 to mm trivial)
        p.NewCell(0, 0, 0, 1000)

        return p

    @property
    def rtn_curve(self):
        # Set up a soil retention curve
        # With parameters as per SISPAT (Braud et al. 2005) Yolo Light clay (Philip, 1957)
        Braud_Ksat = 0.0106272  # in m/d
        Braud_phi = 0.35  # porosity
        Braud_alpha = 1 / 19.3  # is inverse of the air entry potential in :math:`cm^{-1}` Scale value of the water pressure(m)
        Braud_n = 2.22
        Braud_m = 0.099
        Braud_theta_r = 0.01
        Braud_eta = 9.14

        if self.retention == 'VGM':
            r_curve = cmf.VanGenuchtenMualem(Ksat=Braud_Ksat, phi=Braud_phi, alpha=0.005, n=Braud_n,
                                             theta_r=Braud_theta_r)
            # Oversaturation tolerence upto 1% for matrix pot = +1
            r_curve.w0 = 0.9995  # See: https://philippkraft.github.io/cmf/cmf_tut_retentioncurve.html Oversaturation
            r_curve.l = 0.67  # tortuosity
            r_curve.K = r_curve.Ksat * pow(r_curve.w0, Braud_eta)

        elif self.retention == 'BC':
            r_curve = cmf.BrooksCoreyRetentionCurve(ksat=0.0106272, porosity=0.35, _b=1 / 0.22, theta_x=0.3499)
            r_curve.w0 = 0.9995
            # r_curve.residual_theta

        else:
            return NotImplementedError

        return r_curve

    def add_layers(self, c, l_boundaries):

        for d in l_boundaries:
            l = c.add_layer(d, self.rtn_curve)
            # Check tehts exceeds porosity if  l.wetness = 1.0
            # check for wetness = 1, saturated depth = 0

    def add_connections(self, cell, connection=cmf.Richards):

        cell.install_connection(connection)

    def cmf_boundary(self, c, et_pot=0.55):

        vg = c.vegetation

        # vg.Height = 0
        # vg.LAI = 0
        vg.RootDepth = 0.0
        vg.fraction_at_rootdepth = 0.0

        sw = cmf.ShuttleworthWallace
        c.install_connection(sw)

        # Create the boundary condition
        gw = self.project.NewOutlet('groundwater', x=0, y=0, z=-1)
        # Set the potential
        gw.potential = 0
        # Connect the lowest layer to the groundwater using Richards percolation
        gw_flux = cmf.Richards(c.layers[-1], gw)

    def cmf_setup(self, l_boundaries):
        # CMF setup
        p = self.project
        c = p.cells[0]

        # L_boundaries = np.arange(0.01, 1.01, 0.05)  # define layer thicknesses
        self.add_layers(c, l_boundaries)  # define retention curve to all cell layers

        c.saturated_depth = 0  # Set all layers to a potential of 0 m
        c.surfacewater.depth = 0.0  # 0 mm water in the surface water storage for percolation

        self.add_connections(c, connection=cmf.Richards)  # apply layer connection

        return p

    def update_theta(self, c_cmf, c_vapor):

        lcmf, lvp = c_cmf.layers, c_vapor.layers

        # update corresponding head based on theta
        if len(lcmf) == len(lvp):
            for lcm, lv in zip(lcmf, lvp):

                if not lv.theta >= lv.theta_sat:

                    lcm.theta = lv.theta

        else:
            return ValueError("Length of layers and _theta mismatch.")


###### vapor ##########
class vapor_iso(object):

    def __init__(self, cell_cmf):

        self.cell_cmf = cell_cmf
        self.__p = None

    def run(self, initial_dt=1, simulation_time=3600, model='proposed'):

        T = 0
        while T < simulation_time:
            T += initial_dt
            dt = initial_dt
            # print('T: ', T)
            initial_dt = self.project.run(dt=dt, dt_max=60, epsilon=1e-4, model=model,
                                          water_flux=False, heat_flux=True, vapor_flux=False)

    @property
    def project(self):
        return self.__p

    def set_project(self):

        p = vp.project()
        self.__p = p
        self.__p.new_cell(atmosphere=self._atm())

        return self.project

    def _atm(self):

        atm = vp.storage.atmosphere(T=303.1600036621094,  # temperature in the atmosphere in Kelvin
                                    Rh_atmosphere=0.20,  # relative humidity of the atmosphere (-)
                                    Pa_atmosphere=10 ** 5,  # atmospheric pressure (Pa)
                                    initial_wind_speed=2.0,  # wind speed at top of canopy in m/s
                                    R_net=0.0,  # net radiation absorbed
                                    wind_speed=2.0,  # wind speed at top of canopy in m/s
                                    hc=10,  # canopy height [m] (e.g. 40)
                                    d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                                    z0m=0.1 * 10,
                                    # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                                    LAI=1.1,  # leaf area index (e.g. 2.0)
                                    extku=1.5)  # extinction coeff't for windspeed (e.g. 0.5))

        return atm

    def _layers(self, c):

        tmp_ini = 30
        for l in self.cell_cmf.layers:
            new_layer = vp.storage.soil_layer(upper_boundary=l.upper_boundary * 100,  # to cm
                                              lower_boundary=l.lower_boundary * 100,  # to cm
                                              theta=l.theta,  # 0.35,  # 0.2
                                              theta_sat=l.porosity,  # 0.547,
                                              tortuosity=0.67,
                                              rH=0.2,
                                              head=0,
                                              T=30,
                                              T0=30,
                                              K_s=l.Ksat * 100 / 86400,  # m day-1 to cm s-1
                                              # 3.8e-3, #1.23e-5,  # cm /s [braud et.al, cuntz i.e SISPAT / SLI]
                                              phi=-19.3,  # cm  # air entry potential
                                              lam=0.22,  # 0.15313935681470137, #0.22,
                                              eta=9.14,  # braut test cases
                                              # shape coeffecient Brooks-Corey (1964) [value = 0.22 from Han Fu et.al. (2023)
                                              clay=0.5,  # 0.249, # 0.5,
                                              sand=0.2,  # 0.021,  # 0.2,  # clay-sand-silt from sli.utils init params
                                              silt=0.3,  # 0.73,  # 0.3,
                                              som=0.044,
                                              rho_b=1.4,  # 1.4,  # 1.2,  # sli init params
                                              rho_l=1.0
                                              )

            new_layer.head = new_layer.calc_head(l.theta)
            new_layer.T = tmp_ini

            c.add_layer(new_layer)

    def soil_vapor_setup(self):

        p = self.set_project()
        c = p.cells[0]

        self._layers(c)  # add layers

        #####Install connections#######
        c.install_connections(water=False, heat=True, vapor=False)
        c.add_evaporation(top_layer=c.layers[0])

        # boundaries
        c.hleft, c.hright = c.layers[0].head, c.layers[-1].head
        c.Tleft, c.Tright = 30, c.layers[-1].T

    def update_theta(self, c_cmf, c_vapor):

        lcmf, lvp = c_cmf.layers, c_vapor.layers

        # update corresponding head based on theta
        if len(lcmf) == len(lvp):
            for lcm, lv in zip(lcmf, lvp):
                lv.head = lv.calc_head(lcm.theta)
        else:
            return ValueError("Length of layers and _theta mismatch.")

    def update_evap(self, c, qev):

        c.update_evaporation(qev)


####### Isotope #######
class isotope(object):

    def __init__(self, p_iso, p_cmf, p_vapor, testcase=1):

        self.pcmf = p_cmf
        self.pvapor = p_vapor
        self.piso = p_iso

        self.c_cmf = self.pcmf.cells[0]
        self.c_vapor = self.pvapor.cells[0]

        self.test_case = testcase

    @property
    def c_iso(self):
        return self.piso.get_cells()[0]

    def run_iso(self, solutes=["2H", "18O"], dt=60):

        try:
            kwargs = iso_delta.test_case_args(self.test_case)

            for solute in solutes:
                dc = self.piso.run(Isotopologue=solute, delta_time=dt, error_tol=1e-9, **kwargs)

                c = list(np.array(self.c_iso.get_conc_layers(Isotopologue=solute)) + np.array(dc))
                # update iso concentrations to current time step
                self.c_iso.update_c_layers(conc_iso=c, Isotopologue=solute)

        except ValueError:
            raise NotImplementedError

    def iso_soil_setup(self):

        self.piso.new_cell(atmosphere=self.iso_atm(), area=1, x=0, y=0, z=0)  # add new cell
        self.i_layers()  # add layers to the current cell

        #####Install connections#######
        self.c_iso.install_connections()  # install storage connections between the layers

        # boundary connections
        self.c_iso.add_evaporation(layer=self.c_iso.layers[0])

    def iso_atm(self):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        # atmospheric variables
        patm = 1  # from SLI_solve
        Tatm = 30
        Rh_atm = 0.2
        wind_speed = 2  # sli.wind_speed(dt)

        ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=self.test_case)
        ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=self.test_case)
        c_iso_2H = iso_delta.delta_to_concentration(atm_2H, '2H')
        c_iso_18O = iso_delta.delta_to_concentration(atm_18O, '18O')

        atm = iso_atmosphere(conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                             conc_iso_vapor={"2H": c_iso_2H, "18O": c_iso_18O},
                             T=Tatm + Tzero_sli, Rh_atmosphere=Rh_atm,
                             Pa_atmosphere=patm,
                             wind_speed=wind_speed,
                             hc=10,  # canopy height [m] (e.g. 40)
                             d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                             z0m=0.1 * 10,
                             # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                             LAI=1.1,  # leaf area index (e.g. 2.0)
                             extku=1.5)

        return atm

    def i_layers(self):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        l_cmf, l_vap = self.c_cmf.layers, self.c_vapor.layers

        ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=self.test_case)
        ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=self.test_case)
        init_c_iso_2H = iso_delta.delta_to_concentration(ali_2H, '2H')
        init_c_iso_18O = iso_delta.delta_to_concentration(ali_18O, '18O')

        # ignoring boundary layers / right layer / top layer
        for id, (lc, lv) in enumerate(zip(l_cmf, l_vap)):
            new_layer = iso_storages.iso_soil_layer(ID=id,
                                                    upper_boundary=lc.upper_boundary,
                                                    lower_boundary=lc.lower_boundary,
                                                    conc_iso_liquid={"2H": init_c_iso_2H, "18O": init_c_iso_18O},
                                                    theta=lc.theta,
                                                    theta_0=0.01,
                                                    theta_sat=lv.theta_sat,
                                                    tortuosity=lv.tortuosity,
                                                    T=lv.T + Tzero_sli,
                                                    rH=0.9,
                                                    psi=lc.potential
                                                    )

            self.c_iso.add_layer(new_layer)

    def update_iso_storages(self):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        # Need not update atmosphere
        """"
        # atmospheric variables
        patm = 1  # from SLI_solve
        Tatm = sli.Tatm(dt) + Tzero_sli
        Rh_atm = sli.R_humidity_atm(dt)
        c_iso_18O = sli.civa(dt) / sli.cva(dt)
        c_iso_2H = 0.15367056287906838
        wind_speed = sli.wind_speed(dt)
        c.update_atmosphere(c_atm={'2H': c_iso_2H, '18O': c_iso_18O}, T=Tatm, Rh=Rh_atm, Pa=patm, wind_speed=wind_speed)
        """
        l_cmf, l_vap = self.c_cmf.layers, self.c_vapor.layers

        # Update layers
        theta = [l.theta for l in l_cmf]
        T_soil = [l.T + Tzero_sli for l in l_vap]
        rH = None
        psi = [l.potential / 100 for l in l_cmf]

        self.c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)
        # self.c_iso.update_aquifer(c_iso={"2H": 0.0, "18O": 0.0})

    def update_iso_boundaries(self, t):

        # Update states and fluxes for

        # c: cell
        # sli: soil_litter_iso class to import input variables
        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        # Surface variables
        T_surface = self.c_vapor.layers[0].T + Tzero_sli
        q_ev = self.c_cmf.evaporation(t) / 1000 / 86400  # mm day into m s-1
        qv_surface = 0 #-self.c_vapor.connection_evap.q_vapor / 100  # cm s-1 to m s-1
        ql_surface = -q_ev
        self.c_iso.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

        # soil fluxes also ignores the flux from top layer / right layer to the boundary layer
        l_fluxes = [l1.flux_to(l2, t) for l1, l2 in zip(self.c_cmf.layers[:-1], self.c_cmf.layers[1:])]
        ql = np.array(l_fluxes) / 86400000  # m3 day-1 into ms-1
        qv = 0#np.array(self.c_vapor.vapor_fluxes) / 100  # cm s-1 to m s-1
        ql, qv = np.append(ql, 0.0), np.append(qv, 0.0)

        self.c_iso.update_liquid_fluxes(liquid_fluxes=ql)
        self.c_iso.update_vapor_fluxes(vapor_fluxes=None)  # None if qv is computed internally, else list of qv, len = len(layers)


########### cmf #############
_cmf = cmf_iso(r_curve='VGM')
L_boundaries = np.arange(0.05, 1.05, 0.05)
pcmf = _cmf.cmf_setup(l_boundaries=L_boundaries)  # cmf project
ccmf = pcmf.cells[0]  # cmf cell
_cmf.cmf_boundary(c=ccmf)  # 1.005 * 10 ** - 5 kg m-2 s-1 to mm day-1

########### vapor #############
_vapor = vapor_iso(cell_cmf=ccmf)
_vapor.soil_vapor_setup()  # vapor module
pvapor = _vapor.project  # vapor project
cvapor = pvapor.cells[0]  # vapor cell

_v = vp.cmf_vapor_interface  # interface module to add vapor fluxes to cmf
vap = _v.Vaporizer(ccmf)

######### iso ###########
piso = iso_project()
_iso = isotope(p_iso=piso, p_cmf=pcmf, p_vapor=pvapor, testcase=1)
_iso.iso_soil_setup()
ciso = _iso.c_iso

# Define solver
solver = cmf.CVodeBanded(pcmf)
start = datetime.datetime(2024, 1, 1)
end = start + timedelta(days=250)
timestep = timedelta(hours=1)

ev_act, trp_act = [], []
q_liquid, q_vapor = [], []
theta, T = [], []
delta = []
count = 1
f = 1 / 8600 / 1000  # factor to convert [mm day-1] (m3 day-1) to [ms-1]
for t in solver.run(start, end, timestep):

    print(count, ' ', t)

    f = 1 / 8600 / 1000  # factor to convert [mm day-1] (m3 day-1) to [ms-1]

    # run vapor model
    _vapor.update_evap(cvapor, qev=ccmf.evaporation(t) * f * 100)  # since q in vapor model is in cm s-1
    _vapor.update_theta(c_cmf=ccmf, c_vapor=cvapor)  # update vapor with new theta
    _vapor.run(initial_dt=60, simulation_time=timestep.seconds)

    # update and run iso
    _iso.update_iso_storages(), _iso.update_iso_boundaries(t)
    _iso.run_iso(solutes=["2H"], dt=timestep.seconds)

    print(ciso.conc_2H_delta)
    print([l.T - 273.16 for l in ciso.layers])

    qv = np.array(ciso.vapor_fluxes) * 86400000  # ms-1 to m3 day-1
    qvs = ccmf.evaporation(t)
    qv_layers = np.concatenate(([qv[0]], qv[:-1] - qv[1:], [qv[-1]]))
    vap.flux = np.append(qv, 0)

    vap.flux = qv

    count += 1
