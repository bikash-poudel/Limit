'''
Created on 20.12.2024
@author: poudel-b
'''
import numpy as np

# -*- coding: utf-8 -*-

from src import *
from vapor_model import *

from math import sin, pi
import matplotlib.pyplot as plt


class soil_iso(object):

    def __init__(self,
                 soil_project,
                 i_project,
                 iso_testcase=1):

        self.s_project = soil_project
        self.i_project = i_project
        self.test_case = iso_testcase
        self.__qevap = []

    @property
    def c_soil(self):
        return self.s_project.cells[0]

    @property
    def c_iso(self):
        return self.i_project.get_cells()[0]

    @property
    def evap(self):
        return self.__qevap

    """"Functions"""

    def run_soil(self, model='coupled', initial_dt=60, simulation_time=86400 * 250, epsilon=1e-4,
                 run_iso=True, solutes=["2H", "18O"]):

        count = 0
        T = 0
        while T < simulation_time:

            count += 1

            print(count, ' ', 'dt: ', initial_dt, 's ', T / 86400, ': days')
            #print(self.c_iso.conc_2H_delta)
            #print([l.theta for l in self.c_soil.layers])

            T += initial_dt
            dt = initial_dt

            self.update_boundaries(T)
            initial_dt = self.s_project.run(dt=dt, dt_max=60, epsilon=epsilon, model=model)

            if run_iso:
                # update storage states and boundaries to current time
                self.update_iso_storages(), self.update_iso_boundaries()
                if count > 1:  # starting from second time step
                    self.run_iso(solutes=solutes, dt=dt)

    def run_iso(self, solutes=["2H", "18O"], dt=60):

        try:
            kwargs = iso_delta.test_case_args(self.test_case)

            for solute in solutes:
                dc = self.i_project.run(Isotopologue=solute, delta_time=dt, error_tol=1e-9, **kwargs)

                c = list(np.array(self.c_iso.get_conc_layers(Isotopologue=solute)) + np.array(dc))
                self.c_iso.update_c_layers(conc_iso=c,
                                           Isotopologue=solute)  # update iso concentrations to current time step

        except ValueError:
            raise NotImplementedError

    def iso_soil_setup(self):

        self.s_project.new_cell(self._atm())
        self.i_project.new_cell(atmosphere=self.iso_atm(), area=1, x=0, y=0, z=0)  # add new cell

        self._layers(self.c_soil)
        self.i_layers(ciso=self.c_iso, csoil=self.c_soil)  # add layers to the current cell

        #####Install connections#######
        self.c_soil.install_connections()
        self.c_iso.install_connections()  # install storage connections between the layers

        # boundary connections
        self.c_soil.add_evaporation(top_layer=self.c_soil.layers[0])
        self.c_iso.add_evaporation(layer=self.c_iso.layers[0])

        # -9.274162934176326e+03, 8.570390159801457e+05
        self.c_soil.hleft, self.c_soil.hright = 8.570390159801457e+05, -19.3
        self.c_soil.Tleft, self.c_soil.Tright = 30, 30

        aq = iso_storages.iso_aquifer(
            conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # define aquifer as boundary isotope storage
        self.c_iso.add_aquifer(aq, self.c_iso.layers[-1])  # aquifer connected to bottom layer

    def _atm(self):

        atm = storage.atmosphere(T=303.1600036621094,  # temperature in the atmosphere in Kelvin
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

    def iso_atm(self):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        # atmospheric variables
        patm = 1000  # from SLI_solve
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

    def _layers(self, c):

        dx = [0.009999999776482582, 0.009999999776482582, 0.030000001192092896, 0.05000000074505806,
              0.05000000447034836, 0.04999999701976776, 0.04999999701976776, 0.050000011920928955,
              0.09999999403953552, 0.20000001788139343, 0.41999995708465576, 0.24800002574920654,
              0.29499995708465576, 0.35199999809265137, 0.42100000381469727, 0.500999927520752,
              0.5980000495910645, 0.7140002250671387, 0.8509998321533203]

        #lower_boundaries = np.cumsum(dx)*100

        lower_boundaries = np.arange(0.5, 50.5, 0.5)
        #lower_boundaries = np.insert(depth, 0, 0.05)

        theta_init, tmp_ini = 0.35, 30

        id = 0
        upper_boundary = 0
        for lower_boundary in lower_boundaries:
            new_layer = storage.soil_layer(upper_boundary=upper_boundary,
                                           lower_boundary=lower_boundary,
                                           theta=0.5, #0.35,  # 0.2
                                           theta_sat=0.5,  # 0.547,
                                           tortuosity=0.67,
                                           rH=0.2,
                                           head=-19.3,
                                           T=30,
                                           T0=30,
                                           K_s=1.23e-5, # 3.8e-3, #1.23e-5,  # cm /s [braud et.al, cuntz i.e SISPAT / SLI]
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

            new_layer.head = new_layer.calc_head(theta_init)
            new_layer.T = tmp_ini

            id += 1
            upper_boundary = lower_boundary
            c.add_layer(new_layer)

    def i_layers(self, ciso, csoil):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        soil_layers = csoil.layers

        ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=self.test_case)
        ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=self.test_case)
        init_c_iso_2H = [iso_delta.delta_to_concentration(ali_2H, '2H')] * len(soil_layers)
        init_c_iso_18O = [iso_delta.delta_to_concentration(ali_18O, '18O')] * len(soil_layers)

        to_m = 1 / 100  # cm to m
        # ignoring boundary layers / right layer / top layer
        id = 0
        for layer, c2H, c18O in zip(soil_layers, init_c_iso_2H, init_c_iso_18O):
            new_layer = iso_storages.iso_soil_layer(ID=id,
                                                    upper_boundary=layer.upper_boundary * to_m,
                                                    lower_boundary=layer.lower_boundary * to_m,
                                                    conc_iso_liquid={"2H": c2H, "18O": c18O},
                                                    theta=layer.theta,
                                                    theta_0=0.01,
                                                    theta_sat=layer.theta_sat,
                                                    tortuosity=layer.tortuosity,
                                                    T=layer.T + Tzero_sli,
                                                    rH=layer.rH,
                                                    psi=layer.head * to_m
                                                    )
            id += 1
            ciso.add_layer(new_layer)

    def update_boundaries(self, T):

        q_ev = 1e-5 + 5e-6 * sin(2 * pi * T / 86400)   # cm / s
        self.c_soil.update_evaporation(q_ev)
        self.__qevap.append(self.c_soil.q_evap)

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
        layers = self.c_soil.layers

        # Update layers
        theta = [l.theta for l in layers]
        T_soil = [l.T + Tzero_sli for l in layers]
        rH = None
        psi = [l.head / 100 for l in layers]

        self.c_iso.update_layers(theta=theta, T=T_soil, rH=rH, psi=psi)
        self.c_iso.update_aquifer(c_iso={"2H": 0.0, "18O": 0.0})

    def update_iso_boundaries(self):

        c = self.c_iso
        # Update states and fluxes for

        # c: cell
        # sli: soil_litter_iso class to import input variables
        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        # Surface variables
        T_surface = self.c_soil.layers[0].T + Tzero_sli
        q_ev = self.c_soil.q_evap / 100
        ql_surface = -q_ev
        qv_surface = 0
        c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

        # soil fluxes also ignores the flux from top layer / right layer to the boundary layer
        ql, qv = [q / 100 for q in self.c_soil.liquid_fluxes], [q / 100 for q in self.c_soil.vapor_fluxes]  # cm to m
        ql.insert(-1, 0), qv.insert(-1, 0)  # additional flux at left boundary
        c.update_liquid_fluxes(liquid_fluxes=ql)
        c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

        # boundary storage
        c.update_connection_to_aquifer(ql_layer=0)


s = soil_iso(soil_project=project(),
             i_project=iso_project(),
             iso_testcase=1)

s.iso_soil_setup()
s.run_soil(model='proposed', simulation_time=86400 * 250, run_iso=False, solutes=['2H'])

s_layers, i_layers = s.c_soil.layers, s.c_iso.layers
depth = [-l.center.z for l in i_layers]
head = [l.head / 100 for l in s_layers]
theta = [l.theta for l in s_layers]
T = [l.T for l in s_layers]

c2h = s.c_iso.conc_2H_delta
c18O = s.c_iso.conc_18O_delta

ev = s.evap
ql = [q / 100 for q in s.c_soil.liquid_fluxes]
qv = [q / 100 for q in s.c_soil.vapor_fluxes]

plt.plot(c2h, depth, label='2H')
plt.legend()
plt.show()

plt.plot(c18O, depth, label='18O')
plt.legend()
plt.show()

plt.plot(ev, label='evaporation')
plt.legend()
plt.show()

plt.plot(theta, depth, label='theta')
plt.legend()
plt.show()

plt.plot(ql, depth[:-1], label='q_liquid')
plt.plot(qv, depth[:-1], label='q_vapor')
plt.legend()
plt.show()

