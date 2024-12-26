'''
Created on 20.12.2024
@author: poudel-b
'''

# -*- coding: utf-8 -*-

from src import *
from vapor_model import *

from math import sin, pi


class soil_iso(object):

    def __init__(self,
                soil_project,
                i_project,
                iso_testcase=1):

        self.s_project = soil_project
        self.i_project = i_project
        self.test_case = iso_testcase

    @property
    def c_soil(self):
        return self.s_project.cells[0]

    @property
    def c_iso(self):
        return self.i_project.get_cells()[0]


    """"Functions"""

    def iso_soil_setup(self):

        self.s_project.new_cell(self._atm())
        self.i_project.new_cell(atmosphere=self.iso_atm(), area=1, x=0, y=0, z=0)  # add new cell

        self._layers(self.c_soil)
        self.i_layers(ciso=self.c_iso, csoil=self.c_soil)  # add layers to the current cell

        #####Install connections#######
        self.c_soil.install_connections()
        self.c_iso.install_connections()  # install storage connections between the layers

        # boundary connections
        self.c_soil.add_evaporation(top_layer=self.c_soil.layers[-1])
        self.c_iso.add_evaporation(layer=self.c_iso.layers[-1])

        self.c_soil.hleft, self.c_soil.hright = -9.274162934176326e+03, 8.570390159801457e+05
        self.c_soil.Tleft, self.c_soil.Tright = 25, 30

    def run(self, initial_dt=1, simulation_time=86400*250, epsilon=1e-4):

        T = 0
        while T < simulation_time:

            dt = initial_dt
            print(T, ': days')
            q_ev = 1e-7 + 2e-8 * sin(2 * pi * (T + dt) / 86400)  # cm / s
            self.c_soil.q_evap = q_ev

            initial_dt = self.s_project.run(dt=dt, dt_max=60, epsilon=epsilon)

            T += initial_dt

    def _atm(self):

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

    def iso_atm(self):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        # atmospheric variables
        patm = 1  # from SLI_solve
        Tatm = 30
        Rh_atm = 0.2
        wind_speed = 2 # sli.wind_speed(dt)

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

        lower_boundaries = np.arange(2.5, 54, 2.5)

        theta_init = 0.2
        tmp_ini = 25

        id = 0
        upper_boundary = 0
        for lower_boundary in lower_boundaries:
            new_layer = storage.soil_layer(upper_boundary=upper_boundary,
                                           lower_boundary=lower_boundary,
                                           theta=0.2,
                                           theta_sat=0.547,
                                           tortuosity=0.67,
                                           rH=0.2,
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

    def i_layers(self, ciso, csoil):

        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

        ali_2H, atm_2H = iso_delta.delta_testcases('2H', testcase=self.test_case)
        ali_18O, atm_18O = iso_delta.delta_testcases('18O', testcase=self.test_case)
        init_c_iso_2H = [iso_delta.delta_to_concentration(ali_2H, '2H')] * len(csoil.layers)
        init_c_iso_18O = [iso_delta.delta_to_concentration(ali_18O, '18O')] * len(csoil.layers)

        id = 0
        for layer, c2H, c18O in zip(csoil.layers, init_c_iso_2H, init_c_iso_18O):

            new_layer = iso_storages.iso_soil_layer(ID=id,
                                                    upper_boundary=layer.upper_boundary,
                                                    lower_boundary=layer.lower_boundary,
                                                    conc_iso_liquid={"2H": c2H, "18O": c18O},
                                                    theta=layer.theta,
                                                    theta_0=0.01,
                                                    theta_sat=layer.theta_sat,
                                                    tortuosity=layer.tortuosity,
                                                    T=layer.T + Tzero_sli,
                                                    rH=layer.rH,
                                                    psi=layer.head
                                                    )
            id += 1

            ciso.add_layer(new_layer)


s = soil_iso(soil_project=project(),
             i_project=iso_project(),
             iso_testcase=1)

s.iso_soil_setup()
s.run(simulation_time=86400)

