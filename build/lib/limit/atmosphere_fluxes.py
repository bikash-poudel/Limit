# -*- coding=utf-8 -*-
from math import exp, sqrt, log

from . heat_fluxes import heat_node


class flux_atmosphere(heat_node):

    def __init__(self,
                 atmosphere,
                 top_layer
                 ):

        super().__init__(top_layer=top_layer)

        self.atmosphere = atmosphere
        self.qevap = None

        self.surface_layer = top_layer
        self.Ta = self.atmosphere.T
        self.rha = self.atmosphere.Rh
        self.Rnet = self.atmosphere.R_net
        self.wind_speed = self.atmosphere.wind_speed
        self.extku = self.atmosphere.extku
        self.LAI = self.atmosphere.LAI
        self.hc = self.atmosphere.hc
        self.d0 = self.atmosphere.d0
        self.z0m = self.atmosphere.z0m

        self.Mw = 1.80159993E-02  # molecular wt: water (kg/mol)
        self.rgas = 8.31429958  # universal gas const  (J/mol/K)
        self.rmair = 2.89699994E-02  # molecular wt: dry air (kg/mol)
        self.cpa = 1004.6400146484375  # air spec. heat capacity (J/kg/K)

    def T_surface(self):

        # soil-atm interface temperature
        Tzero = 273.16000366210938
        dz = 0.5 * self.surface_layer.thickness
        dh = self.dh()
        gh = self.gh()

        T0 = (- dz * self.lE0() +
             dz * self.Rnet +
             dh * (self.surface_layer.T - Tzero) +
             dz * gh * (self.Ta - Tzero)) / (dh + dz * gh)

        return T0

    def q_evaporation(self):
        return self.lE0() / self.lambdav(self.surface_layer) / 1000

    def lE0(self):

        return min(self.ev_pot(), self.E_total())
        # return self.qevap * self.lambdav(self.surface_layer) * 1000

    def E_total(self):

        return self.E_vapor() + self.E_liquid()

    def E_vapor(self):

        return self.qv_evap() * self.lambdav(self.surface_layer)

    def E_liquid(self):

        return self.ql_evap() * self.lambdav(self.surface_layer)

    def qv_evap(self):

        rh = self.rh(self.surface_layer)
        csat = self.csat(self.surface_layer.T)
        cva = self.cva() * 1000
        rbw = self.rbw()
        dv = self.dv(self.surface_layer)
        dz = 0.5 * self.surface_layer.thickness

        if not dv == 0:
            return (rh * csat - cva) / (rbw + dz / dv)
        else:
            return 0

    def ql_evap(self):

        phimin = self.phiz(self.surface_layer, h=self.hmin)
        phi = self.phi(self.surface_layer)

        if self.repeat:
            K = self.KZ(self.surface_layer, h=self.hmin)
        else:
            K = self.K(self.surface_layer)

        dz = 0.5 * self.surface_layer.thickness

        return ((phi - phimin) / dz - K ) * 1000

    def ev_pot(self):

        lambdav = self.lambdav(self.surface_layer)
        rhocp = self.rhocp()
        gamma = self.gamma(lambdav=lambdav)
        s = self.s_esat(self.Ta)
        ea = self.esat(self.Ta) * self.rha
        k = self.kth(self.surface_layer)
        rbh = self.rrc()
        rbw = self.rbw()
        dz = 0.5 * self.surface_layer.thickness
        Da = ea / self.rha - ea

        E = (rhocp * (Da * (k * rbh + dz * rhocp) +
                      rbh * s * (dz * self.Rnet + k * (- self.Ta + self.surface_layer.T)))) /\
            (gamma * rbw * (k * rbh + dz * rhocp) + dz * rbh * rhocp * s)

        return E

    def dEdTsoil(self):

        rhocp = self.rhocp()
        rbh = self.rrc()
        rbw = self.rbw()
        s = self.s_esat(self.Ta)
        lambdav = self.lambdav(self.surface_layer)
        gamma = self.gamma(lambdav=lambdav)
        k = self.kth(self.surface_layer)
        dz = 0.5 * self.surface_layer.thickness

        return - ((k * rbh * rhocp * s) / (gamma * k * rbh * rbw +
                                           dz * gamma * rbw * rhocp + dz * rbh * rhocp * s))

    def dE_vapT1(self):

        dz = self.surface_layer.thickness * 0.5
        dv = self.dv(self.surface_layer)
        lambdav = self.lambdav(self.surface_layer)
        rh = self.rh(self.surface_layer)
        sc = self.s_csat(self.surface_layer.T)

        return (rh * sc) / (1 / self.gv() + dz / dv) * lambdav

    def cva(self):

        return self.rha * self.esat(self.Ta) * 0.018 / (1000 * 8.314 * self.Ta)

    def s_csat(self, T):

        Tzero = 273.16000366210938
        csat = self.csat(T)

        return csat * 17.270000457763672 * 237.30000305175781 / (((T - Tzero) + 237.30000305175781) ** 2)

    def s_esat(self, T):

        Tzero = 273.16000366210938

        return self.esat(T) * 17.270000457763672 * 237.30000305175781 / (((T - Tzero) + 237.30000305175781) ** 2)

    def csat(self, T):

        Mw = 1.8015999346971512E-002

        return self.esat(T) * Mw / self.rgas / T

    def esat(self, T):

        Tzero = 273.16000366210938

        esat = 610.59999465942383 * exp(17.270000457763672 * (T - Tzero) / ((T - Tzero) + 237.30000305175781))

        return esat

    def G0(self):

        Tzero = 273.16000366210938
        dz = 0.5 * self.surface_layer.thickness

        return self.dh() / dz * (self.T_surface() - (self.surface_layer.T - Tzero))

    def dh(self):
        return self.kth(self.surface_layer)

    def gh(self):

        rhocp = 1189.8000488281250    # cpa * rhoa at std(25 degC) [J / m3K]

        return 1 / self.rrc() * rhocp

    def rrc(self):

        # combined convective and radiative resistance to heat transfer
        return 1 / self.grc()

    def grc(self):

        # combined convective and radiative conductivity to heat transfer
        rhocp = 1189.8000488281250

        return 1 / self.rbh() + self.gr() / rhocp

    def gr(self):

        # radiation conductance Wm-2K-1
        Tzero = 273.16000366210938
        e_soil = 0.97000002861022949  # soil emissivity
        return 4 * e_soil * (self.Ta - Tzero + 273.15) ** 3 * 5.67e-8

    def rbh(self):

        # resistance to heat transfer beween surface and lowest atmospheric layer
        return self.rbw()

    def rbw(self):

        # esistance to water vapour transfer beween surface and lowest atmospheric layer
        return self.ra() + self.rs()

    def rs(self):

        # soil laminar boundary layer resistance
        return 1 / (0.004 + 0.012 * self.us())

    def us(self):

        # windspeed below canopy
        return self.wind_speed * exp(-self.extku * self.LAI)

    def ra(self):

        # aerodyamic resistance( from z0m to hc)
        return log((self.hc - self.d0) / self.z0m) ** 2 / 0.41 ** 2 / self.wind_speed

    def rh(self, node):

        gravity = 9.8000001907348633  # gravity acceleration (m/s2)
        rhmin = 5.0000000745058060E-002  # min relative humidity SLI_utils:L1409
        Mw = 1.8015999346971512E-002

        if node.S < 1:
            hr = max(exp(Mw * gravity * self.h(node) / (self.rgas * node.T)), rhmin)
        else:
            hr = 1

        return hr

    def gv(self):
        return 1 / self.rbw()

    def rhocp(self):

        rmair = 0.028969999382339233  # molecular wt: dry air (kg/mol)
        R = 8.3142995834350586  # universal gas constant (j / mol / k)

        return rmair * 101325 / self.Ta * self.cpa / R
        # return 1169.9788141999604

    def gamma(self, lambdav):

        rmair = 0.02896999916806393
        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)

        return 101325 * self.cpa / lambdav / (Mw / rmair)
        # return 67.287029307028789

    def qyb_qTb(self):

        dz = 0.5 * self.surface_layer.thickness
        lambdav = self.lambdav(self.surface_layer)

        if self.ev_pot() > self.E_total():
            qTb = - self.dE_vapT1() / (1000 * lambdav)
            qyb = - (self.phiS(self.surface_layer) / dz - self.KS(self.surface_layer))
        else:
            qTb = - self.dEdTsoil() / (1000 * lambdav)
            qyb = 0.0

        return qyb, qTb

    def qlyb(self):

        if self.ev_pot() > self.E_total():

            return - (self.phiS(self.surface_layer) / (0.5 * self.surface_layer.thickness)
                      - self.KS(self.surface_layer))

        else:

            return 0

    def qvyb(self):

            return 0

    def qlTb(self):
        return 0


