
'''
Created on 06.08.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-
from math import exp

import iso_cell
import iso_sli
import iso_storages

from iso_fluxes import liquid_diffusion_base_class


# TODO: recheck this test
class Barnes_Alison(object):

    def __init__(self,
                 cell):

        """

        Barnes and Allison (1983) solution for steady state analytical solution:
        for isothermal saturated case
        for non-iso thermal unsaturated case

        """
        self.cell = cell
        self.__q_evap = None

    @property
    def layers(self):
        return self.cell.layers

    @property
    def connection_evap(self):

        return self.cell.connection_evap

    @property
    def q_evap(self):
        return self.connection_evap.q_evap

    def isothermal(self, Isotopologue):

        ld = liquid_diffusion_base_class()

        rhow = 1000  # kg m-3
        c_s = self.connection_evap.c_iso_liq_surface(Isotopologue=Isotopologue, BA83=True)  # surface concentration

        ci_liq = []
        for l in self.layers:

            c_ali = l.get_conc_iso_liquid(Isotopologue)  # conc alimentation water in soil
            #dl_i = ld.dl_i(l.T, Isotopologue, l.theta, l.tortuosity)
            dl_i = ld.dl_i_self_diffusivity(l.T, Isotopologue=Isotopologue)
            z = l.center.z

            c = c_ali + (c_s - c_ali) * exp(- self.q_evap * z / dl_i)
            ci_liq.append(c)

        return ci_liq


def update_boundaries(c, sli, dt):
    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # Surface variables
    T_surface = sli.Ts(dt) + Tzero_sli
    q_ev = sli.qevap(dt)
    ql_surface = sli.ql0(dt)
    qv_surface = sli.qv0(dt)
    c.update_evaporation(q_ev=q_ev, T_surface=T_surface, ql_surface=ql_surface, qv_surface=qv_surface)

    # soil fluxes
    ql = sli.qlsig(dt)[1:]
    qv = sli.qvsig(dt)[1:]
    c.update_liquid_fluxes(liquid_fluxes=ql)
    c.update_vapor_fluxes(vapor_fluxes=qv)  # None if qv is computed internally, else list of qv, len = len(layers)

    # boundary storage
    c.update_connection_to_aquifer()


def setup_iso(sli, testcase):

    # Define boundary storages
    location = iso_storages.Point(x=0, y=0, z=0)
    atm = iso_sli._atm(sli, testcase=7)  # atmosphere
    aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # define aquifer as boundary isotope storage

    c = iso_cell.iso_cell(location, atm)
    iso_sli._layers(c, sli, testcase=testcase),  # add layers to the current cell

    #####Install connections#######
    c.install_connections()  # install storage connections between the layers
    # boundary connections
    c.add_evaporation(), c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

    return run_BA84(sli, c)


def run_BA84(sli, c):

    solutes = ["2H", "18O"]
    BA84 = Barnes_Alison(cell=c)

    c_iso, c_iso_delta = {'2H': [], '18O': []}, {'2H': [], '18O': []}
    for dt in range(1, len(sli.get_in_soil()) - 1):

        print(dt)
        c_iso["2H"].append(c.conc_2H), c_iso["18O"].append(c.conc_18O)
        c_iso_delta["2H"].append(c.conc_2H_delta), c_iso_delta["18O"].append(c.conc_18O_delta)

        update_boundaries(c, sli, dt)
        for solute in solutes:

            c_liq = BA84.isothermal(Isotopologue=solute)
            c.update_c_layers(conc_iso=c_liq, Isotopologue=solute)  # update iso concentrations to current time step

    return c_iso_delta


#test = 7
#sli = iso_sli.get_sli(testcase=test)
#delta = setup_iso(sli, testcase=test)


