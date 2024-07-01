'''
Created on 03.06.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

import os
import numpy as np

import Sli
import iso_project
import iso_storages
import iso_fluxes


def _layers(c, sli):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    dt = 0  # initial states
    lower_boundaries = np.cumsum(sli.dx(dt))
    init_c_iso = sli.ciso(dt)[1:] # ignoring the litter layer concentration
    theta = sli.theta(dt)
    theta_r = sli.thetar(dt)
    theta_sat = sli.thetasat(dt)
    tortuosity = sli.tortuosity(dt)
    T_soil = sli.T_soil0(dt)
    dT = sli.deltaT(dt)
    rH = sli.R_humidity(dt)
    psi = sli.matric_pot(dt)

    id = 0
    upper_boundary = 0
    for lower_boundary, ciso, th, tr, tsat, tor, T, d_T, r_H, PSI in\
            zip(lower_boundaries, init_c_iso, theta, theta_r, theta_sat, tortuosity, T_soil, dT, rH, psi):
        new_layer = iso_storages.iso_soil_layer(ID=id,
                                                upper_boundary=upper_boundary,
                                                lower_boundary=lower_boundary,
                                                conc_iso_liquid={"2H": 1.0, "18O": ciso},
                                                theta=th,
                                                theta_0=tr,
                                                theta_sat=tsat,
                                                tortuosity=tor,
                                                T=T+ Tzero_sli,
                                                dT=d_T,
                                                rH=r_H,
                                                psi=PSI
                                                )
        id += 1
        upper_boundary = lower_boundary
        c.add_layer(new_layer)

    return c


def _atm(sli):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    dt = 0  # initial states
    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = sli.Tatm(dt) + Tzero_sli
    Rh_atm = sli.R_humidity_atm(dt)
    c_iso_atm = sli.civa(dt) / sli.cva(dt)
    wind_speed = sli.wind_speed(dt)

    atm = iso_storages.iso_atmosphere(conc_iso_liquid={"2H":1.0, "18O":1.0},conc_iso_vapor={"2H":1.0, "18O":c_iso_atm},
                                      T=Tatm, Rh_atmosphere=Rh_atm,
                                      Pa_atmosphere=patm,
                                      wind_speed=wind_speed,
                                      hc=10,  # canopy height [m] (e.g. 40)
                                      d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                                      z0m=0.1 * 10,# roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                                      LAI=1.1,  # leaf area index (e.g. 2.0)
                                      extku=1.5)
    return atm


def update_storages(c, sli, dt):

    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # atmospheric variables
    patm = 1  # from SLI_solve
    Tatm = sli.Tatm(dt) + Tzero_sli
    Rh_atm = sli.R_humidity_atm(dt)
    c_iso_atm = sli.civa(dt) / sli.cva(dt)
    wind_speed = sli.wind_speed(dt)
    c.update_atmosphere(c_atm={'2H': 1.0, '18O': c_iso_atm}, T=Tatm, Rh=Rh_atm, Pa=patm, wind_speed=wind_speed)

    # Update layers
    theta = sli.theta(dt)
    T_soil = np.array(sli.T_soil0(dt)) + Tzero_sli
    dT = sli.deltaT(dt)[1:]  # ignoring litter layer
    rH = sli.R_humidity(dt)
    psi = sli.matric_pot(dt)
    c.update_layers(theta=theta, T=T_soil, dT=dT, rH=rH, psi=psi)

    # boundary storage
    h0 = 0.02
    c.update_pond(pond_height=h0)

    c_aquifer = sli.cali(dt)
    c.update_aquifer(c_iso={"2H": 1.0, "18O": c_aquifer})


def update_boundaries(c, sli, dt):

    # Update states and fluxes for

    # c: cell
    # sli: soil_litter_iso class to import input variables
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision

    # precipitation
    q_prec = sli.qprec(dt)
    c_prec = sli.cprec(dt)
    c.update_precipitation(q_prec=q_prec, c_prec={"2H": 1.0, "18O": c_prec})

    #runoff
    qrunoff = sli.qrunoff(dt)
    c.update_runoff(q_runoff=qrunoff)

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
    c.update_vapor_fluxes(vapor_fluxes=qv)

    # transpiration
    ql_trans = sli.qex(dt)
    c.update_transpiration(q_trans=ql_trans)

    # boundary storage
    c.update_connection_to_aquifer()

    return c


ignore = {'ignoredvi': True, 'ignoredli': True, 'ignorealphai': True, 'ignorealphaik': True}  # Testcases: Mathieu and Bariac (1996)

pth = os.getcwd()
sli = Sli.SlI(pth + '/sli_label3/iso_variables')  # imports all the variable files /variables folder: testcase-1, sig=1

solute = '18O'  #len(sli.get_in_soil())-1  #len(sli.get_in_soil())-1   # dt is the nth time step solution from sli

atm = _atm(sli)  # atmosphere
p = iso_project.iso_project()  # create a project
p.new_cell(atmosphere=atm, area=1)  # add new cell

c = p.get_cells()[0]
_layers(c, sli)  # add soil layers to the new cell
c.install_connections()  # flux connections within the layers

pd = iso_storages.iso_pond(pond_height=0.02)
c.add_pond(pd)

#boundary connections
c.add_evaporation()
c.add_transpiration()
c.add_surface_runoff()
c.add_precipitation()

aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 0.0, "18O": 0.0})  # aquifer as boundary isotope storage
c.add_aquifer(aq, c.layers[-1])  # aquifer connected to bottom layer

c2H_delta = [iso_storages.flux_node.concentration_to_delta(c_iso, solute) for c_iso in c.conc_2H]
c18O_delta = [iso_storages.flux_node.concentration_to_delta(c_iso, solute) for c_iso in c.conc_18O]

c_iso = {'2H': [c.conc_2H], '18O': [c.conc_18O]}
c_iso_delta = {'2H': [c2H_delta], '18O': [c18O_delta]}

#for dt in range(len(sli.get_in_soil())):
dt = len(sli.get_in_soil()) - 1
#print(dt)

delta_t = sli.dt(dt)

# update storage states and boundaries to current time
update_storages(c, sli, dt)
update_boundaries(c, sli, dt)
c.update_c_layers(conc_iso=sli.ciso(dt)[1:], Isotopologue=solute)

current_c_2H = c.conc_2H
current_c_18O = c.conc_18O

dc_18O = p.run(Isotopologue=solute,
               delta_time=delta_t,
               delta_cv=sli.deltacv(dt),
               deltaS_liq=sli.deltaSliq(dt),
               **ignore)

c18O = list(np.array(current_c_18O) + np.array(dc_18O))
cdelta_18O = [iso_storages.flux_node.concentration_to_delta(c_iso, solute) for c_iso in c18O]

#c.update_c_layers(conc_iso=c18O, Isotopologue=solute)

c_iso[solute].append(c18O)
c_iso_delta[solute].append(cdelta_18O)

