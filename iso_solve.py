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

pth = os.getcwd()
sli = Sli.SlI(pth + '/sli_label3/iso_variables')  # imports all the variable files /variables folder: testcase-1, sig=1

ignore = {'ignoredvi': True, 'ignoredli': True, 'ignorealphai': True, 'ignorealphaik': True}
Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli for floating point precision


Ciso = []

##### Soil Littre iso: SUBROUTINE (isotope_vap)

solute, dt = '18O', len(sli.get_in_soil())-1  #len(sli.get_in_soil())-1   # dt is the nth time step solution from sli

delta_t = sli.dt(dt)  # current time step

# atmospheric variables
patm = 1  # from SLI_solve
Tatm = sli.Tatm(dt)
Rh_atm = sli.R_humidity_atm(dt)
c_iso_atm = sli.civa(dt) / sli.cva(dt)
wind_speed = sli.wind_speed(dt)

# Surface variables
T_surface = sli.Ts(dt)
q_ev = sli.qevap(dt)
ql_surface = sli.ql0(dt)
qv_surface = sli.qv0(dt)

ql_trans = sli.qex(dt)
q_prec = sli.qprec(dt)
c_prec = sli.cprec(dt)
qrunoff = sli.qrunoff(dt)

# Soli layer variables
sig = sli.sig(dt)  # implicit/explicit time stepping constant
n_layers = sli.dx(dt)
lower_boundaries = np.cumsum(sli.dx(dt))
init_c_iso = sli.ciso(dt)[1:]  # ignoring litter layer in sli
theta, theta_r, theta_sat = sli.theta(dt), sli.thetar(dt), sli.thetasat(dt)
tortuosity = sli.tortuosity(dt)
r_humidity = sli.R_humidity(dt)
pot = sli.matric_pot(dt)

del_cv = sli.deltacv(dt)
c_v = sli.var_cv(dt)

# Temperartures
dT = sli.deltaT(dt)[1:]  # change in temperature from previous time
Tsoil0 = sli.T_soil0(dt)
# Tsoil = np.array(Tsoil0) + (sig - 1) * np.array(dT)  # SLi: checked

# Saturations
del_sliqice, del_sliq = sli.deltaSliqice(dt), sli.deltaSliq(dt)
sliqice, sliq = sli.Sliqice(dt), sli.Sliq(dt)
S = np.array(sliqice) + np.array(del_sliqice) * (sig - 1)
Sliq = np.array(sliq) + np.array(del_sliq) * (sig - 1)

# soil fluxes
ql = sli.qlsig(dt)[1:]
qv = sli.qvsig(dt)[1:]

# boundary storage
c_aquifer = sli.cali(dt)

####### Isotope Module #######
# Define Atmosphere
atm = iso_storages.iso_atmosphere(T=Tatm + Tzero_sli,
                                  Rh_atmosphere=Rh_atm,
                                  Pa_atmosphere=patm,
                                  wind_speed=wind_speed)
atm.set_conc_iso_vapor(c_iso_atm, solute)

# Define layers
my_layer = []
id = 0
upper_boundary = 0
for lower_boundary, ciso, th, tr, tsat, tor, T, d_T, psi in \
        zip(lower_boundaries, init_c_iso, theta, theta_r, theta_sat, tortuosity, Tsoil0, dT, pot):
    new_layer = iso_storages.iso_soil_layer(ID=id,
                                            upper_boundary=upper_boundary,
                                            lower_boundary=lower_boundary,
                                            conc_iso_liquid={"2H": 1.0, "18O": ciso},
                                            theta=th,
                                            theta_0=tr,
                                            theta_sat=tsat,
                                            tortuosity=tor,
                                            T=T + Tzero_sli,
                                            dT=d_T,
                                            psi=psi,
                                            )
    id += 1
    upper_boundary = lower_boundary
    my_layer.append(new_layer)

# Define aquifer as isotope storage
aq = iso_storages.iso_aquifer(conc_iso_liquid={"2H": 1.0, "18O": c_aquifer})

# Define project
p = iso_project.iso_project()
p.new_cell(atmosphere=atm, area=1)  # add new cell

c = p.get_cells()[0]

# add layers to ta cell
for lr in my_layer:
    c.add_layer(lr)

c.add_evaporation(q_ev=q_ev, Ts=T_surface + Tzero_sli, ql_surface=ql_surface, qv_surface=qv_surface)
c.add_precipitation(q_prec=q_prec, c_prec=c_prec)
c.add_surface_runoff(q_runoff=qrunoff)

c.add_liquid_fluxes(liquid_fluxes=ql)
c.add_vapor_fluxes(vapor_fluxes=qv)
c.install_connections()

# add transpiration to all layers
for l in range(len(my_layer)):
    c.add_transpiration(soil_layer=my_layer[l], ql_transpiration=ql_trans[l])

c.add_aquifer(soil_layer=my_layer[-1], aquifer=aq)

dciso = p.run(Isotopologue=solute, delta_time=delta_t, delta_cv=del_cv, deltaS_liq=del_sliq, **ignore)
c_iso_dt = np.array(init_c_iso) + dciso
Ciso.append(c_iso_dt)



""""

l = 0
for layer in my_layer:

    Tsoil = layer.T + (sig - 1) * layer.dT
    dcv, dsliq = del_cv[l] * (sig - 1), del_sliq[l] * (sig - 1)

    alpha = layer.alpha_i(Isotopologue=solute, T=Tsoil, ignore_alpha_i=ignorealphai)
    d_beta = layer.d_beta(Isotopologue=solute, T=Tsoil, dT=layer.dT, ignore_alpha_i=ignorealphai)
    beta = alpha + sig * d_beta

    cv = layer.cv  # c_v[l]  # instead use: cv = layer.cv
    cvs = cv + dcv

    sat_liq = layer.get_liq_saturation(cv) + dsliq
    s_eff = layer.eff_saturation(cvs, sat_liq, beta)

    if sat_liq < 1:

        del_s_eff = del_sliq[l] + beta * del_cv[l] + cvs * d_beta \
                    - (sat_liq * beta * del_cv[l] + cvs * beta * del_sliq[l] + cvs * sat_liq * d_beta)
    else:
        del_s_eff = del_sliq[l]

    # liquid diffusivity in soil
    Dl0 = lD.dl_i(T=Tsoil, Isotopologue=solute, tortuosity=layer.tortuosity, ignore_dl_i=ignoredl)
    Dl = Dl0 * (min(sat_liq, 1) * ((layer.theta_sat - layer.theta_0) + layer.theta_0)) # according to SLI_solve:

    # vapour diffusivity in soil
    # dv with the soil temperature before time stepping obtained from hydraulic module
    dv = vD.dv_soil_air(T=layer.T, theta=layer.theta, theta_sat=layer.theta_sat, tortuosity=layer.tortuosity)
    Dvi = vD.dv_i(dv=dv, Isotopologue=solute, ignore_dv_i=ignoredvi)
    Dv = Dvi * cvs

    # Surface Fluxes
    if l == 0:

        nksli = eV.nk_sli_solve(thetasat_surface=layer.theta_sat, Sl=sat_liq, ignore_dv_i=ignoredvi)

        # sat_liq, alpha and Dl depend upon Tsoil + (sig - 1), hence passed as argument instead of using as the layer property
        ci_surface = eV.c_iso_liq_surface(Isotopologue=solute, alphai=alpha, dl=Dl, nk=nksli,
                                          ignore_dv_i=ignoredvi,
                                          ignore_dl_i=ignoredl,
                                          ignore_alpha_i=ignorealphai,
                                          ignore_alpha_i_k=ignorealphaik)

        cev_in = eV.get_conc_evapin(Isotopologue=solute, nk=nksli)
        cev_out = eV.get_conc_evapout(Isotopologue=solute, c_iso_surface=ci_surface, nk=nksli,
                                      ignore_alpha_i=ignorealphai)
        qev_in, qev_out = eV.q_evaporation_in, eV.q_evaporation_out



    l += 1

"""
