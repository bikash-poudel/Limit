'''
Created on 15.05.2023

@author: poudel-b
'''
# -*- coding: utf-8 -*-

from math import exp, log
from scipy.sparse import lil_matrix
import numpy as np
from scipy.sparse.linalg import spsolve

# import os

# TODO: diffusion time scale

# Function to run 1D vertical transport

def run_1D_model(atmosphere, layers, Q, BC, time, solutes,
                 ignore_alpha_i=False,
                 ignore_alpha_i_k=False,
                 ignore_dl_i=False,
                 ignore_dv_i=False):

    #n_timesteps = time.time_steps
    ql_up, ql_down, qv_up, qv_down, theta = Q
    theta_t0 = theta[0]  # theta - current timesteps
    theta_t1 = theta[1]  # theta - next timesteps
    Pa = atmosphere.Pa   # atmospheric pressure  Pa = 1 as per SLI

    C = {}
    for solute in solutes:

        A1_ij = []
        A2_ij = []
        A3_ij = []
        Bij = []

        l = 0
        for layer in layers:

            T = layer.T
            porosity = layer.porosity
            tortuosity = layer.tortuosity

            # TODO: Check in SLI where Dv & Dl comes from??? cable_sli_solve.f90::L2291:L2292
            vapor_diffusivity = dv(T=T, Pa=Pa) #* time.total_seconds()  # since diffusion is in m2 / s
            vapor_dv_i = dv_i(T, solute, Pa, ignore_dv_i) #* time.total_seconds()
            n_D = nD(theta=theta_t0[l], theta_0=layer.theta_0)
            eff_vapor_diffusivity = dv_i_eff(porosity=porosity, theta=theta_t0[l], nD=n_D, dv=vapor_diffusivity, dv_i=vapor_dv_i, tortuosity=tortuosity)

            liquid_diffusivity = dl_i(T, solute, tortuosity, theta_t0[l], ignore_dl_i) #* time.total_seconds()
            eff_liquid_diffusivity = dl_i_eff(dl_i=liquid_diffusivity, q_l=ql_up[l])

            alpha = alpha_i(T, solute, ignore_alpha_i)
            beta = beta_i(alpha_i=alpha, density_h2o_vapour=0.0822)    #TODO: check for vapor density

            Dlv_upper = D_lv_eff(eff_liquid_diffusivity, eff_vapor_diffusivity, beta)
            Dlv_current = D_lv_eff(eff_liquid_diffusivity, eff_vapor_diffusivity, beta)
            Dlv_lower = D_lv_eff(eff_liquid_diffusivity, eff_vapor_diffusivity, beta)

            Dlv_up = (Dlv_upper + Dlv_current) / 2
            Dlv_down = (Dlv_current + Dlv_lower) / 2

            Dv_upper = eff_vapor_diffusivity
            Dv_current = eff_vapor_diffusivity
            Dv_lower = eff_vapor_diffusivity

            current_beta_upper = beta
            current_beta = beta
            current_beta_lower = beta

            dz = layer.thickness

            Qlv_up = ql_up[l] + qv_up[l] * (current_beta_upper + current_beta) / 2 \
                     - (Dv_upper + Dv_current) / 2 * (current_beta - current_beta_upper) / dz

            Qlv_down = ql_down[l] + qv_down[l] * (current_beta_lower + current_beta) / 2 \
                       - (Dv_lower + Dv_current) / 2 * (current_beta_lower - current_beta) / dz

            theta_eff_t0 = theta_t0[l] + (layer.porosity - theta_t0[l]) * current_beta
            theta_eff_t1 = theta_t1[l] + (layer.porosity - theta_t1[l]) * current_beta

            delta_z = layer.thickness
            delta_t = time.dt * time.total_seconds()

            dzdt = delta_z / delta_t
            current_co = layer.c_solutes[solute]

            if layer == layers[0]:  # upper boundary
                if BC.upper_boundary_type == 'dirichlet':

                    A2 = 1
                    A2_ij.append(A2)

                    A3 = 0
                    A3_ij.append(A3)

                    B = BC.upper_boundary_content  # conc for n days
                    Bij.append(B)

                if BC.upper_boundary_type == 'neuman':

                    A2 = dzdt * theta_eff_t1 + Dlv_down / dz + Qlv_down / 2
                    A2_ij.append(A2)

                    A3 = - Dlv_down / dz + Qlv_down / 2
                    A3_ij.append(A3)

                    B = dzdt * theta_eff_t0 * current_co + BC.upper_boundary_content   # qi:neuman isotope flux
                    Bij.append(B)

                if BC.upper_boundary_type == 'atmosphere':

                    A2 = dzdt * theta_eff_t1 + Dlv_down / dz + Qlv_down / 2
                    A2_ij.append(A2)

                    A3 = - Dlv_down / dz + Qlv_down / 2
                    A3_ij.append(A3)

                    nk = nK_MathieuBariac(theta_surface=theta_t0[0],
                                          theta_surface_sat=layer.theta_sat,
                                          theta_surface_0=layer.theta_0)

                    alpha_ik = alpha_i_k(dv=vapor_diffusivity,
                                         dv_i=vapor_dv_i,
                                         nK_MathieuBariac=nk,
                                         formulation="MathieuBariac",
                                         ignore_alpha_i_k=ignore_alpha_i_k)

                    Rl_soil = concentration_to_ratio(c_solutes_i=current_co,
                                                     solute_i=solute)
                    Rv_atm = concentration_to_ratio(c_solutes_i=atmosphere.c_atmosphere[solute],
                                                    solute_i=solute)
                    #hs = 0.2  # soil surface relative humidity (-)   TODO:  check the values
                    rh_normalized = rh_atm_normalized(ha=atmosphere.Rh_atmosphere,
                                                      T_atmosphere=atmosphere.T_atmosphere,
                                                      T_soil=layer.T)

                    #rH_soil = theta_t0[0] / layer.porosity  # wetnes(humidity)s of the layer
                    rH_soil = 0.20

                    qi_evap = q_i_Evaporation_Braud_old(E_pot=BC.upper_boundary_content,
                                                       alpha_i=alpha,
                                                       alpha_i_k=alpha_ik,
                                                       Rl_i_Soil=Rl_soil,
                                                       Rv_i_Atmosphere=Rv_atm,
                                                       rH_Soil=rH_soil,
                                                       rH_normalized_Atmosphere=rh_normalized,
                                                       solute_i=solute)
                    """
                    qi_evap_SLI = q_i_Evaporation(qev=-BC.upper_boundary_content,
                                                  solute=solute,
                                                  cv_a=atmosphere.cv_a(Rh=atmosphere.Rh_atmosphere,
                                                                       T=atmosphere.T_atmosphere),
                                                  c_atmosphere=atmosphere.c_atmosphere,
                                                  c_soil=layer.c_solutes,
                                                  ram=ram(wind_speed=atmosphere.wind_speed,
                                                                          hc=40,
                                                                          d0=0.7*40,
                                                                          z0m=0.1*40),
                                                  rs=rs(wind_speed=atmosphere.wind_speed,
                                                                        extku=0.0,
                                                                        LAI=0.5),
                                                  alpha_i=alpha,
                                                  alpha_i_k=alpha_ik)
                    """

                    B = dzdt * theta_eff_t0 * current_co + qi_evap
                    Bij.append(B)

            if layer != layers[0] and layer != layers[-1]:  # intermediate layers

                A1 = - Dlv_up / dz - Qlv_up / 2
                A1_ij.append(A1)

                A2 = dzdt * theta_eff_t1 + Dlv_down / dz + Qlv_down / 2 \
                     + Dlv_up / dz - Qlv_up / 2
                A2_ij.append(A2)

                A3 = - Dlv_down / dz + Qlv_down / 2
                A3_ij.append(A3)

                B = dzdt * theta_eff_t0 * current_co
                Bij.append(B)

            if layer == layers[-1]:  # lower boundary

                if BC.lower_boundary_type == 'dirichlet':

                    A1 = 0
                    A1_ij.append(A1)

                    A2 = 1
                    A2_ij.append(A2)

                    B = BC.lower_boundary_content
                    Bij.append(B)
                if BC.lower_boundary_type == 'neuman':

                    A1 = - Dlv_up / dz - Qlv_up / 2
                    A1_ij.append(A1)

                    A2 = 2 * dzdt * theta_eff_t1 + Dlv_up / dz - Qlv_up / 2
                    A2_ij.append(A2)

                    B = 2 * dzdt * theta_eff_t0 * current_co - BC.lower_boundary_content  # qi = neuman flux isotope
                    Bij.append(B)

            l += 1

        n_layers = len(layers)
        A = lil_matrix((n_layers, n_layers))
        A.setdiag(A1_ij, k=-1)     # TODO: Check for the diagonal-1 matrix starts from - A[1,0]
        A.setdiag(A2_ij, k=0)
        A.setdiag(A3_ij, k=1)

        B = np.asarray(Bij)
        A = A.tocsr()
        Ciso = spsolve(A, B).tolist()
        #C = np.array(Ci_layers).transpose()

        C[solute] = Ciso

    return C


def dv_i_eff(porosity, theta, nD, dv, dv_i, tortuosity):
    """
    Calculates the total (effective) vapour diffusivity for isotope i (m**2/s) in the soil air spaces

    Description
    ===========

    porosity= soil porosity (m**3/m**3 )
    theta= liquid phase (m**3/m**3 )
    tortuosity = tortuosity (m/m)
    dv = vapour diffusivity of water in air (m**2/s)
    dv_i = vapour diffusivity of isotope i in air (m**2/s)
    nD = exponent that describes the isotopic fractioning due to diffusion (-)


    @param porosity: soil porosity (m**3/m**3 )
    @type porosity: Float

    @param theta: liquid phase (m**3/m**3 )
    @type theta: Float

    @param tortuosity: tortuosity (m/m)
    @type tortuosity: Float

    @param dv: vapour diffusivity of water in air (m**2/s)
    @type dv: Float

    @param dv_i: vapour diffusivity of isotope i in air (m**2/s)
    @type dv_i: Float

    @param dv_i: exponent that describes the isotopic fractioning due to diffusion (-)
    @type dv_i: Float
    """
    try:
        # SLI: cable_sli_solve.f90::L2398-2417
        dv_i_eff = (porosity - theta) * tortuosity * dv * (dv_i / dv) ** nD

        return dv_i_eff

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dv_i(T, solute,  Pa=10 ** 5, ignore_dv_i=False):
    """
    Calculates the vapour diffusivity of a solute in free air (m**2/s)

    Description
    ===========

    @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
    @type T: Float

    @param Pa: air pressure in the atmosphere (1bar = 1 atm = 10**5 Pa) [Pa]
    @type Pa: Float

    @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
    @type solute: Sting

    @param ignore_dv_i:  If set to "True" vapour diffusivity of a solute in free air will be set to dv_i = dv --> the same as normal water
    @type ignore_dv_i: Boolean

    SLI: sli_init.f90 Lines: 2366 - L2367 <--checked and own implementation is okay
    """
    try:
        # SLI: cable_sli_solve.f90::L2366 - L2367
        if "18O" == solute:
            b_i = 1.0251  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        elif "2H" == solute:
            b_i = 1.0285  # HDO diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        else:
            raise NotImplementedError

        d_v = dv(T=T, Pa=Pa)

        if ignore_dv_i == True:
            return d_v
        else:
            # SLI: cable_sli_solve.f90::L2373 - L2374
            dv_i = d_v / b_i  # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff

            return dv_i

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dv(T=283.15, Pa=10 ** 5):
    """
    Calculates the vapour diffusivity of water in free air (m**2/s)

    Description
    ===========

    @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
    @type T: Float

    @param Pa: air pressure in the atmosphere (1bar = 1 atm = 10**5 Pa) [Pa]
    @type Pa: Float

    TODO: SLI uses 1 for Pa and not 10.000 why????

    """
    try:
        dva = 2.17e-5  # vapour diffusivity of water in air at 0 degC [m2/s]
        dv = dva * 1e5 / Pa * (T / 273.16) ** 1.88
        # SLI: Dva       = 2.17e-5   ! vapour diffusivity of water in air at 0 degC [m2/s]
        # SLI: Dvs  = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffuxivity of water in air (m2s-1)

        return dv

    except ValueError as err:
        print(err)
        raise NotImplementedError


def nD(theta, theta_0=0.05):
    """
    Calculates the exponent that describes the isotopic fractioning due to diffusion (-) Melayah et al. (1996)


    Description
    ===========


    @param theta: liquid phase (m**3/m**3 )
    @type theta: Float

    @param theta_0: defined residual water content (water remaining at high suctions) (m**3/m**3 ) --> should be 0
    @type theta_0: Float


    """
    try:
        nD = 0.67 + 0.33 * exp(1 - theta / theta_0)

        return nD

    except ZeroDivisionError as err:
        return 0.67


def dl_i_eff(dl_i, q_l=None, hydrodynamic_dispersivity=0.0):
    """
    Calculates the total (effective) liquid diffusivity for isotope i in the soil liquide phase (m**2/s)

    Description
    ===========

    Under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones and are frequently set to zero (Auriault and Adler, 1995)

    @param theta: liquid phase (m**3/m**3 )
    @type theta: Float

    @param tortuosity: tortuosity (m/m)
    @type tortuosity: Float

    @param q_l: liqiud water flux (m/s)
    @type q_l: Float

    @param hydrodynamic_dispersivity: per default set to 0 (m) because under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones
    @type hydrodynamic_dispersivity: Float
    """
    try:
        if q_l == None or hydrodynamic_dispersivity == 0.0:
            dli_eff = dl_i
        elif q_l != None and hydrodynamic_dispersivity != 0.0:
            dli_eff = dl_i + hydrodynamic_dispersivity * abs(q_l)

        return dli_eff

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dl_i(T, solute, tortuosity, theta, ignore_dl_i=False):
    """
    Calculates the liquid diffusivity for isotope i in the soil liquide phase (m**2/s)

    Description
    ===========

    @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
    @type T: Float

    @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
    @type solute: Sting

    @param theta: liquid phase (m**3/m**3 )
    @type theta: Float

    @param tortuosity: tortuosity (m/m)
    @type tortuosity: Float

    @param ignore_dl_i: If set to "True" liquid diffusivity for isotope will be set to 0 = no liquid diffusivity
    @type ignore_dl_i: Boolean

    """
    try:
        if ignore_dl_i == False:
            # SLI: cable_sli_solve.f90::L2381-2390
            dl_i = dl_i_self_diffusivity(T, solute, formulation="Cuntz") * theta * tortuosity
        else:
            dl_i = 0.0

        return dl_i

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dl_i_self_diffusivity(T, solute, formulation="Cuntz"):
    """
    Calculates the liquide self diffusivity of an given isotope in pure water (m**2/s) either after Melayah et al.  (1996) or Cuntz et al. (2007)

    Description
    ===========

    @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
    @type T: Float

    @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
    @type solute: Sting

    @param formulation: Identifier for the formulation to be used to calculated dl_i_self_diffusivity (currently supported 'Melayah'and 'Cuntz')
    @type formulation: Sting

    """
    try:

        if formulation == "Melayah":
            if "18O" == solute:
                a_i = 0.9669
            elif "2H" == solute:
                a_i = 0.9833
            else:
                raise NotImplementedError

            dl_i = a_i * 10 ** -9 * exp(-(535400 / T ** 2) + (1393.3 / T) + 2.1876)

        # SLI: cable_sli_solve.f90::L2381-2390
        elif formulation == "Cuntz":
            if "18O" == solute:
                a_i = 1 / 1.026
            elif "2H" == solute:
                a_i = 1 / 1.013
            else:
                raise NotImplementedError

            dl_i = a_i * 1.0e-7 * exp(-577 / (T - 145))  # TODO: reviewed calculaton from SLI: L2275-L2280

        else:
            raise NotImplementedError

        return dl_i

    except ValueError as err:
        print(err)
        raise NotImplementedError


def D_lv_eff(liquid_diffusivity, vapor_diffusivity, beta):
    return liquid_diffusivity + vapor_diffusivity * beta


def alpha_i(T, solute, ignore_alpha_i=False):
    """
    Calculates the liquid-vapour isotopic frctionation factor at equilibrium given by Majoube (1971) as a function of temperature T (K)

    Description
    ===========

    @param T:  liquid-vapour isotopic frctionation factor at equilibrium given by Majoube (1971) as a function of temperature T (K) at soil depth
    @type T: Float

    @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
    @type solute: Sting

    @param ignore_alpha_i:  If set to "True" isotopic frctionation factor will be set to 1 = no fractionation
    @type ignore_alpha_i: Boolean

    Control status: Checked on 16.05.2013 --> Results are the same as for SLI
    """
    try:
        # SLI: cable_sli_solve.f90::L2275 - L2305
        if ignore_alpha_i == False:
            if "18O" == solute:
                alpha_i = exp(-(24844 / T ** 2 + (-76.248) / T + 0.052612))
            elif "2H" == solute:
                alpha_i = exp(-(1137 / T ** 2 + (-0.4156) / T - 0.0020667))
            else:
                raise NotImplementedError
            return alpha_i
        else:
            return 1.0

    except ValueError as err:
        print(err)
        raise NotImplementedError


def alpha_i_k(dv=None, dv_i=None, nK_MathieuBariac=None, nK_Brutsaert=None, rbh=None, ram=None,
              formulation="Melayah", ignore_alpha_i_k=False):
    """
    Calculates isotopic kinetic fractionation factor under non-saturated and non-isothermal conditions at the soil surface

    Description
    ===========

    Follogwing formulations are available:

    BarnesAllsion: Based on Barnes and Allsion (1983) --> alpha_i_k = dv /dv_i
    MathieuBariac: Based on Mathieu and Bariac (1996) --> alpha_i_k = (dv /dv_i)**nK_MathieuBariac
    Melayah: Based on Melayah et al. (1996) --> alpha_i_k = 1
    Brutsaert: Based on Brutsaert (1982) -->alpha_i_k =  1 + nK_Brutsaert*(dv/dv_i-1)*(r_am/r_a)  <-- with aerodynamic resistance (r_a) = turbulent resistance (r_aT) + molecular resistance (r_am); note nK_Brutsaert != nk MathieuBariac

    @param dv: vapour diffusivity of water in air (m**2/s)
    @type dv: Float

    @param dv_i: vapour diffusivity of isotope i in air (m**2/s)
    @type dv_i: Float

    @param nK_MathieuBariac: exponent that describes the isotopic kinetic fractioning at the soil surface used by the formulation of 'MathieuBariac'
    @type nK_MathieuBariac: float

    @param nK_Brutsaert: exponent that describes the isotopic kinetic fractioning at the soil surface used by the formulation of 'Brutsaert' depending on friction velcity, roughness length and kinematic velocity
    @type nK_Brutsaert: float
    nK_Brutsaert

    @param rbh: rbh =  ra + rs resistance to water vapour transfer beween surface and lowest atmospheric layer (s*m**-1)  used by the formulation of 'MathieuBariac'
    @type rbh: float

    @param ram: molecular resistance between level  hc the soil surface  z0m (s*m**-1)  used by the formulation of 'MathieuBariac'
    @type ram: float

    @param formulation: Identifier for the formulation to be used to calculated alpha_i_k (currently supported 'Melayah', 'BarnesAllsion', 'MathieuBariac' and 'Brutsaert')
    @type formulation: Sting

    @param ignore_alpha_i_k:  If set to "True" kinetic isotopic frctionation factor will be set to 1 = no kinetic fractionation
    @type ignore_alpha_i_k: Boolean
    """
    try:
        if ignore_alpha_i_k == False:
            # SLI: cable_sli_solve.f90::L2275 - L2305
            if formulation == "Melayah":
                alpha_i_k = 1  # SLI: cable_sli_solve.f90::L2379
            elif formulation == "BarnesAllsion":
                if dv != None and dv_i != None:
                    alpha_i_k = dv / dv_i
                else:
                    print("The formualtion of alpha_i_k based on Barnes and Allsion (1983) requires dv  & dv_i")
                    raise NotImplementedError
            elif formulation == "MathieuBariac":
                if nK_MathieuBariac != None and dv != None and dv_i != None:
                    alpha_i_k = (dv / dv_i) ** nK_MathieuBariac  ##SLI: cable_sli_solve.f90::L2377
                else:
                    print(
                        "The formualtion of alpha_i_k based on Mathieu and Bariac (1996) requires nK_MathieuBariac, dv  & dv_i")
                    raise NotImplementedError
            elif formulation == "Brutsaert":
                if nK != None and dv != None and dv_i != None and ram != None and rbh != None:
                    alpha_i_k = 1 + nK_Brutsaert * (dv / dv_i - 1) * (ram / rbh)  ##SLI: cable_sli_solve.f90::L2378
                else:
                    print(
                        "The formualtion of alpha_i_k based on Mathieu and Bariac (1996) requires nK, dv, dv_i, r_a & r_am")
                    raise NotImplementedError
            else:
                raise NotImplementedError
        else:
            alpha_i_k = 1.0

        return alpha_i_k

    except ValueError as err:
        print(err)
        raise NotImplementedError


def beta_i(alpha_i, density_h2o_vapour=0.0822, density_h2o_liquide=1000.00):
    """
    Calculates the factor relating liquid and vapour isotope concentration for species i (-).

    Description
    ===========

    @param alpha_i: liquid-vapour isotopic frctionation factor at equilibrium given by Majoube (1971) as a function of temperature T (K) for 2H or 18O
    @type alpha_i: Float

    @param density_h2o_vapour: Density of vapor at a certain temperature and pressure [kg/m3]
    @type density_h2o_vapour: Float

    @param density_h2o_liquide: Density of water at a certain temperature and pressure [kg/m3]
    @type density_h2o_liquide: Float

    TODO: Check formulation -> according to Braud et al. 2005 Eq. 14 beta_i = alpha_i * density_h2o_vapour / density_h2o_liquide in SLI::cable_sli_solve.f90::2308 beta_i = alpha_i

    """
    try:
        # return alpha_i
        return alpha_i * density_h2o_vapour / density_h2o_liquide

    except ValueError as err:
        print(err)
        raise NotImplementedError

    # beta_i = alpha_i according to SLI


def nK_MathieuBariac(theta_surface, theta_surface_sat, theta_surface_0=0.05):
    """
    Calculates the exponent that describes the isotopic kinetic fractioning at the soil surface (nK = 1 --> dry soils ; nK = 0.5 --> saturated soils) based on Mathieu and Bariac (1996)

    Description
    ===========

    @param theta_surface: liquid phase at soil surface (m**3/m**3 )
    @type theta_surface: Float

    @param theta_surface_sat: water content at saturation (water remaining at high suctions) (m**3/m**3 )
    @type theta_surface_sat: Float

    @param theta_surface_0: defined residual water content (water remaining at high suctions) (m**3/m**3 ) --> should be 0
    @type theta_surface_0: Float

    SLI: cable_sli_solve.f90::L2375 <--checked and own implementation is the same
    """
    try:
        # SLI: cable_sli_solve.f90::L2375
        n_a = 0.5
        n_s = 1.0

        nK = ((theta_surface - theta_surface_0) * n_a + (theta_surface_sat - theta_surface) * n_s) / (
                    theta_surface_sat - theta_surface_0)
        return nK

    except ValueError as err:
        return 0.67


def q_i_Evaporation_Braud_old(E_pot,  # Potential Evaporation (kg m**2 s*1)
                              alpha_i,  # isotopic equilibrium fractionation factor for species i (-)
                              alpha_i_k,  # isotopic kinetic fractionation factor for species i (-)
                              Rl_i_Soil,  # soil surface liquid isotopic ratio of isotopic species i (-)
                              Rv_i_Atmosphere,  # atmospheric vapour isotopic ratio of isotopic species i (-)
                              rH_Soil,  # soil surface relative humidity (-)
                              rH_normalized_Atmosphere,
                             # air surface relative humidity normalised at the soil surface temperature (-)
                             solute_i,  # identifier for the isotope i
                             M_w=0.018,  # molar mass of water (kg)
                             M_i={"2H": 0.019, "18O": 0.020}):  # molar mass of isotope species i (kg)
    """
    Upper boundary condition: q_i_Evaporation

    Calculates the isotope extraction from the top soil layer due to evaporation. Based on Braud et al. 2005

    Description
    ===========

    """

    q_i_Evaporation = (E_pot / alpha_i_k) * (M_i[solute_i] / M_w) * (
                (alpha_i * Rl_i_Soil * rH_Soil) - (rH_normalized_Atmosphere * Rv_i_Atmosphere)) / (
                                  1 - rH_normalized_Atmosphere)

    return q_i_Evaporation


def q_i_Evaporation(qev, solute, cv_a, c_atmosphere, c_soil, ram, rs, alpha_i, alpha_i_k):
    """
    Upper boundary condition: q_i_Evaporation

    Returns the isotope extraction from the top soil layer due to evaporation at the current time step in kg.

    Description
    ===========
    @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
    @type solute: Sting

    @param q_Evaporation: Flux of precipitation in m**3/d
    @type q_Evaporation: float

    @param cv_a: concentration of atmospheric water vapour (m3 H2O(l)/ m3 (air))    Flux of precipitation in m**3/d
    @type cv_a: float

    @param c_atmosphere: Concentration of isotope species i in the atmospheric vapor in kg/m**3 (currently supported "2H" and/or "18O")
    @type c_atmosphere: dictionary containing two entries {"2H" : 1.0,"18O" : 1.0}

    @param c_soil: Concentration of isotope species i in the soil water liquide phase in kg/m**3 (currently supported "2H" and/or "18O")
    @type c_soil: dictionary containing two entries {"2H" : 1.0,"18O" : 1.0}

    @param ram: Aerodyamic resistance (from z0m to hc)
    @type ram: float

    @param rs: Soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).
    @type rs: float

    @param alpha_i: liquid-vapour isotopic fractionation factor at equilibrium given by Majoube (1971) as a function of temperature T (K) in the layer at the interface (litter, top soil layer or pond)
    @type alpha_i: float

    @param alpha_i_k: isotopic kinetic fractionation factor under non-saturated and non-isothermal conditions at the soil surface
    @type alpha_i_k: float

    @param delta_time: Time step
    @type delta_time: datetime.timedelta

    @return: Returns the isotope extraction from the top soil layer due to evaporation at the current time step in kg.

    TODO Check each value
    """
    #TODO: check ci_va is multiplied with cv_a in SLI
    ci_va = c_atmosphere[solute] * cv_a  # concentration of isotope species i in the atmosphere  TODO: check the concentration neecds o be again multplied with vapor concentration
    ci_soil = c_soil[solute]  # concentration of isotope species i in the soil

    #ql_Evaporation = self.ql_evaporation(layer, delta_time)  # get liquid evaporation flux in kg for the current time step

    r_b_h = rbh(ram, rs)  # resistance to water vapour transfer beween surface and lowest atmospheric layer

    # print ("this needs to account for the area!")
    cv_s = cv_a + qev * (ram + r_b_h)  # cv_s = concentration of water vapour at soil/air interface (m3 (H2O liq)/ m3 (air))

    if cv_s < 0.0:
        print("cv_s this can not be come less then zero")
        raise NotImplementedError

    # cv_a Flaeche??Volumen??

    qevapin = ci_va / (ram + r_b_h)  # amount of water (vapour) entering the soil
    qevapout = cv_s / (ram + r_b_h)  # amount of water leaving the soil due to evaporation

    cevapin = ci_va / cv_a * alpha_i_k  # concentration of isotope species i entering the top soil layer
    cevapout = alpha_i_k * alpha_i * ci_soil  # concentration of isotope species i leaving the top soil layer

    q_i_Evaporation = qevapout * cevapout - qevapin * cevapin

    return q_i_Evaporation


def rbh(ram, rs):
    """
    Returns the boundary layer resistance. Resistance to water vapour transfer beween surface and lowest atmospheric layer.

    Description
    ===========
    @param ram: Aerodyamic resistance (from z0m to hc)
    @type ram: float

    @param rs: Soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).
    @type rs: float

    @return: Returns the boundary layer resistance. Resistance to water vapour transfer beween surface and lowest atmospheric layer.

    Control status: Checked on 15.05.2013 --> Results are the same as for SLI
    """
    rbh = ram + rs  # resistance to water vapour transfer beween surface and lowest atmospheric layer

    return rbh


def rs(wind_speed, extku, LAI):
    """
    Returns the soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).

    Description
    ===========
    @param wind_speed: wind speed at top of canopy [m/s]
    @type wind_speed: float

    @param LAI: leaf area index (e.g. 2.0)
    @type LAI: float

    @param extku: extinction coeff't for windspeed (e.g. 0.5)
    @type extku: float

    @return:  Returns the soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).

    Control status: Checked on 15.05.2013 --> Results are the same as for SLI

    """
    us = wind_speed * exp(-extku * LAI)  # us = windspeed below canopy [m/s]
    rs = 1. / (0.004 + 0.012 * us)

    return rs


def ram(wind_speed, hc, d0, z0m):
    """
    Returns the aerodyamic resistance (from z0m to hc).

    Description
    ===========
    @param wind_speed: wind speed at top of canopy [m/s]
    @type wind_speed: float

    @param hc: canopy height [m] (e.g. 40)
    @type hc: float

    @param d0: displacement height (e.g. 0.7)
    @type d0: float

    @param z0m: roughness height for momentum (e.g. 0.1)
    @type z0m: float

    @return: Returns the aerodyamic resistance (from z0m to hc).

    Control status: Checked on 15.05.2013 --> Results are the same as for SLI
    """
    ram = log((hc - d0) / z0m) ** 2 / 0.41 ** 2 / wind_speed

    return ram


def update_c_i(c_i_layers, solute, layers):
    """
    @param c_i_layers:  List of new isotope concentrations of "solute". Negative liquid flux (m3/d) leaving the cell or positive entering the cell from/to the upper cell
    @type c_i_layers: List

    @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
    @type solute: Sting
    """
    try:
        n_layers = len(layers)
        assert (len(c_i_layers) == n_layers), "len of list theta_layers must be of length(layers)"

        updated_layer = []
        for i_layer, i_c_i in zip(layers, c_i_layers):
            i_layer.c_solutes[solute] = i_c_i
            updated_layer.append(i_layer)

    except AssertionError as err:
        print(err)
        raise NotImplementedError

    return updated_layer


def soil_surface_Rh(Rh_atmosphere, T_soil, T_atmosphere):
    """
    Returns the relative humidity at the soil surface

    Rh_atmosphere = relative humidity in the atmosphere
    T_soil = Temperature in Kelvin at the first soil layer ???
    pv_sat_soil = saturated water vapor volumetric mass of the soil (kg/m**3)
    T_atmosphere = Temperature in Kelvin in the atmosphere
    pv_sat_atmosphere = saturated water vapor volumetric mass of the atmosphere (kg/m**3)
    """

    soil_surface_Rh = Rh_atmosphere * pv_sat(T=T_atmosphere) / pv_sat(T=T_soil)

    return soil_surface_Rh


def delta_to_concentration(delta_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}, M_w=0.018,
                           M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    R_i = delta_to_ratio(delta_i, solute_i, R_ref)  # isotopic ratio of isotopic species i in sample
    c_solutes = R_i * density_H2O * M_i[solute_i] / M_w
    return c_solutes


def delta_to_ratio(delta_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    R_i = ((delta_i) / 1000 + 1) * R_ref[solute_i]  # isotopic ratio of isotopic species i in sample
    return R_i


def ratio_to_delta(R_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    delta_i = (R_i / R_ref[solute_i] - 1) * 1000
    return delta_i


def concentration_to_delta(c_solute_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}, M_w=0.018,
                           M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
    """
    Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

    delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    R_i = concentration_to_ratio(c_solute_i, solute_i, M_w, M_i,
                                      density_H2O)  # isotopic ratio of isotopic species i

    delta_i = ratio_to_delta(R_i, solute_i, R_ref)

    return delta_i


def concentration_to_ratio(c_solutes_i, solute_i, M_w=0.018, M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
    """
    Converts the given concentration of isotope species i  (kg/m**3) into the ratio of isotope species i

    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """

    R_i = c_solutes_i / (density_H2O * M_i[solute_i] / M_w)

    return R_i


def rh_atm_normalized(ha, T_soil, T_atmosphere):

    rh = ha * pv_sat(T_atmosphere) / pv_sat(T_soil)

    return rh


def pv_sat(T):
    """
    Returns the saturated water vapor volumetric mass (kg/m**3)

    T_atmosphere = Temperature in Kelvin

    pv_sat = saturated water vapor volumetric mass (kg/m**3)

    TODO: Check for which temperature range it is valid

    Control status: Checked on 15.05.2013 --> Results are the same as for SLI
    """

    # SLI: sli_main.f90::L205
    pv_sat = 0.002166 * p_sat(T) / (T + 273.16)  # kg/m3

    return pv_sat


def p_sat(T):
    """
    Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch and valid over the range 0 to 60 Degree celsius

    T_atmosphere = Temperature in Kelvin

    p_sat = saturated water vapor volumetric mass (kg/m**3)

    SLI: cable_sli_utils.f90 Lines: 2513-2525 <--checked and own implementation is okay for temp. above 273.15 Kelvin
    """

    # SLI: cable_sli_utils.f90::L2513-2525

    p_sat = 610.6 * exp(17.27 * (T - 273.15) / ((T - 273.15) + 237.3))  # Saturation vapour pressure, ps, in pascal

    if T < 273.16:  # below 0 Degree celsius
        p_sat = -4.86 + 0.855 * p_sat + 0.000244 * p_sat ** 2  # ps for ice

    return p_sat


