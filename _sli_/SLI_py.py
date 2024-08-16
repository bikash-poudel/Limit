import numpy as np
from math import exp, log

from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

import Sli

import os

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
            if "2H" == solute:
                alpha_i = exp(-(24844 / T ** 2 + (-76.248) / T + 0.052612))
            elif "18O" == solute:
                alpha_i = exp(-(1137 / T ** 2 + (-0.4156) / T - 0.0020667))
            else:
                raise NotImplementedError
            return alpha_i
        else:
            return 1.0

    except ValueError as err:
        print(err)
        raise NotImplementedError


def alpha_plus(T, solute, ignore_alpha_i=False):
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
        if not ignore_alpha_i:
            if "2H" == solute:
                alpha_i = 1 / exp(24844 / T ** 2 + (-76.248) / T + 0.052612)
            elif "18O" == solute:
                alpha_i = 1 / exp(1137 / T ** 2 + (-0.4156) / T - 0.0020667)
            else:
                raise NotImplementedError
            return alpha_i
        else:
            return 1.0

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dalphplus_dT(T, solute, ignore_alpha_i=False):
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
        if not ignore_alpha_i:
            if "2H" == solute:
                dalphaplusdT = (2 * 24844 / T ** 3 + (-76.248) / T ** 2) / exp(
                    24844 / T ** 2 + (-76.248) / T + 0.052612)
            elif "18O" == solute:
                dalphaplusdT = (2 * 1137 / T ** 3 + (-0.4156) / T ** 2) / exp(1137 / T ** 2 + (-0.4156) / T - 0.0020667)
            else:
                raise NotImplementedError
            return dalphaplusdT
        else:
            return 0

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dv_i(T, solute, Pa=10 ** 5, ignore_dv_i=False):
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
        if "2H" == solute:
            b_i = 1.0251  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        elif "18O" == solute:
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
        Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli
        dva = 2.1699999706470408E-005  # precision value from sli, 2.17e-5  # vapour diffusivity of water in air at 0 degC [m2/s]
        dv = dva * 1e5 / Pa * (T / Tzero_sli) ** 1.88
        # SLI: Dva       = 2.17e-5   ! vapour diffusivity of water in air at 0 degC [m2/s]
        # SLI: Dvs  = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffuxivity of water in air (m2s-1)

        return dv

    except ValueError as err:
        print(err)
        raise NotImplementedError


def n_k(thetasat, S):
    ## Different formulation in SLI: Appendix : B.7
    
    # Testcase 7 or 8, nk=1
    nk = ((thetasat * min(S, 1) - 0) * 0.5 + (thetasat * (1 - min(S, 1)))) \
            / (thetasat - 0)

    return nk
    

def alpha_k(dvs, divs, nk, ignore_alpha_k=True):
    
    if ignore_alpha_k:
        alphak = 1
    else:
        
        alphak = 1 / ((dvs / divs) ** nk)

    return alphak


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


def dl_i(T, solute, tortuosity, ignore_dl_i=False):
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
            dl_i = dl_i_self_diffusivity(T, solute, formulation="Cuntz") * tortuosity
        else:
            dl_i = 0.0

        return dl_i

    except ValueError as err:
        print(err)
        raise NotImplementedError


def dv_sli(d_v, solute, ignore_dv_i=False):
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
        if "2H" == solute:
            b_i = 1.0251  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        elif "18O" == solute:
            b_i = 1.0285  # HDO diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        else:
            raise NotImplementedError

        if ignore_dv_i:
            return d_v
        else:
            # SLI: cable_sli_solve.f90::L2373 - L2374
            dv_i = d_v / b_i  # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff

            return dv_i

    except ValueError as err:
        print(err)
        raise NotImplementedError


def alphak_vdiff(solute, ignore_dv_i=False):
    try:
        # SLI: cable_sli_solve.f90::L2366 - L2367
        if "2H" == solute:
            b_i = 1.0251  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        elif "18O" == solute:
            b_i = 1.0285  # HDO diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
        else:
            raise NotImplementedError

        if ignore_dv_i:
            return 1
        else:
            # SLI: cable_sli_solve.f90::L2373 - L2374
            alphak_diff = 1 / b_i  # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff

            return alphak_diff

    except ValueError as err:
        print(err)
        raise NotImplementedError


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


def concentration_to_ratio(c_solutes_i, solute_i, M_w=0.018, M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
    """
    Converts the given concentration of isotope species i  (kg/m**3) into the ratio of isotope species i

    solute_i = sting identifying the isotope species currently supported '2H' or '18O'
    M_w = molar mass of water (kg)  0.018 kg for H_2 O
    M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
    mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

    SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
    """
    #M_w = 1.7999999225139618E-002  # Value from SLI due to floating point error
    R_i = c_solutes_i / (density_H2O * M_i[solute_i] / M_w)

    return R_i


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


def solve_iso(sli, dt, ignore):

    
    ###### Soil Littre iso: SUBROUTINE (isotope_vap)

    ignorealphai, ignorealphaik, ignoredl, ignoredvi = ignore
    
    if sli.isotopologue(dt) == '1':        
        solute = '2H'     
    elif sli.isotopologue(dt) == '2':
        solute = '18O'
    else:
        raise ValueError
    

    sig = sli.sig(dt)
    Tzero_sli = 273.16000366210938  # [k] 0 celcius in kelvin, value taken from sli
    deltaz0 = 0.5 * sli.dx(dt)[0]
    patm = 1
    Tsoil = np.array(sli.T_soil0(dt)) + (sig - 1) * np.array(sli.deltaT(dt)[1:])  # SLi: checked

    # SlI has slightly different formulation in Appendix B.7    
        
    alphaplus_s = alpha_plus(sli.Ts(dt) + Tzero_sli, solute, ignorealphai)  # soil surface ,  SLi: checked
    cvsig_all = np.array(sli.var_cv(dt)) + np.array(sli.deltacv(dt)) * (sig - 1)  # SLI: checked

    Beta = []
    Sicesig, Seff, deltaSeff = [], [], []

    Dv, Dl, wl, wv = [], [], [], []

    dcevapoutdciso = 0
    qevapin, qevapout, cevapin, cevapout = 0, 0, 0, 0

    cql, cqv = [], []
    dcqldca, dcqldcb, dcqvdca, dcqvdcb = [], [], [], []
    betaqv, dbetaqv = [], []

    for l in range(sli.n(dt)):  # soil layers
       
            
        alphaplpus = alpha_plus(Tsoil[l] + Tzero_sli, solute, ignore_alpha_i=ignorealphai)  # SLi: checked
        dalphplusdT = dalphplus_dT(Tsoil[l] + Tzero_sli, solute, ignore_alpha_i=ignorealphai)  # SLi: checked
        
        _beta = alphaplpus
        deltabeta = dalphplusdT * sli.deltaT(dt)[l + 1]  # ignoring litter layer, Sli: checked

        beta = _beta + sig * deltabeta  # Sli: checked

        ## ice , SLI: checked
        S = sli.Sliqice(dt)[l] + sli.deltaSliqice(dt)[l] * (sig - 1)
        Sliqsig = sli.Sliq(dt)[l] + sli.deltaSliq(dt)[l] * (sig - 1)  # Sliqsig including pond
        _Sicesig = sli.Sice(dt)[l] + sli.deltaSice(dt)[l] * (sig - 1)
        deltathetaice = sli.deltaSice(dt)[l] * (sli.thetasat(dt)[l] - sli.thetar(dt)[l])

        # h0 = sli.h0new(dt) + sli.dh0(dt) * (sig - 1)
        cvsig = sli.var_cv(dt)[l] + sli.deltacv(dt)[l] * (sig - 1)
        # check formulation for _Seff
        if S < 1:  # SLI: checked
            _Seff = Sliqsig + cvsig * beta - cvsig * S * beta + (sli.thetar(dt)[l] / sli.thetasat(dt)[l])
            _deltaSeff = sli.deltaSliq(dt)[l] + beta * sli.deltacv(dt)[l] + cvsig * deltabeta \
                         - (S * beta * sli.deltacv(dt)[l] + cvsig * beta * sli.deltaSliqice(dt)[l]
                            + cvsig * S * deltabeta)
        else:
            _Seff = Sliqsig + (sli.thetar(dt)[l] / sli.thetasat(dt)[l])
            _deltaSeff = sli.deltaSliq(dt)[l]

        ## liquid diffusivity in soil      
        D_l = dl_i(Tsoil[l] + Tzero_sli, solute, sli.tortuosity(dt)[l], ignore_dl_i=ignoredl)
        _Dl = D_l * (min(S, 1) * (sli.thetasat(dt)[l] - sli.thetar(dt)[l]) + sli.thetar(dt)[
            l])  ## Recheck the value of DL from SLI

        ## vapour diffusivity in soil       
        D_v = dv_sli(sli.var_dv(dt)[l], solute, ignore_dv_i=ignoredvi)
        _Dv = D_v * cvsig  ## SLI: checked

        if l < sli.n(dt) - 1:  ## SLI: checked

            if abs((cvsig_all[l] - cvsig_all[l + 1])) > 1.e-8:

                Dv_eff = sli.qvsig(dt)[l + 1] / \
                         (cvsig_all[l] - cvsig_all[l + 1]) * sli.deltaz(dt)[l] * alphak_vdiff(solute,
                                                                                              ignore_dv_i=ignoredvi)

            else:
                Dv_eff = dv_sli(sli.var_dv(dt)[l], solute, ignore_dv_i=ignoredvi)

        if sli.ql0(dt) > 0:
            w1 = 0
        else:
            w1 = 1

        if l == 0:  # 1st layer,  Sli: checked            
                
            ## diffusional fractionation factor in air
            if sli.testcase(dt) == '7' or sli.testcase(dt) == '8':                
                nk = 1
            elif sli.testcase(dt) in ['1', '2', '3', '4', '5', '6']:
                nk = n_k(sli.thetasat(dt)[l], S)
            else:
                raise ValueError
            
            Divs = dv_i(sli.Ts(dt) + Tzero_sli, solute, Pa=patm, ignore_dv_i=ignoredvi)
            Dvs = dv(sli.Ts(dt) + Tzero_sli, Pa=patm)
            alphak = alpha_k(Dvs, Divs, nk, ignore_alpha_k=ignorealphaik)           

            ## Upper boundary condition
            cvs = sli.cva(dt) + sli.qevap(dt) * (sli.ram(dt) + sli.rbh(dt))  # Sli: checked

            if sli.var_dv(dt)[l] != 0:  # Sli: checked
                cv1 = - sli.qv0(dt) * (0.5 * sli.dx(dt)[l]) / sli.var_dv(dt)[l] + cvs
            else:
                cv1 = sli.var_cv(dt)[l]

            # surface concentration of minor isotopologue at surface  #SLI checked
            alphaplus_s = alpha_plus(sli.Ts(dt) + Tzero_sli, solute, ignore_alpha_i=ignorealphai)  # surface
            Dv_eff_0 = sli.var_dv(dt)[0]

            num = alphak_vdiff(solute, ignore_dv_i=ignoredvi) * Dv_eff_0 / (0.5 * sli.dx(dt)[0]) \
                  * cv1 * alphaplpus * sli.ciso(dt)[1] \
                  - sli.ql0(dt) * sli.ciso(dt)[1] * w1 \
                  + sli.civa(dt) * alphak / (sli.ram(dt) + sli.rbh(dt)) \
                  + _Dl * sli.ciso(dt)[1] / (0.5 * sli.dx(dt)[0])

            den = alphak * cvs * alphaplus_s / (sli.ram(dt) + sli.rbh(dt)) \
                  + sli.ql0(dt) * (1 - w1) \
                  + alphaplpus * alphak_vdiff(solute, ignore_dv_i=ignoredvi) * cvs * Dv_eff_0 / (0.5 * sli.dx(dt)[0]) \
                  + _Dl / (0.5 * sli.dx(dt)[0])

            cisos = num / den
            dcevapout_dciso = alphak * alphaplus_s
            dcevapoutdciso = dcevapout_dciso * (alphak_vdiff(solute, ignore_dv_i=ignoredvi) * Dv_eff_0
                                                / (0.5 * sli.dx(dt)[0]) * sli.var_cv(dt)[0] * alphaplpus - sli.ql0(dt)
                                                + _Dl / (0.5 * sli.dx(dt)[0])) / den

            # SlI: checked
            qevapin = sli.cva(dt) / (sli.ram(dt) + sli.rbh(dt))
            qevapout = cvs / (sli.ram(dt) + sli.rbh(dt))
            cevapin = sli.civa(dt) / sli.cva(dt) * alphak
            cevapout = alphak * alphaplus_s * cisos 

        ## all output variables checked with SLI
        Seff.append(_Seff), deltaSeff.append(_deltaSeff)
        Sicesig.append(_Sicesig)

        Beta.append(beta)
        Dv.append(_Dv), Dl.append(_Dl)

        wl.append(sli.dx(dt)[l]), wv.append(sli.dx(dt)[l])

    for l in range(sli.n(dt)):  # soil layers

        ### Concentration of advective fluxes and corresponding partial derivs wrt ciso [formulation = 2]
        if l != sli.n(dt) - 1:  # not including ciso[0] ignoring litter later, SLI: checked

            _cql = 0.5 * (sli.ciso(dt)[l + 1] + sli.ciso(dt)[l + 2])
            _dcqldca = 0.5
            _dcqldcb = 0.5

            _cqv = 0.5 * (sli.ciso(dt)[l + 1] * Beta[l] + sli.ciso(dt)[l + 2] * Beta[l + 1])
            _dcqvdca = 0.5 * Beta[l]
            _dcqvdcb = 0.5 * Beta[l + 1]

        else:

            _cql = sli.ciso(dt)[l + 1]
            _dcqldca = 1
            _dcqldcb = 0

            _cqv = sli.ciso(dt)[l + 1] * Beta[l]
            _dcqvdca = 1 * Beta[l]
            _dcqvdcb = 0

        _betaqv = alphak_vdiff(solute, ignore_dv_i=ignoredvi)
        _dbetaqv = 0

        # all output variables checked with SLI
        cql.append(_cql), dcqldca.append(_dcqldca), dcqldcb.append(_dcqldcb)
        cqv.append(_cqv), dcqvdca.append(_dcqvdca), dcqvdcb.append(_dcqvdcb)

        betaqv.append(_betaqv), dbetaqv.append(_dbetaqv)

    # mean diffusivities
    Dl_mean = (np.array(Dl[0:-1]) * np.array(wl[0:-1]) + np.array(Dl[1:]) * np.array(wl[1:])) / (
                np.array(wl[0:-1]) + np.array(wl[1:]))
    Dv_mean = (np.array(Dv[0:-1]) * np.array(wv[0:-1]) + np.array(Dv[1:]) * np.array(wv[1:])) / (
                np.array(wv[0:-1]) + np.array(wv[1:]))
    Dl_mean = np.append(Dl_mean, 0)
    Dv_mean = np.append(Dv_mean, 0)

    # Coefficient of tridiagonal matrix
    aa, bb, cc, dd = [], [], [], []

    for l in range(sli.n(dt)):

        # # Coefficient of tridiagonal matrix
        if l == 0:  # first layer

            b = - Seff[0] * (sli.thetasat(dt)[0] - sli.thetar(dt)[0]) * sli.dx(dt)[0] / sig / sli.dt(dt) \
                - qevapout * dcevapoutdciso \
                - sli.qlsig(dt)[1] * dcqldca[0] - sli.qvsig(dt)[1] * betaqv[0] * dcqvdca[0]  \
                - sli.qvsig(dt)[1] * dbetaqv[0] * Dv_mean[0] \
                - Dl_mean[0] / sli.deltaz(dt)[0] \
                - Dv_mean[0] / sli.deltaz(dt)[0] * Beta[0] \
                - (sli.qex(dt)[0] + sli.qrunoff(dt))

            c = - sli.qlsig(dt)[1] * dcqldcb[0] - sli.qvsig(dt)[1] * betaqv[0] * dcqvdcb[0] \
                - sli.qvsig(dt)[1] * dbetaqv[0] * dcqvdcb[0] * Dv_mean[0] \
                + Dl_mean[0] / sli.deltaz(dt)[0] \
                + Dv_mean[0] / sli.deltaz(dt)[0] * Beta[1]

            d = (sli.thetasat(dt)[0] - sli.thetar(dt)[0]) * sli.dx(dt)[0] / sig / sli.dt(dt) \
                * (sli.ciso(dt)[1] * deltaSeff[0] + sli.cisoice(dt)[1] * sli.deltaSice(dt)[0]) \
                - sli.qprec(dt) * sli.cprec(dt) / sig - qevapin * cevapin / sig \
                + qevapout * cevapout / sig \
                + sli.qlsig(dt)[1] * cql[0] / sig + sli.qvsig(dt)[1] * betaqv[0] * cqv[0] / sig \
                + sli.qvsig(dt)[1] * dbetaqv[0] * cqv[0] * Dv_mean[0] / sig \
                + Dl_mean[0] / sli.deltaz(dt)[0] / sig * (sli.ciso(dt)[1] - sli.ciso(dt)[2]) \
                + Dv_mean[0] / sli.deltaz(dt)[0] / sig * (sli.ciso(dt)[1] * Beta[0] - sli.ciso(dt)[2] * Beta[1]) \
                + (sli.qex(dt)[0] + sli.qrunoff(dt)) * sli.ciso(dt)[1] / sig

            aa.append(0), bb.append(b), cc.append(c), dd.append(d)

        elif l < sli.n(dt) - 1:  # Intermediate layers

            # SLI: checked
            a = sli.qlsig(dt)[l] * dcqldca[l - 1] + sli.qvsig(dt)[l] * betaqv[l - 1] * dcqvdca[l - 1] \
                + sli.qvsig(dt)[l] * dbetaqv[l - 1] * dcqvdca[l - 1] * Dv_mean[l - 1] \
                + Dl_mean[l - 1] / sli.deltaz(dt)[l - 1] + Dv_mean[l - 1] / sli.deltaz(dt)[l - 1] * Beta[l - 1]

            b = - Seff[l] * (sli.thetasat(dt)[l] - sli.thetar(dt)[l]) * sli.dx(dt)[l] / sig / sli.dt(dt) \
                + sli.qlsig(dt)[l] * dcqldcb[l - 1] + sli.qvsig(dt)[l] * betaqv[l - 1] * dcqvdcb[l - 1] \
                + sli.qvsig(dt)[l] * dbetaqv[l - 1] * dcqvdcb[l - 1] * Dv_mean[l - 1] \
                - sli.qlsig(dt)[l + 1] * dcqldca[l] - sli.qvsig(dt)[l + 1] * betaqv[l] * dcqvdca[l] \
                - sli.qvsig(dt)[l + 1] * dbetaqv[l] * dcqvdca[l] * Dv_mean[l] \
                - (Dl_mean[l - 1] / sli.deltaz(dt)[l - 1] + Dl_mean[l] / sli.deltaz(dt)[l]) \
                - (Dv_mean[l - 1] / sli.deltaz(dt)[l - 1] + Dv_mean[l] / sli.deltaz(dt)[l]) * Beta[l] \
                - sli.qex(dt)[l]

            c = - sli.qlsig(dt)[l + 1] * dcqldcb[l] - sli.qvsig(dt)[l + 1] * betaqv[l] * dcqvdcb[l] \
                - sli.qvsig(dt)[l + 1] * dbetaqv[l] * dcqvdcb[l] * Dv_mean[l] \
                + Dl_mean[l] / sli.deltaz(dt)[l] \
                + Dv_mean[l] / sli.deltaz(dt)[l] * Beta[l + 1]

            d = (sli.thetasat(dt)[l] - sli.thetar(dt)[l]) * sli.dx(dt)[l] / sig / sli.dt(dt) \
                * (sli.ciso(dt)[l + 1] * deltaSeff[l] + sli.cisoice(dt)[l + 1] * sli.deltaSice(dt)[l]) \
                - sli.qlsig(dt)[l] * cql[l - 1] / sig \
                - sli.qvsig(dt)[l] * betaqv[l - 1] * cqv[l - 1] / sig \
                - sli.qvsig(dt)[l] * dbetaqv[l - 1] * cqv[l - 1] * Dv_mean[l - 1] / sig \
                + sli.qlsig(dt)[l + 1] * cql[l] / sig \
                + sli.qvsig(dt)[l + 1] * betaqv[l] * cqv[l] / sig \
                + sli.qvsig(dt)[l + 1] * dbetaqv[l] * cqv[l] * Dv_mean[l] / sig \
                - Dl_mean[l - 1] / sli.deltaz(dt)[l - 1] / sig * (sli.ciso(dt)[l] - sli.ciso(dt)[l + 1]) \
                - Dv_mean[l - 1] / sli.deltaz(dt)[l - 1] / sig * (
                            sli.ciso(dt)[l] * Beta[l - 1] - sli.ciso(dt)[l + 1] * Beta[l]) \
                + Dl_mean[l] / sli.deltaz(dt)[l] / sig * (sli.ciso(dt)[l + 1] - sli.ciso(dt)[l + 2]) \
                + Dv_mean[l] / sli.deltaz(dt)[l] / sig * (
                            sli.ciso(dt)[l + 1] * Beta[l] - sli.ciso(dt)[l + 2] * Beta[l + 1]) \
                + sli.qex(dt)[l] * sli.ciso(dt)[l + 1] / sig

            aa.append(a), bb.append(b), cc.append(c), dd.append(d)

        else:  # last layer

            # SLI: checked
            a = sli.qlsig(dt)[l] * dcqldca[l - 1] + sli.qvsig(dt)[l] * betaqv[l - 1] * dcqvdca[l - 1] \
                + sli.qvsig(dt)[l] * dbetaqv[l - 1] * dcqvdca[l - 1] * Dv_mean[l - 1] \
                + Dl_mean[l - 1] / sli.deltaz(dt)[l - 1] + Dv_mean[l - 1] / sli.deltaz(dt)[l - 1] * Beta[l - 1]

            b = - Seff[l] * (sli.thetasat(dt)[l]  - sli.thetar(dt)[l]) * sli.dx(dt)[l] / sig / sli.dt(dt) \
                + sli.qlsig(dt)[l] * dcqldcb[l - 1] + sli.qvsig(dt)[l] * betaqv[l - 1] * dcqvdcb[l - 1] \
                + sli.qvsig(dt)[l] * dbetaqv[l - 1] * dcqvdcb[l - 1] * Dv_mean[l - 1] \
                - sli.qlsig(dt)[l + 1] * dcqldca[l] \
                - Dl_mean[l - 1] / sli.deltaz(dt)[l - 1] \
                - Dv_mean[l - 1] / sli.deltaz(dt)[l - 1] * Beta[l] \
                - sli.qex(dt)[l]

            d = (sli.thetasat(dt)[l] - sli.thetar(dt)[l]) * sli.dx(dt)[l] / sig / sli.dt(dt) \
                * (sli.ciso(dt)[l + 1] * deltaSeff[l] + sli.cisoice(dt)[l + 1] * sli.deltaSice(dt)[l]) \
                - sli.qlsig(dt)[l] * cql[l - 1] / sig \
                - sli.qvsig(dt)[l] * betaqv[l - 1] * cqv[l - 1] / sig \
                - sli.qvsig(dt)[l] * dbetaqv[l - 1] * cqv[l - 1] * Dv_mean[l - 1] / sig \
                + sli.qlsig(dt)[l + 1] * cql[l] / sig \
                - Dl_mean[l - 1] / sli.deltaz(dt)[l - 1] / sig * (sli.ciso(dt)[l] - sli.ciso(dt)[l + 1]) \
                - Dv_mean[l - 1] / sli.deltaz(dt)[l - 1] / sig * (
                            sli.ciso(dt)[l] * Beta[l - 1] - sli.ciso(dt)[l + 1] * Beta[l]) \
                + sli.qex(dt)[l] * sli.ciso(dt)[l + 1] / sig         
            
            aa.append(a), bb.append(b), cc.append(0), dd.append(d)


    if sli.cali(dt) > 0 or sli.testcase(dt) == '7' or sli.testcase(dt) == '8':
        bb[-1] = bb[-1] + sli.qlsig(dt)[-1] * dcqldca[-1]
        dd[-1] = dd[-1] + sli.qlsig(dt)[-1] * (sli.cali(dt) - cql[-1] / sig)

    n_layers = sli.n(dt)
    A = lil_matrix((n_layers, n_layers))
    A.setdiag(aa[1:], k=-1)  # TODO: Check for the diagonal-1 matrix starts from - A[1,0]
    A.setdiag(bb, k=0)
    A.setdiag(cc[:-1], k=1)
    A = A.tocsr()
    B = np.asarray(dd)

    dc = spsolve(A, B).tolist()

    # Isotope fluxes
    qiso_in = sli.qprec(dt) * sli.cprec(dt) + qevapin * cevapin
    qiso_out = qevapout * (cevapout + sig * dc[0] * dcevapoutdciso)
    qiso_evap = qevapout * (cevapout + sig * dc[0] * dcevapoutdciso) - qevapin * cevapin
    qiso_trans = 0

    qiso_liq_diff, qiso_vap_diff, qiso_liq_adv, qiso_vap_adv = [], [], [], []

    for l in range(sli.n(dt)):

        if l != sli.n(dt) - 1:

            qiso_liq_diff.append(Dl_mean[l] * ((sli.ciso(dt)[l + 1] + sig * dc[l])
                                               - (sli.ciso(dt)[l + 2] + sig * dc[l + 1])) / sli.deltaz(dt)[l])
            qiso_vap_diff.append(Dv_mean[l] * ((sli.ciso(dt)[l + 1] + sig * dc[l]) * Beta[l]
                                               - (sli.ciso(dt)[l + 2]
                                                  + sig * dc[l + 1]) * Beta[l + 1]) / sli.deltaz(dt)[l])

            qiso_liq_adv.append(sli.qlsig(dt)[l + 1] * (cql[l] + sig * dc[l] * dcqldca[l]
                                                        + sig * dc[l + 1] * dcqldcb[l]))
            qiso_vap_adv.append(sli.qvsig(dt)[l + 1] * betaqv[l] * (cqv[l] + sig * dc[l] * dcqvdca[l]
                                                                    + sig * dc[l + 1] * dcqvdcb[l]))
        else:

            qiso_liq_adv.append(sli.qlsig(dt)[-1] * (cql[-1] + sig * dc[-1]))
            qiso_vap_adv.append(0)

    
    #####################################################################
    # Mass Balance
    LHS = []
    RHS = []
    for l in range(sli.n(dt)):

        # # Coefficient of tridiagonal matrix
        lhs = (sli.thetasat(dt)[l] - sli.thetar(dt)[l]) * sli.dx(dt)[l] / sli.dt(dt) * (dc[l] *Seff[l] + sli.ciso(dt)[l+1] * deltaSeff[l])
        LHS.append(lhs)
        
        if l == 0:  # first layer
            
            r = sli.qprec(dt) * sli.cprec(dt) - qevapout * (cevapout + dc[l] * dcevapoutdciso) + qevapin * cevapin \
                - sli.qlsig(dt)[l+1] * (cql[l] + dc[l] * dcqldca[l] + dc[l+1] * dcqldcb[l]) \
                - sli.qvsig(dt)[l+1] * betaqv[l] * (cqv[l] + dc[l] * dcqvdca[l] + dc[l+1] * dcqvdcb[l]) \
                - Dl_mean[l] *((sli.ciso(dt)[l+ 1] + dc[l]) - (sli.ciso(dt)[l+2] + dc[l+1])) / sli.deltaz(dt)[l] / sig \
                - Dv_mean[l] * ((sli.ciso(dt)[l+1] + dc[l]) * Beta[l] - (sli.ciso(dt)[l+2] + dc[l+1]) * Beta[l+1]) / sli.deltaz(dt)[l] / sig \
                - (sli.qex(dt)[l] + sli.qrunoff(dt)) * (sli.ciso(dt)[l+1] + dc[l]) / sig     
                
            RHS.append(r)
            
        elif l < sli.n(dt) - 1:  # Intermediate layers

            r = sli.qlsig(dt)[l] * (cql[l - 1] + dc[l - 1] * dcqldca[l -1] + dc[l] * dcqldcb[l-1]) / sig \
                + sli.qvsig(dt)[l] * betaqv[l - 1] * (cqv[l - 1] + dc[l - 1] * dcqvdca[l -1] + dc[l] * dcqvdcb[l-1])/ sig \
                - sli.qlsig(dt)[l + 1] * (cql[l] + dc[l] * dcqldca[l] + dc[l+1] * dcqldcb[l]) / sig \
                - sli.qvsig(dt)[l + 1] * betaqv[l] * (cqv[l] + dc[l] * dcqvdca[l] + dc[l+1] * dcqvdcb[l]) / sig \
                + Dl_mean[l - 1] * ((sli.ciso(dt)[l] + dc[l-1]) - (sli.ciso(dt)[l + 1] + dc[l])) / sli.deltaz(dt)[l - 1] / sig  \
                + Dv_mean[l - 1] * ((sli.ciso(dt)[l] + dc[l-1]) * Beta[l - 1] 
                                    - (sli.ciso(dt)[l + 1] + dc[l]) * Beta[l]) / sli.deltaz(dt)[l - 1] / sig  \
                - Dl_mean[l] * ((sli.ciso(dt)[l + 1] + dc[l]) - (sli.ciso(dt)[l + 2] + dc[l+1] )) / sli.deltaz(dt)[l] / sig  \
                - Dv_mean[l] * ((sli.ciso(dt)[l + 1] + dc[l]) * Beta[l] 
                                - (sli.ciso(dt)[l + 2] + dc[l+1]) * Beta[l + 1])/ sli.deltaz(dt)[l] / sig \
                + sli.qex(dt)[l] * (sli.ciso(dt)[l + 1] + dc[l]) / sig

            RHS.append(r)

            
        else:  # last layer            

            r =  sli.qlsig(dt)[l] * (cql[l - 1] + dc[l-1] * dcqldca[l-1] + dc[l] * dcqldcb[l-1]) / sig \
                + sli.qvsig(dt)[l] * betaqv[l - 1] * (cqv[l - 1] + dc[l - 1] * dcqvdca[l-1] + dc[l] * dcqvdcb[l-1] ) / sig \
                - sli.qlsig(dt)[l + 1] * (cql[l] +dc[l])/ sig \
                + Dl_mean[l - 1] * ((sli.ciso(dt)[l] + dc[l-1]) - (sli.ciso(dt)[l + 1] + dc[l])) / sli.deltaz(dt)[l - 1] / sig  \
                + Dv_mean[l - 1] * ((sli.ciso(dt)[l] + dc[l-1]) * Beta[l - 1] 
                                    - (sli.ciso(dt)[l + 1] + dc[l]) * Beta[l]) / sli.deltaz(dt)[l - 1] / sig \
                - sli.qex(dt)[l] * (sli.ciso(dt)[l + 1] + dc[l]) / sig         

            RHS.append(r)

    if sli.cali(dt) > 0 or sli.testcase(dt) == '7' or sli.testcase(dt) == '8':
        RHS[-1] = RHS[-1] - sli.qsig(dt)[-1] * (sli.cali(dt) - (cql[-1] + dc[-1])) 
        
    mass = max(abs(np.array(LHS) - np.array(RHS)))
    print(sli.testcase(dt), sli.dt(dt))
    
    Ciso = np.array(sli.ciso(dt)[1:]) + np.array(dc)
    c_delta = [concentration_to_delta(c_iso, solute) for c_iso in Ciso]

    return dc, Ciso, c_delta


pth = os.getcwd()
path = os.path.abspath(os.path.join(pth, "..", ".."))

ignorealphai, ignorealphaik, ignoredl, ignoredvi = False, False, False, False
ignore = [ignorealphai, ignorealphaik, ignoredl, ignoredvi]
sli = Sli.SlI(path + '\_sli_\sli_label3\iso_variables_7')  # imports all the variable files /variables folder: testcase-1, sig=1

d_c = []
ciso = []
cdelta = []
for dt in range(len(sli.get_in_soil())-1):
    
    print(dt)
   
    dc, c, d = solve_iso(sli, dt, ignore)
    #print(c.tolist())
    d_c.append(dc)
    ciso.append(c)
    cdelta.append(d)    
    #print(d, c.tolist())

print(ciso)

