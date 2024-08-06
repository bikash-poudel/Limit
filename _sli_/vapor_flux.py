# -*- coding=utf-8 -*-
from math import exp, sqrt

from iso_fluxes import vapor_diffusion_base_class


class vapor_flux(object):
    def __init__(self):
        """
        Class to calculate the vapor flux based on a given temprature gradient and hydraulic head (psi in m) gradient.
        Based on Haverd&Cuntz (2005) Eq. A.12
        """
        pass

    def qv_vapor(self, left_node, right_node):

        vd = vapor_diffusion_base_class()

        distance = left_node.center.distance(to_Point=right_node.center)

        left_node_dv_soil = vd.dv_soil_air(T=left_node.T,
                                           theta=left_node.theta,
                                           theta_sat=left_node.theta_sat,
                                           tortuosity=left_node.tortuosity)

        right_node_dv_soil = vd.dv_soil_air(T=right_node.T,
                                            theta=right_node.theta,
                                            theta_sat=right_node.theta_sat,
                                            tortuosity=right_node.tortuosity)

        mean_dv_soil = (left_node_dv_soil * left_node.thickness + right_node_dv_soil * right_node.thickness) / distance

        # cv_sat
        mean_cv_sat = (left_node.cv_sat * left_node.thickness + right_node.cv_sat * right_node.thickness) / distance

        # hr
        left_node_hr = left_node.relative_humidity(psi=left_node.psi, T=left_node.T)
        right_node_hr = right_node.relative_humidity(psi=right_node.psi, T=right_node.T)

        mean_hr = (left_node_hr * left_node.thickness + right_node_hr * right_node.thickness) / distance

        q_vh = mean_dv_soil * mean_cv_sat * (left_node_hr - right_node_hr) / distance

        q_vT = mean_dv_soil * mean_hr * (left_node.cv_sat - right_node.cv_sat) / distance

        q_v = q_vh + q_vT  # whole vapour flux has one part from humidity (q_vh) and one part from temp diff (q_vT) ; concersion from sec to days

        return q_v

    def q_vapor(self, left_node, right_node):
        """
        Function to calculate the vapor flux (m/day) based on a given temprature gradient and hydraulic head (psi in m) gradient.
        Based on Haverd&Cuntz (2005) Eq. A.12. as implemented in SLI::cable_sli_utils.f90 Lines 626-636


        qvh(i) = ((var(i)%Dv+var(i+1)%Dv)/two)* ((var(i)%cvsat+var(i+1)%cvsat)/two)*(var(i)%rh-var(i+1)%rh)/dz(i)
        qvT(i) = (Tsoil(i)-Tsoil(i+1))*(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
        qv(i)  = qvh(i) + qvT(i) ! whole vapour flux has one part from humidity (qvh) and one part from temp diff (qvT)

        var%lambdav   = 1.91846e6_r_2*((Tsoil+Tzero)/((Tsoil+Tzero)-33.91_r_2))**2

        var%kE     = var%Dv*var%rh*var%sl*thousand*var%lambdav*var%eta_th
        var%sl = slope_esat(Tsoil) * Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1
        var%eta_th = one
        """

        Tzero = 273.16000366210938

        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
        R = 8.3142995834350586  # universal gas constant (j / mol / k)

        eta_th = 1

        dz = left_node.center.distance(right_node.center)

        # slope of psat
        s_psat_left = left_node.p_sat * 17.270000457763672 * 237.30000305175781 / \
                  ((left_node.T - Tzero) + 237.30000305175781) ** 2

        s_psat_right = right_node.p_sat * 17.270000457763672 * 237.30000305175781 / \
                   ((left_node.T - Tzero) + 237.30000305175781) ** 2

        # var sl sli_utils, l:1405
        sl_left = s_psat_left * Mw / 1000 / R / left_node.T
        sl_right = s_psat_right * Mw / 1000 / R / right_node.T

        vd = vapor_diffusion_base_class()

        dv_left = vd.dv_soil_air(left_node.T, left_node.theta, left_node.theta_sat, left_node.tortuosity)
        dv_right = vd.dv_soil_air(right_node.T, right_node.theta, right_node.theta_sat, right_node.tortuosity)
        mean_dv = (dv_left + dv_right) / 2

        rh_left = left_node.relative_humidity(left_node.psi, left_node.T)
        rh_right = right_node.relative_humidity(right_node.psi, right_node.T)
        mean_rh = (rh_left - rh_right) / dz

        mean_cv_sat = (left_node.cv_sat + right_node.cv_sat) / 2

        # latent heat of sublimation (vaporization)
        # sli_utils, l:1306
        lambdav_left = 1.91846e6 * (left_node.T / (left_node.T - 33.91)) ** 2
        lambdav_right = 1.91846e6 * (right_node.T / (right_node.T - 33.91)) ** 2

        kE_left = dv_left * rh_left * sl_left * lambdav_left * 1000 * eta_th
        kE_right = dv_right * rh_right * sl_right * lambdav_right * 1000 * eta_th

        qv_h = mean_dv * mean_cv_sat * mean_rh
        qv_T = (left_node.T - right_node.T) * (kE_left + kE_right) / 1000 / lambdav_left / 2 / dz

        """
        dv_left = self.dv_soil_air(left_node.T, left_node.theta, left_node.theta_sat, left_node.tortuosity)
        dv_right = self.dv_soil_air(right_node.T, right_node.theta, right_node.theta_sat, right_node.tortuosity)

        hr_left = left_node.relative_humidity(left_node.psi, left_node.T)
        hr_right = right_node.relative_humidity(right_node.psi, right_node.T)

        mean_dv = (dv_left * left_node.thickness + dv_right * right_node.thickness) / dz
        mean_hr = (hr_left * left_node.thickness + hr_right * right_node.thickness) / dz
        mean_cvsat = (left_node.cv_sat * left_node.thickness + right_node.cv_sat * right_node.thickness) / dz

        q_vh = mean_dv * mean_cvsat * (hr_left - hr_right) / dz
        q_vT = mean_dv * mean_hr * (left_node.cv_sat - right_node.cv_sat) / dz
        """

        q_v = qv_h + qv_T

        return q_v

    def qvl_vapor(self, left_node, right_node):

        """"Returns the vapor flux between two nodes
        # Fayer: 2000, Saito et al.:2006,

        """
        dz = left_node.center.distance(right_node.center)

        #  Isothermal vapor conductivity
        kvh_mean = (self.K_vh(left_node) * left_node.thickness + self.K_vh(right_node) * right_node.thickness) / dz

        # Thermal vapor conductivity
        kvT_mean = (self.K_vT(left_node) * left_node.thickness + self.K_vT(right_node) * right_node.thickness) / dz

        flux = kvh_mean + kvT_mean

        return flux

    def K_vh(self, node):
        """"Returns isothermal vapor conductivity"""
        # Fayer: 2000, Saito et al.:2006,

        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
        R = 8.3142995834350586  # universal gas constant (j / mol / k)
        g = 9.81  # m s-1
        rhow = 1000  # density of liquid water (kg m-3)

        vd = vapor_diffusion_base_class()
        dv = vd.dv_soil_air(node.T, node.theta, node.theta_sat, node.tortuosity)
        rh = node.relative_humidity(node.psi, node.T0)

        rho_vsat = self.rhow_v_sat(node.T0)
        kvh = - dv * Mw * g * rho_vsat / rhow / R / node.T0 * rh

        return kvh

    def K_vT(self, node):
        """"Returns thermal vapor conductivity"""
        # Fayer: 2000, Saito et al.:2006

        rhow = 1000  # density of liquid water (kg m-3)

        vd = vapor_diffusion_base_class()
        dv = vd.dv_soil_air(node.T, node.theta, node.theta_sat, node.tortuosity)
        rh = node.relative_humidity(node.psi, node.T0)
        n = self.n(node)  # enhancement factor

        kvT = - dv / rhow * n * rh * self.d_dT_rhow_vsat(node.T0)

        return kvT

    def rhow_v_sat(self, T):
        """"Returns saturated vapour density [kg m-3]"""
        # Fayer: 2000, Saito et al.:2006

        rhow_sv = exp(31.3716 - (6014.79 / T) - 7.92495 * 10 ** -3 * T) / T / 1000

        return rhow_sv

    def d_dT_rhow_vsat(self, T):
        """"Returns d/dT(rhow_sat_v) saturated vapour density """
        # Fayer: 2000, Saito et al.:2006
        d_dT = self.rhow_v_sat(T) * (6014 / T - 7.92495 * 10 ** -3 * T - 1)

        return d_dT

    def n(self, node, fc=0.3):
        """"Enhancement factor

        equation for the enhancement factor was derived by Cass et al. (1984)
        Saito et al.(2006)

        @fc: caly fraction [-]

        """
        n = 9.5 + 3 * (node.theta_sat - node.theta / node.theta_sat) \
            - 8.5 * exp(-((1 + 2.6 / sqrt(fc)) * (node.theta_sat - node.theta / node.theta_sat)) ** 4)

        return n

