# -*- coding=utf-8 -*-

from iso_fluxes import vapor_diffusion


class vapor_flux(object):
    def __init__(self):
        """
        Class to calculate the vapor flux based on a given temprature gradient and hydraulic head (psi in m) gradient.
        Based on Haverd&Cuntz (2005) Eq. A.12
        """
        pass

    @classmethod
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

        vd = vapor_diffusion(left_node, right_node)

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

