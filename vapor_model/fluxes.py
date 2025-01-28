
from math import sqrt, exp


class flux_connection(object):

    def __init__(self, left_node, right_node):
        self.left_node = left_node
        self.right_node = right_node

        # register connection
        left_node.RegisterConnection(newConnection=self)
        right_node.RegisterConnection(newConnection=self)


class water_flux(flux_connection):

    def __init__(self, left_node, right_node):
        flux_connection.__init__(self, left_node, right_node)

    def hy_conduct_potential(self):

        return abs(sqrt(self.left_node.kml * self.right_node.kml))

    def q_pot_water(self):
        return self.hy_conduct_potential() * (self.right_node.head - self.left_node.head)

    def hy_conduct_temp(self):

        return abs(sqrt(self.left_node.dtl * self.right_node.dtl))

    def q_tmp_water(self):
        return self.hy_conduct_temp() * (self.right_node.T - self.left_node.T)

    def q_water(self):
        return self.q_pot_water() + self.q_tmp_water()

    def ql_sli(self):

        w = 0.5  # weight
        gf = 1  # gravity factor / direction of gravity
        mpf = 1  # macropore factor

        q = (self.right_node.ph_i - self.left_node.ph_i) * self.rdz() + \
            gf * (w * self.right_node.K * mpf + w * self.left_node.K * mpf)

        return q

    def qya(self):

        gf = 1
        w = 0.5
        mpf = 1
        if self.left_node.S < 1:
            return self.left_node.phi_S * self.rdz() + gf * w * self.left_node.KS * mpf
        elif self.left_node.S >= 1:
            return self.rdz()
        else:
            raise NotImplementedError

    def qyb(self):

        gf = 1
        w = 0.5
        mpf = 1
        if self.right_node.S < 1:
            return -self.right_node.phi_S * self.rdz() + gf * w * self.right_node.KS * mpf
        elif self.right_node.S >= 1:
            return -self.rdz()
        else:
            raise NotImplementedError

    def rdz(self):
            return 1 / self.left_node.distance(self.right_node)


class heat_flux(flux_connection):

    def __init__(self, left_node, right_node):
        flux_connection.__init__(self, left_node, right_node)

    def th_conductivity(self):

        return abs(sqrt(self.left_node.lambda_bulk * self.right_node.lambda_bulk))

    def hy_conduct_potential(self):

        return abs(sqrt(self.left_node.kml * self.right_node.kml))

    def th_conductivity_pot(self):

        return self.hy_conduct_potential() * (self.right_node.head - self.left_node.head)

    def th_hy_conductivity_correction(self):

        rho_w = 1  # g cm-3 density of water
        c_l = 4.187  # Specific Heat at Constant Pressure [J g-1 K-1]

        return self.th_conductivity_pot() * rho_w * c_l

    def th_hy_conductivity_T0(self):

        T0 = self.left_node.T0

        return self.th_hy_conductivity_correction() * T0

    def q_cond(self):

        return self.th_conductivity() * (self.right_node.T - self.left_node.T)


class vapor_flux(flux_connection):

    def __init__(self, left_node, right_node):
        flux_connection.__init__(self, left_node, right_node)

    def v_cond_pot(self):

        return abs(sqrt(self.left_node.Dmv * self.right_node.Dmv))

    def v_cond_temp(self):

        return abs(sqrt(self.left_node.Dtv * self.right_node.Dtv))

    def q_vapor_pot(self):

        return self.v_cond_pot() * (self.right_node.head - self.left_node.head)

    def q_vapor_tmp(self):

        return self.v_cond_temp() * (self.right_node.T - self.left_node.T)

    def q_vapor(self):

        return self.q_vapor_pot() + self.q_vapor_tmp()

    def q_vapor_latent_right(self):

        rho_l = 1  # g cm-3 density of water
        L_0 = 2270  # [J g-1]
        c_v = 1.864  # Specific Heat at Constant Temperature [J g-1 K-1]

        return rho_l * (L_0 + c_v * (self.right_node.T - self.left_node.T0))


class evaporation(flux_connection):

    def __init__(self, atmosphere, soil_layer, q=0):

        flux_connection.__init__(self, atmosphere, soil_layer)
        self.__q_evap = q

    @property
    def q_evap(self):

        return self.__q_evap

    @q_evap.setter
    def q_evap(self, qevap):
        self.__q_evap = qevap

    def ev_potential(self):

        Tzero = 273.15

        # sli Utils
        rmair = 0.02897  # kg / mol
        rgas = 8.3143  # J/mol/K
        capp = 1004.64  # air spec.heat capacity(J / kg / K)
        cpa = capp  # specific heat capacity of dry air at 0-40 degC [J/kgK]
        lambdav = 1.91846e6*((self.right_node.T +Tzero)/((self.right_node.T +Tzero)-33.91)) ** 2  ##  Henderson-Sellers, QJRMS, 1984
        rmh2o = 0.018016   # molecular wt: water(kg / mol)
        tetena = 6.106
        esata = tetena * 100.0  # constants for saturated vapour pressure calculation
        esatb= 17.27
        esatc = 237.3

        rhocp = rmair * 101325 / rgas / (self.left_node.T + Tzero) * cpa
        gamma = 101325 * cpa / lambdav / (rmh2o / rmair)
        esat = esata * exp(esatb*self.right_node.T * (self.right_node.T + esatc))
        slope_esat = esat * esatb * esatc/(self.right_node.T + esatc) ** 2
        ea = esata * exp(esatb * self.left_node.T * (self.left_node.T + esatc))
        Da = ea / self.left_node.Rh - ea

