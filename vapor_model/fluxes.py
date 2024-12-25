
from math import sqrt


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

    def q_liquid(self):
        pass

    def hy_conduct_potential(self):

        return abs(sqrt(self.left_node.kml * self.right_node.kml))

    def hy_conduct_temp(self):

        return abs(sqrt(self.left_node.dtl * self.right_node.dtl)) * (self.left_node.T - self.right_node.T)


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
