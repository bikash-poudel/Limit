# -*- coding=utf-8 -*-
from math import exp, sqrt


class heat_storage(object):

    def __init__(self, node):

        # copy all attributes from node
        for key, value in vars(node).items():
            setattr(self, key, value)

        self.node = node
        self.head = node.psi * 100  # into cm
        self.rho_b = node.rho * 0.001  # soil bulk density  into g cm-3
        self.rho_l = 1.0  # density of water g cm-3

    @property
    def T1(self):

        """ Thermal_1
        - we defined here as T1
        - unit J cm-3 K-1
        - T1 is not as same as the coefficient of the pTpt.
        """

        # preliminary parameters

        L0 = 2500 # 2270  # J g-1
        cv = 1.864  # J g-1 K-1

        theta_a = self.theta_sat - self.theta
        Tt = self.T

        rho_s = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # g cm-3
        HR = exp(2.1238e-4 * self.head / Tt)
        # rho_v = rho_s * HR
        prho_SpT = rho_s * 4975.9 / (Tt ** 2)

        pHRpT = -HR * 2.1238e-4 * self.head / (Tt ** 2)
        prho_vpT = rho_s * pHRpT + prho_SpT * HR

        t1 = self.cs + (L0 + cv * (self.T - self.T0)) * theta_a * prho_vpT

        return t1

    @property
    def cs(self):

        """soil volumetric head capacity
        - unit J cm-3 K-1
        """

        # preliminary parameters for each partitions
        c_o = 1.900  # "J g^-1 K-1"
        c_m = 0.730  # "J g^-1 K-1"
        c_v = 1.864  # "J g^-1 K-1"
        c_l = 4.187  # "J g^-1 K-1"

        # weight of each partition
        phi_o = self.som  # soil organic matter
        phi_m = self.sand + self.silt + self.clay

        # potential - temperature correction
        ht = self.head
        Tt = self.T
        theta_a = self.theta_sat - self.theta
        rho_vs = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # gm cm-3
        HR = exp(2.1238e-4 * ht / Tt)
        rho_v = rho_vs * HR  # gm cm-3

        c = self.rho_b * (c_o * phi_o + c_m * phi_m) + c_v * rho_v * theta_a + c_l * self.rho_l * self.theta

        return c

    @property
    def lambda_bulk(self):

        """Thermal conductivity (lambda)
        - unit "w cm-1 K-1"
        """

        lambda_dry = - 0.56 * self.theta_sat + 0.51
        alpha = 0.67 * self.clay + 0.24
        beta = 1.97 * self.sand + 1.87 * self.rho_b - 1.36 * self.sand * self.rho_b - 0.95

        lam = 0.01 * (lambda_dry + exp(beta - self.theta ** (-alpha)))  # convert to

        return lam

    @property
    def latent_liquid(self):

        """ J cm-3 """

        rho_l = 1  # g cm-3 density of water
        L_0 = 2270  # [J g-1]
        c_l = 4.187  # Specific Heat at Constant Temperature [J g-1 K-1]

        return rho_l * c_l * (self.T - self.T0)

    @property
    def latent(self):

        """ J cm-3 """

        rho_l = 1  # g cm-3 density of water
        L_0 = 2270  # [J g-1]
        c_v = 1.864  # Specific Heat at Constant Temperature [J g-1 K-1]

        return rho_l * (L_0 + c_v * (self.T - self.T0))


class heat_fluxes(object):

    def __init__(self, left_node, right_node):
        self.left_node = left_node
        self.right_node = right_node

        # register connection
        left_node.RegisterConnection(newConnection=self)
        right_node.RegisterConnection(newConnection=self)


class heat_flux_connection(object):

    def __init__(self, left_node, right_node):
        self.left_node = left_node
        self.right_node = right_node

    def th_conductivity(self):
        return abs(sqrt(self.left_node.lambda_bulk * self.right_node.lambda_bulk))

    def q_cond(self):

        """ W cm-2 [J cm-2 s-1] """

        # dz = self.left_node.node.center.distance(self.right_node.node.center) * 100  # into cm

        return self.th_conductivity() * (self.left_node.T - self.right_node.T)  # / dz


class heat_solve(object):

    def __init__(self, storages):

        self.storages = storages

    def solve(self, ql, qv, qh, dt):

        ql_latent = [s.latent_liquid for s in self.storages]
        qv_latent = [s.latent for s in self.storages]

        # first layer / top layer
        T = [self.storages[0].T]
        for i, s in enumerate(self.storages[1:-1], start=1):

            ql_left, ql_right = ql[i - 1] * dt, ql[i] * dt
            qv_left, qv_right = qv[i - 1] * dt, qv[i] * dt
            q_cond_left, q_cond_right = qh[i - 1] * dt, qh[i] * dt

            # latent liquid
            if ql_right > 0:
                ql_latent_right = ql_latent[i + 1] * ql_right
            elif ql_right < 0:
                ql_latent_right = ql_latent[i] * ql_right
            else:
                ql_latent_right = 0

            if ql_left > 0:
                ql_latent_left = ql_latent[i] * ql_left
            elif ql_left < 0:
                ql_latent_left = ql_latent[i - 1] * ql_left
            else:
                ql_latent_left = 0

            # latent vapor
            if qv_right > 0:
                qv_latent_right = qv_latent[i + 1] * qv_right
            elif qv_right < 0:
                qv_latent_right = qv_latent[i] * qv_right
            else:
                qv_latent_right = 0

            if qv_left > 0:
                qv_latent_left = qv_latent[i] * qv_left
            elif qv_left < 0:
                qv_latent_left = qv_latent[i - 1] * qv_left
            else:
                qv_latent_left = 0

            tmp = s.T + (q_cond_right - q_cond_left
                         + ql_latent_right - ql_latent_left
                         + qv_latent_right - qv_latent_left) / s.T1

            T.append(tmp)

        T.append(self.storages[-1].T)  # last layer / bottom layer

        return T

