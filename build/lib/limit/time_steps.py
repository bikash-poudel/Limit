
import numpy as np

from . heat_fluxes import heat_solver


class timeSteps(heat_solver):

    def __init__(self, layers, ts=0, tfin=1800, dt_max=3600):

        super().__init__()

        self.layers = layers
        self.dtmax = dt_max
        self.tfin = tfin
        self.ts = ts
        self.dsmax = 0.1
        self.dtlmax = 0.1
        self.T_step0 = 0
        self.T_step = 0
        self.nsat_0 = 0
        self.nsat = 0

    def set_nsat(self):
        nsat = sum([l.S >= 1 for l in self.layers])
        self.nsat = nsat

    def repeat(self):
        if self.nsat > 0:
            return True
        else:
            return False

    def dt(self, dS_max, dT_max):

        try:
            if dS_max > 0:

                dt = min(self.dsmax / dS_max, self.dtlmax / dT_max,  self.dtmax)

            else:  # steady state flow

                dt = min(self.tfin - self.ts, self.dtlmax / dT_max, self.dtmax)

            if dt > self.dtmax:
                dt = self.dtmax

            return dt

        except ValueError:
            raise NotImplementedError

    def dt_moisture(self, q):

        """ [m / s] """

        if not len(q) + 1 == len(self.layers):
            raise 'check length of q'

        dm = np.zeros(len(self.layers))
        S = [l.S < 1 for l in self.layers]
        thre_dx = [(l.theta_sat - l.theta_0) * l.thickness for l in self.layers]

        qs = np.abs(np.array(q[1:]) - np.array(q[:-1])) / np.array(thre_dx)
        dm[S] = qs[S]

        return np.concatenate(([0], dm))

    def dt_temperature(self, qh):
        """ [K / s] """

        if not len(qh) + 1 == len(self.layers):
            raise 'check length of qh'

        csoil = [self.csoileff(l) * l.thickness for l in self.layers]
        qhs = np.abs(np.array(qh[1:]) - np.array(qh[:-1])) / np.array(csoil)

        return np.concatenate(([0], qhs))

    def dS_max(self, dt_moisture):

        return max(dt_moisture)

    def dT_max(self, dt_temperature):

        return max(dt_temperature)


