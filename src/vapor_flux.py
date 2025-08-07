# -*- coding=utf-8 -*-
from math import exp, sqrt, log
import numpy as np

from .iso_fluxes import vapor_diffusion_base_class


class flux_node(vapor_diffusion_base_class):

    def __init__(self, top_node):

        super().__init__()

        self.top_layer = top_node
        self.hmin = -1e6   # minimum matric head h
        
        self.repeat = False
        self.dy = None
        self.nodes = None

    def dv(self, node):

        return self.dv_soil_air(T=node.T, theta=node.theta,
                                theta_sat=node.theta_sat,
                                tortuosity=node.tortuosity)

        # return self.dv_free_air(T=node.T, Pa=10**5)

    def kE(self, node):

        Tzero = 273.16000366210938

        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
        R = 8.3142995834350586  # universal gas constant (j / mol / k)

        psat_0 = self.p_sat_0(node)  # Saturation vapour pressure, ps, in pascal

        s_psat = psat_0 * 17.270000457763672 * 237.30000305175781 / ((node.T - Tzero) + 237.30000305175781) ** 2

        dv = self.dv_soil_air(T=node.T, theta=node.theta, theta_sat=node.theta_sat, tortuosity=node.tortuosity)
        rh = node.relative_humidity(node.psi, node.T)
        sl = s_psat * Mw / 1000 / R / node.T
        lambdav = self.lambdav(node)
        eta_th = 1

        ke = dv * rh * sl * 1000 * lambdav * eta_th

        return ke

    def lambdav(self, node):
            return 1.91846e6 * (node.T / (node.T - 33.91)) ** 2

    def phivS(self, node):

        dva = 2.1699999706470408E-005  # vapour diffusivity of water in air at 0 degC [m2/s]
        thre = node.theta_sat - node.theta_0

        dhds = self.dhds(node)
        c = self.c(node)
        h = self.h(node)

        ph = dhds * dva * node.tortuosity * \
             (thre * c * exp(c * h) - thre * c * (h / node.he) ** (-node.lamda) * exp(c * h))

        return ph

    def c(self, node):

        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
        R = 8.3142995834350586  # universal gas constant (j / mol / k)
        g = 9.8000001907348633  #

        return Mw * g / R / node.T

    def dhds(self, node):

        try:
            if node.S < 1:  # unsaturated

                return - node.he / node.lamda * node.S ** (-1 / node.lamda - 1)

            elif node.S >= 1:  # saturated

                return 0

            else:
                raise ValueError

        except ValueError:
            raise NotImplementedError

    def phiS(self, node):

        try:

            if node.S >= 1:
                return (node.eta - 1 / node.lamda) * node.phie
            elif node.S < 1:
                return (node.eta - 1 / node.lamda) * self.phi(node) / node.S
            else:
                raise ValueError

        except ValueError:
            raise NotImplementedError

    def KS(self, node):

        try:

            if node.S >= 1:
                return node.eta * node.ksat
            elif node.S < 1:
                v4 = exp(node.eta * log(node.S))
                K = node.ksat * v4
                return node.eta * K / node.S
            else:
                raise ValueError

        except ValueError:
            raise NotImplementedError

    def p_sat_0(self, node):

        Tzero = 273.16000366210938

        p_sat_0 = 610.59999465942383 * exp(
            17.270000457763672 * (node.T - Tzero) /
            ((node.T - Tzero) + 237.30000305175781))  # Saturation vapour pressure, ps, in pascal

        return p_sat_0

    def cvsat_0(self, node):

        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
        R = 8.3142995834350586

        cv = self.p_sat_0(node) * Mw / R / node.T / 1000  # m3/m3

        return cv

    def adjust_phi(self, node):

        idx = self.nodes.index(node)
        dy = self.dy[idx]

        phi = (node.phie * self.v3(node) * self.v4(node)) + dy

        if phi < node.phie:
            return node.phie
        else:
            return phi
        
    def phi(self, node):

        if not self.repeat:
            return node.phie * self.v3(node) * self.v4(node)
        else:
            if node.S >= 1:  # Saturated
                return self.adjust_phi(node)
            else:
                return node.phie * self.v3(node) * self.v4(node)
            
    def v4(self, node):
        return exp(node.eta * log(node.S))

    def v3(self, node):
        return exp(-log(node.S) / node.lamda)

    def K(self, node):

        try:
            if self.repeat and node == self.top_layer:
                return self.KZ(node, self.hmin)
            else:
                if node.S < 1:  # unsaturated
                    return node.ksat * self.v4(node)
                elif node.S >= 1:  # saturated
                    return node.ksat
                else:
                    raise ValueError

        except ValueError:
            raise NotImplementedError

    def KZ(self, node, h):

        return node.ksat * exp(self.a(node) * log(h / node.he))

    def phiz(self, node, h):

        return self.KZ(node, h) * h / (1 + self.a(node))

    def a(self, node):

        return - node.lamda * node.eta

    def h(self, node):

        try:

            if self.S0(node) >= 1:  # saturated

                return node.he

            elif self.S0(node) < 1:  # unsaturated

                return node.he * self.v3(node)
            else:
                raise ValueError

        except ValueError:
            raise NotImplementedError

    def S0(self, node):

        return (node.theta - node.theta_0) / (node.theta_sat - node.theta_0)


class vapor_flux(flux_node):

    def __init__(self, left_node, right_node, top_layer):
        """
        Class to calculate the vapor flux based on a given temprature gradient and hydraulic head (psi in m) gradient.
        Based on Haverd&Cuntz (2005) Eq. A.12
        """

        super().__init__(top_layer)

        self.left_node = left_node
        self.right_node = right_node
        self.gf = 1  # gravity factor i.e. flow downwards

    def rdz(self):

        return self.left_node.center.distance(self.right_node.center)

    def qvapor(self):

        dy_left = self.left_node.theta - self.left_node.theta_t0
        dy_right = self.right_node.theta - self.right_node.theta_t0

        # dqv = self.qvya() * dy_left + self.qvyb() * dy_right
        # dqt = self.qTa() * self.left_node.dT + self.qTb() * self.right_node.dT

        return self.q_v()  # + dqv + dqt

    def q(self):
        return self.q_v() + self.q_l()

    def q_l(self):

        gf = 1
        macropore_factor = 1
        hmin = -1.0e6

        w = self.weight(gf=gf)
        if self.left_node == self.top_layer:
            K1 = self.KZ(self.left_node, h=hmin)
        else:
            K1 = self.K(self.left_node)

        K2 = self.K(self.right_node)
        ph1 = self.phi(self.left_node)
        ph2 = self.phi(self.right_node)

        return (ph1 - ph2) / self.rdz() +\
            gf * (w * K1 * macropore_factor + (1 - w) * K2 * macropore_factor)

    def q_v(self):

        dv_left, dv_right = self.dv(self.left_node), self.dv(self.right_node)

        try:

            if dv_left == 0 or dv_right == 0:
                qv = 0.0
            else:
                qv = self.qvh() + self.qvT()

        except ValueError:
            raise NotImplementedError

        return qv

    def qvT(self):

        q = (self.left_node.T - self.right_node.T) * (self.kE(self.left_node) + self.kE(self.right_node)) \
            / 1000 / self.lambdav(self.top_layer) / 2 / self.rdz()

        return q

    def qvh(self):

        node_l, node_r = self.left_node, self.right_node
        cv_l, cv_r = node_l.cv_sat, node_r.cv_sat

        rh_l = node_l.relative_humidity(node_l.psi, node_l.T)
        rh_r = node_r.relative_humidity(node_r.psi, node_r.T)

        mean_dv = (self.dv(self.left_node) + self.dv(self.right_node)) / 2
        mean_cvsat = (cv_l + cv_r) / 2
        mean_rh = (rh_l + rh_r) / 2

        dz = self.rdz()

        return mean_dv * mean_cvsat * (rh_l - rh_r) / dz

    def qya(self):
        return self.qlya() + self.qvya()

    def qyb(self):
        return self.qlyb() + self.qvyb()

    def qTa(self):

        qta = self.qta()

        q = qta + (self.kE(self.left_node) + self.kE(self.right_node)) / 1000 / \
            self.lambdav(self.top_layer) / 2 / self.rdz()

        return q

    def qTb(self):

        qtb = self.qtb()

        q = qtb - (self.kE(self.left_node) + self.kE(self.right_node)) / 1000 / \
            self.lambdav(self.top_layer) / 2 / self.rdz()

        return q

    def qvya(self):

        try:

            dv_left, dv_right = self.dv(self.left_node), self.dv(self.right_node)

            if dv_left == 0 or dv_right == 0:
                q = 0.0
            else:

                S = self.left_node.S
                if S < 1:  # unsaturated

                    Tzero = 273.16000366210938

                    phivs = self.phivS(self.left_node)
                    dz = self.rdz()

                    q = phivs / dz * (((self.left_node.T / Tzero) ** 1.88 + (self.right_node.T / Tzero) ** 1.88) / 2) * \
                        ((self.cvsat_0(self.left_node) + self.cvsat_0(self.right_node)) / 2)

                elif S >= 1:  # saturated

                    q = 0.0
                else:
                    raise ValueError

            return q

        except ValueError:
            raise NotImplementedError

    def qvyb(self):

        try:

            dv_left, dv_right = self.dv(self.left_node), self.dv(self.right_node)

            if dv_left == 0 or dv_right == 0:
                q = 0.0
            else:

                S = self.left_node.S

                if S < 1:  # unsaturated

                    Tzero = 273.16000366210938

                    phivs = self.phivS(self.right_node)
                    dz = self.rdz()

                    q = - phivs / dz * (((self.left_node.T / Tzero) ** 1.88 + (self.right_node.T / Tzero) ** 1.88) / 2) * \
                        ((self.cvsat_0(self.left_node) + self.cvsat_0(self.right_node)) / 2)

                elif S >= 1:  # saturated
                    q = 0.0
                else:
                    raise ValueError

            return q

        except ValueError:
            raise NotImplementedError

    def qlya(self):

        try:
            node = self.left_node
            S = node.S

            if S < 1:  # unsaturated

                gf = 1
                macropore_factor = 1

                w = self.weight(gf=gf)

                qya = self.phiS(node) / self.rdz() + gf * w * self.KS(node) * macropore_factor

            elif S >= 1:  # saturated

                qya = 1 / self.rdz()

            else:
                raise ValueError

            return qya

        except ValueError:
            raise NotImplementedError

    def qta(self):

        try:
            S = self.left_node.S

            if S < 1:  # unsaturated
                gf = 1
                macropore_factor = 1
                phiT = 0
                kT = 0
                w = self.weight(gf=gf)

                qTa = phiT / self.rdz() + gf * w * kT * macropore_factor

            elif S >= 1:  # saturated

                qTa = 0

            else:
                raise ValueError

            return qTa

        except ValueError:
            raise NotImplementedError

    def qlyb(self):

        try:
            node = self.right_node
            S = node.S
            if S < 1:  # unsaturated

                gf = 1
                macropore_factor = 1
                w = self.weight(gf=gf)

                qyb = -self.phiS(node) / self.rdz() + gf * (1 - w) * self.KS(node) * macropore_factor

            elif S >= 1:  # saturated

                qyb = - 1 / self.rdz()

            else:
                raise ValueError

            return qyb

        except ValueError:
            raise NotImplementedError

    def qtb(self):

        try:

            S = self.left_node.S

            if S < 1:  # unsaturated
                gf = 1
                macropore_factor = 1
                phiT = 0
                kT = 0
                w = self.weight(gf=gf)

                qTb = -phiT / self.rdz() + gf * (1 - w) * kT * macropore_factor

            elif S >= 1:  # saturated

                qTb = 0

            else:
                raise ValueError

            return qTb

        except ValueError:
            raise NotImplementedError

    def weight(self, gf=1):

        sl, sr = self.left_node.S, self.right_node.S
        h = self.h(self.right_node)
        dz = self.rdz()
        f = h - gf * dz

        try:
            if (sl == 1 and sr == 1) or f >= self.left_node.he:
                return 0.0
            else:
                return self.wt(node=self.right_node, gf=gf)

        except ValueError:
            raise NotImplementedError

    def wt(self, node, gf=1):

        try:
            h = self.h(node)
            K = self.K(node)
            dz = self.rdz()
            hz = h - gf * dz

            if h < node.he:

                a = node.lamda * node.eta
                x = - gf * dz / h

                if a <= 3.0 or x * (a - 3.0) < 4.0:
                    # predetermined approximation
                    w = (60 + x * (70 + 10 * a + x * (16 + a * (5 + a)))) /\
                        (120 + x * (120 + x * (22 + 2 * a ** 2)))

                else:
                    # more accurate slower method

                    w = -((self.phiz(node, hz) - self.phi(node)) / (gf * dz) + K) /\
                        (self.KZ(node, hz) - K)

            else:

                w = -((self.phiz(node, hz) - self.phi(node)) / (gf * dz) + K) / \
                    (self.KZ(node, hz) - K)

        except ValueError:
            raise NotImplementedError

        return min(max(w, 0), 1)


class surface_fluxes(vapor_flux):

    def __init__(self, pond, top_layer):

        super().__init__(left_node=pond, right_node=top_layer, top_node=top_layer)

        self.pond = pond
        self.top_layer = top_layer

    def set_ponding(self, n_nodes, nsat):

        node = self.top_layer

        if node.pond is None:
            h0 = 0
        elif node.pond is not None:
            h0 = node.pond.pond_height
        else:
            raise ValueError

        if self.phi(node) <= self.phip(node) and h0 <= 0 and nsat < n_nodes or node.S0 < 1:
            self.ponding = False
        else:
            self.ponding = True

    def phip(self, node):

        # mfp pond

        x1 = node.phie - node.he * node.ksat
        x2 = (1 + 1e-5) * node.phie

        return max(x1, x2)

    def relative_humidity(self):
        return 1

    def set_theta(self):

        self.pond.theta = self.top_layer.theta
        self.pond.theta_sat = self.top_layer.theta_sat


class vapor_solve(object):

    def delta_cc(self, node, dt):

        if node.S < 1:  # unsaturated

            return (node.theta_sat - node.theta_0) * node.thickness * 1 / dt

        elif node.S >= 1:  # saturated
            return 0
        else:
            raise ValueError

    def thomas_scalar(self, a, b, c, d):
        """
        Solves tridiagonal system using Thomas algorithm.
        a: sub-diagonal (length n-1)
        b: diagonal (length n)
        c: super-diagonal (length n-1)
        d: right-hand side (length n)
        Returns: x (solution vector, length n)
        """
        n = len(b)
        cp = np.zeros(n - 1)
        dp = np.zeros(n)

        # Forward sweep
        cp[0] = c[0] / b[0]
        dp[0] = d[0] / b[0]
        for i in range(1, n - 1):
            denom = b[i] - a[i - 1] * cp[i - 1]
            cp[i] = c[i] / denom
            dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / denom
        dp[n - 1] = (d[n - 1] - a[n - 2] * dp[n - 2]) / (b[n - 1] - a[n - 2] * cp[n - 2])

        # Back substitution
        x = np.zeros(n)
        x[-1] = dp[-1]
        for i in range(n - 2, -1, -1):
            x[i] = dp[i] - cp[i] * x[i + 1]

        return x

    def solve_sparse(self, X):
        """
        Solves a scalar tridiagonal system for dy.
        Inputs:
            aa, bb: sub-diagonal terms (length n-1)
            cc: diagonal terms (length n)
            ee, ff: super-diagonal terms (length n-1)
            gg: right-hand side (length n)
        Returns:
            dy: solution vector (length n)
        """
        aa, bb, cc, dd, ee, ff, gg = X
        n = len(cc)

        # Input checks (optional)
        assert len(aa) == n - 1
        assert len(bb) == n - 1
        assert len(ee) == n - 1
        assert len(ff) == n - 1
        assert len(gg) == n

        # Tridiagonal matrix coefficients (combine multiple inputs if needed)
        a = aa + bb  # sub-diagonal
        b = cc  # diagonal
        c = ee + ff  # super-diagonal
        d = gg  # RHS

        # Solve
        dy = self.thomas_scalar(a, b, c, d)
        return dy

