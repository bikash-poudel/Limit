from math import exp, log

import numpy as np

# from .iso_fluxes import vapor_diffusion_base_class
from .vapor_flux import vapor_flux, flux_node


class heat_node(flux_node):

    def __init__(self, top_layer):

        flux_node.__init__(self, top_layer)

    def kth(self, node):

        return self.kH(node) + self.kE(node)

    def kH(self, node):

        A = 0.65 - 0.78 * node.rho / 1000 + 0.6 * (node.rho / 1000) ** 2
        C1 = 1 + 2.6 * node.clay ** (-0.5)
        B = 2.8 * (1 - (node.theta_sat - node.theta_0)) * node.theta
        D = 0.03 + 0.7 * (1 - (node.theta_sat - node.theta_0)) ** 2
        E = 4

        kh = (A + B * node.theta - (A - D) * exp(- (C1 * node.theta) ** E))

        return kh

    def phiv(self, node):

        dva = 2.1699999706470408E-005  # vapour diffusivity of water in air at 0 degC [m2/s]
        thre = node.theta_sat - node.theta_0
        c = self.c(node)
        h = self.h(node)
        int = self.intt(node)

        return dva * node.tortuosity * (thre * exp(c * h) - thre * int)

    def intt(self, node):

        c = self.c(node)
        h = self.h(node)
        a, x = 1 - node.lamda, -c * h
        igamma = self.digamma(a, x)

        return (-c * node.he) ** node.lamda * igamma

    def digamma(self, a, x):

        gln = self.dgammln(a)

        if x < a + 1.0:

            igamma = 1.0 - self.dgser(a, x)
        else:
            igamma = self.dgcf(a, x)

        return igamma * exp(gln)

    def dgammln(self, z):

        lnsqrt2pi = 0.9189385332046727
        a = [0.9999999999995183, 676.5203681218835, -1259.139216722289, 771.3234287757674, -176.6150291498386,
             12.50734324009056, -0.1385710331296526, 0.9934937113930748E-05, 0.1659470187408462E-06]

        tmpgammln = 0
        tmp = z + 7.0

        for j in reversed(a[1:]):
            tmpgammln += j / tmp
            tmp -= 1

        tmpgammln = tmpgammln + a[0]
        gammln = log(tmpgammln) + lnsqrt2pi - (z + 6.5) + (z - 0.5) * log(z + 6.5)

        return gammln

    def dgser(self, a, x):

        eps = np.finfo(type(x)).eps

        if x != 0.0:

            ITmax = 100

            ap = a
            summ = 1.0 / a
            dell = summ

            for i in range(1, ITmax + 1):

                ap += 1
                dell = dell * x / ap
                summ += dell

                if abs(dell) < abs(summ) * eps:
                    break

        elif x == 0.0:
            return 0.0
        else:
            raise ValueError

        return summ * exp(-x + a * log(x) - self.dgammln(a))

    def dgcf(self, a, x):

        eps = np.finfo(type(x)).eps  # Equivalent to EPSILON(x)
        fpmin = np.finfo(type(x)).tiny / eps  # Equivalent to TINY(x) / EPS

        itmax = 100

        if x != 0.0:

            b = x + 1 - a
            c = 1 / fpmin
            d = 1.0 / b
            h = d

            for i in range(1, itmax + 1):
                an = -i * (i - a)
                b = b + 2
                d = an * d + b
                if abs(d) < fpmin:
                    d = fpmin

                c = b + an / c
                if abs(c) < fpmin:
                    c = fpmin

                d = 1 / d
                dell = d * c
                h = h * dell
                if abs(dell - 1) <= eps:
                    break
        elif x == 0:
            return 1.0

        else:
            raise ValueError

        return exp(-x + a * log(x) - self.dgammln(a)) * h


class heat_flux(heat_node, vapor_flux):

    def __init__(self, left_node, right_node, top_layer, advection=False):
        """
        Class to calculate the heat flux based on a given temprature gradient and hydraulic head (psi in m) gradient.
        Based on Haverd&Cuntz (2005) Eq. A.12
        """

        heat_node.__init__(self, top_layer=top_layer)
        vapor_flux.__init__(self, left_node=left_node, right_node=right_node, top_layer=top_layer)

        self.advection = advection

    def qh(self):

        Tzero = 273.16000366210938
        ln, rn = self.left_node, self.right_node

        q = self.keff() * (ln.T - rn.T) + \
            (self.phiv(ln) - self.phiv(rn)) * self.lambdav(ln) * 1000 / self.rdz() * \
            (((ln.T / Tzero) ** 1.88 + (rn.T / Tzero) ** 1.88) / 2) * \
            ((self.cvsat_0(ln) + self.cvsat_0(rn)) / 2)

        return q

    def qhya(self):

        try:
            ln, rn = self.left_node, self.right_node
            S = ln.S

            if S < 1:  # unsaturated

                Tzero = 273.16000366210938

                phivs = self.phivS(ln)
                lambdav = self.lambdav(ln)
                dz = self.rdz()

                q = 1 / dz * lambdav * 1000 * phivs * \
                    (((ln.T / Tzero) ** 1.88 + (rn.T / Tzero) ** 1.88) / 2) * \
                    ((self.cvsat_0(ln) + self.cvsat_0(rn)) / 2)

            elif S >= 1:  # saturated

                q = 0
            else:
                raise ValueError

            return q

        except ValueError:
            raise NotImplementedError

    def qhyb(self):

        try:
            ln, rn = self.left_node, self.right_node
            S = ln.S

            if S < 1:  # unsaturated

                Tzero = 273.16000366210938

                phivs = self.phivS(rn)
                lambdav = self.lambdav(ln)
                dz = self.rdz()

                q = -1 / dz * lambdav * 1000 * phivs * \
                    (((ln.T / Tzero) ** 1.88 + (rn.T / Tzero) ** 1.88) / 2) * \
                    ((self.cvsat_0(ln) + self.cvsat_0(rn)) / 2)

            elif S >= 1:  # saturated

                q = 0
            else:
                raise ValueError

            return q

        except ValueError:
            raise NotImplementedError

    def qhTa(self):
        return self.keff()

    def qhTb(self):
        return - self.keff()

    def keff(self):

        ln, rn = self.left_node, self.right_node

        k = 2 * (self.kth(ln) * self.kth(rn)) / (self.kth(ln) * ln.thickness + self.kth(rn) * rn.thickness)

        return k


class heat_solver(object):

    def __int__(self):
        pass

    def csoil(self, node):

        rhow = 1000  # density of water
        cswat = 4.218e3  # sp. heat capacity of water []
        csice = 4.218e3  # sp. heat capacity of ice
        thetaI = 0  # ice content

        return node.css * node.rho + rhow * cswat * node.theta + rhow * csice * thetaI

    def csoileff(self, node):

        rhow = 1000  # density of water
        lambdaf = 335000  # latent heat of fusion [J kg-1]
        c_soil = self.csoil(node)
        dthetadT = 0  # no ice

        return c_soil + rhow * lambdaf * dthetadT

    def iice(self, node):

        return 0  # no ice

    def delta_ddh(self, node, dt):

        imelt = 0
        dthetadT = 0  # no ice
        rhow = 1000  # density of water
        cswat = 4.218e3  # sp. heat capacity of water []
        csice = 4.218e3  # sp. heat capacity of ice

        r_dt = 1 / dt
        dx = node.thickness
        csoil_eff = self.csoileff(node)
        c_soil = self.csoil(node)

        d = csoil_eff * dx * r_dt * (1 - imelt) + c_soil * dx * r_dt * imelt + \
            (cswat - csice) * dx * dthetadT * rhow * node.T * r_dt * (1 - imelt) * self.iice(node)

        return d

    def delta_cc(self, node, dt):

        if node.S < 1:  # unsaturated

            return (node.theta_sat - node.theta_0) * node.thickness * 1 / dt

        elif node.S >= 1:  # saturated
            return 0
        else:
            raise ValueError

    def delta_cch(self, node, dt):

        imelt = 0
        rhow = 1000  # density of water
        lambdaf = 335000.   # latent heat of fusion(J kg - 1)
        r_dt = 1 / dt
        dx = node.thickness
        thre = node.theta_sat - node.theta_0

        if node.S < 1:  # unsaturated
            isat = 0
        elif node.S >= 1:  # Saturated
            isat = 1
        else:
            return ValueError

        d = (1 - imelt) * self.iice(node) * (1 - isat) * rhow * lambdaf * thre * dx * r_dt

        return d

    def delta_ggh(self, node, dt):

        imelt = 0

        rhow = 1000  # density of water
        cswat = 4.218e3  # sp. heat capacity of water []
        csice = 4.218e3  # sp. heat capacity of ice
        lambdaf = 335000.  # latent heat of fusion(J kg - 1)
        thetai = 0  # theta ice

        r_dt = 1 / dt
        dx = node.thickness

        d = rhow * (lambdaf + (cswat - csice) * node.T) * thetai * dx * r_dt * imelt

        return d

    def solve_sparse(self, x, X):

        n = len(x[2])

        allmat = self.build_allmat(x, X)

        lS, lT, cS, cT = self.conditioning(allmat, n)

        lc = [lS, lT, cS, cT]

        A, B, C, d = self.block_matrix(x, X, n, lc)

        x = self.generic_thomas(n, A, B, C, d)

        dy = x[:, 0] * cS
        dT = x[:, 1] * cT

        return dy, dT

    def build_allmat(self, x, X):

        aa, bb, cc, dd, ee, ff = x[0], x[1], x[2], x[3], x[4], x[5]
        aah, bbh, cch, ddh, eeh, ffh = X[0], X[1], X[2], X[3], X[4], X[5]

        n = len(cc)

        """Constructs the allmat matrix without mimicking Fortran indexing."""

        # Initialize 2D blocks
        A = np.zeros((n, 2, 2))
        B = np.zeros((n, 2, 2))
        C = np.zeros((n - 1, 2, 2))  # Last row doesn't need C

        # Fill matrices
        A[1:, 0, 0] = aa
        A[1:, 0, 1] = bb
        A[1:, 1, 0] = aah
        A[1:, 1, 1] = bbh

        B[:, 0, 0] = cc
        B[:, 0, 1] = dd
        B[:, 1, 0] = cch
        B[:, 1, 1] = ddh

        C[:, 0, 0] = ee
        C[:, 0, 1] = ff
        C[:, 1, 0] = eeh
        C[:, 1, 1] = ffh

        # Build block tridiagonal matrix
        allmat = np.block([
            [np.diag(B[:, 0, 0]), np.diag(B[:, 0, 1])],
            [np.diag(B[:, 1, 0]), np.diag(B[:, 1, 1])]
        ])

        for i in range(n - 1):
            allmat[2 * i + 2:2 * i + 4, 2 * i:2 * i + 2] = C[i]  # Lower diagonal
            allmat[2 * i:2 * i + 2, 2 * i + 2:2 * i + 4] = A[i + 1]  # Upper diagonal

        return allmat

    def conditioning(self, allmat, n):

        eps = np.finfo(allmat.dtype).eps * 1000

        # cond == 3, rows and columns:

        lST = np.max(np.abs(allmat), axis=1)
        lST[lST < eps] = 1.0
        lST = 1.0 / lST
        lS = lST[0:2 * n:2]
        lT = lST[1:2 * n:2]
        allmat_scaled = lST[:, np.newaxis] * allmat
        cST = np.max(np.abs(allmat_scaled), axis=0)
        cST[cST < eps] = 1.0
        cST = 1.0 / cST
        cS = cST[0:2 * n:2]
        cT = cST[1:2 * n:2]

        return lS, lT, cS, cT

    def block_matrix(self, x, X, n, lc):

        aa, bb, cc, dd, ee, ff, gg = x[0], x[1], x[2], x[3], x[4], x[5], x[6]
        aah, bbh, cch, ddh, eeh, ffh, ggh = X[0], X[1], X[2], X[3], X[4], X[5], X[6]

        lS, lT, cS, cT = lc[0], lc[1], lc[2], lc[3]

        # Assemble the 2x2 block matrices A, B, and C.
        # A, B, C have shape (n, 2, 2)
        A_block = np.zeros((n, 2, 2))
        B_block = np.zeros((n, 2, 2))
        C_block = np.zeros((n, 2, 2))

        # For the subdiagonal (A): For i = 1 to n-1 (Python index 1 to n-1)
        for i in range(1, n):
            A_block[i, 0, 0] = aa[i - 1] * lS[i] * cS[i - 1]
            A_block[i, 0, 1] = bb[i - 1] * lS[i] * cT[i - 1]
            A_block[i, 1, 0] = aah[i - 1] * lT[i] * cS[i - 1]
            A_block[i, 1, 1] = bbh[i - 1] * lT[i] * cT[i - 1]
        # The first block row of A remains zeros.

        # Main diagonal (B): For i = 0 to n-1
        for i in range(n):
            B_block[i, 0, 0] = cc[i] * lS[i] * cS[i]
            B_block[i, 0, 1] = dd[i] * lS[i] * cT[i]
            B_block[i, 1, 0] = cch[i] * lT[i] * cS[i]
            B_block[i, 1, 1] = ddh[i] * lT[i] * cT[i]

        # Superdiagonal (C): For i = 0 to n-2
        for i in range(n - 1):
            C_block[i, 0, 0] = ee[i] * lS[i] * cS[i + 1]
            C_block[i, 0, 1] = ff[i] * lS[i] * cT[i + 1]
            C_block[i, 1, 0] = eeh[i] * lT[i] * cS[i + 1]
            C_block[i, 1, 1] = ffh[i] * lT[i] * cT[i + 1]
        # The last block row of C remains zeros.

        # Assemble right-hand side vector d with shape (n,2)
        d_rhs = np.zeros((n, 2))
        for i in range(n):
            d_rhs[i, 0] = gg[i] * lS[i]
            d_rhs[i, 1] = ggh[i] * lT[i]

        return A_block, B_block, C_block, d_rhs

    def generic_thomas(self, n, A, B, C, d):

        # Make copies so we do not alter originals:
        B_mod = B.copy()
        d_mod = d.copy()

        # Forward elimination:
        for i in range(1, n):
            # Solve B_mod[i-1] * y = d_mod[i-1] and invert B_mod[i-1]
            invB_prev = np.linalg.inv(B_mod[i - 1])
            # Multiplier M: A[i] * inv(B[i-1])
            M = A[i] @ invB_prev
            B_mod[i] = B_mod[i] - M @ C[i - 1]
            d_mod[i] = d_mod[i] - M @ d_mod[i - 1]

        # Back substitution:
        x = np.zeros_like(d_mod)
        x[-1] = np.linalg.solve(B_mod[-1], d_mod[-1])
        for i in range(n - 2, -1, -1):
            x[i] = np.linalg.solve(B_mod[i], d_mod[i] - C[i] @ x[i + 1])

        return x

    def solve_tridiagonal_dT(aa, bb, cc, dd, ee, ff, gg):
        """
        Solves a tridiagonal system for dT only:
            lower_diag[i] * dT[i-1] + diag[i] * dT[i] + upper_diag[i] * dT[i+1] = ggh[i]

        This merges the humidity-only components.

        Parameters:
            aah, bbh: lower-diagonal components (length n-1)
            cch, ddh: main-diagonal components (length n)
            eeh, ffh: upper-diagonal components (length n-1)
            ggh: right-hand side vector (length n)

        Returns:
            dT: solution vector (length n)
        """

        # Convert all input to numpy arrays (in case they're passed as lists)
        aah = np.asarray(aa)
        bbh = np.asarray(bb)
        cch = np.asarray(cc)
        ddh = np.asarray(dd)
        eeh = np.asarray(ee)
        ffh = np.asarray(ff)
        ggh = np.asarray(gg)

        n = len(cc)
        if not (len(dd) == len(cc) == len(gg) and
                len(aa) == len(bb) == n - 1 and
                len(ee) == len(ff) == n - 1):
            raise ValueError("Input array dimensions do not match expected sizes.")

        # Construct diagonals
        a = aa + bb  # lower diagonal (length n-1)
        b = cc + dd  # main diagonal (length n)
        c = ee + ff  # upper diagonal (length n-1)
        d = gg.copy()  # RHS (length n)

        # Forward elimination
        for i in range(1, n):
            m = a[i - 1] / b[i - 1]
            b[i] -= m * c[i - 1]
            d[i] -= m * d[i - 1]

        # Back substitution
        dT = np.zeros(n)
        dT[-1] = d[-1] / b[-1]
        for i in range(n - 2, -1, -1):
            dT[i] = (d[i] - c[i] * dT[i + 1]) / b[i]

        return dT









