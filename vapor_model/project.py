'''
Created on 11.12.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

from . import storage
from . import cell
from . import fluxes

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


class project(object):

    def __init__(self):

        self.__cells = []

    @property
    def cells(self):
        return self.__cells

    def get_layers(self):

        return self.cells[0].layers

    def get_flux_connections(self):

        return self.cells[0].flux_connections

    def new_cell(self, atmosphere):
        """
        Creates a new cell and registers it in the project

        Parameters
        ----------
        atmosphere
        """
        # i_point = storages.Point(x, y, z)
        c = cell.cell(atmosphere=atmosphere)
        self.__cells.append(c)

    def run(self, dt=1, dt_max=60, epsilon=1e-4):

        self.water_equation(delta_t=dt, epsilon=epsilon)

        self.heat_equation(delta_t=dt, epsilon=epsilon)

        dt = self.vapor_equation(delta_t=dt, dt_max=dt_max, epsilon=epsilon)

        return dt

    def water_equation(self, delta_t, epsilon=1e-4):

        self.update_boundaries()

        headt = self.solve(delta_t, equation='w')
        headt[headt > -0.13] = -0.13

        water_fail = 1
        p = 0
        head = None
        while p < 5:
            p += 1
            head = headt
            self.cells[0].update_head(head)
            self.update_boundaries()

            headt = self.solve(delta_t=delta_t, equation='w')
            headt[headt > -0.13] = -0.13

            max1, max2 = np.max(np.abs(headt - head) / np.max(np.abs(head))), 0.00

            if np.maximum(max1, max2) < epsilon:
                p, head, water_fail = 6, headt, 0

            if p == 5:
                delta_t, water_fail = delta_t / 2, 1

        self.cells[0].update_head(head)

    def heat_equation(self, delta_t, epsilon=1e-4):

        self.update_boundaries()
        tempt = self.solve(delta_t, equation='h')

        heat_fail = 1
        p = 0
        temp = None
        while p < 5:
            p += 1
            temp = tempt
            self.cells[0].update_temperature(temp), self.update_boundaries()

            tempt = self.solve(delta_t=delta_t, equation='h')

            max1, max2 = 0.0, np.max(np.abs(tempt - temp) / np.max(np.abs(temp)))

            if np.maximum(max1, max2) < epsilon:
                p, temp, heat_fail = 6, tempt, 0

            if p == 5:
                delta_t, heat_fail = delta_t / 2, 1

        self.cells[0].update_temperature(temp)

    def vapor_equation(self, delta_t, dt_max, epsilon=1e-4):

        headt = np.array([l.head for l in self.cells[0].layers])
        tempt = np.array([l.T for l in self.cells[0].layers])

        self.update_boundaries()
        head, T = self.solve(delta_t, equation='v')

        headss = headt + head
        tempss = T + tempt
        headss[headss > -0.13] = -0.13

        p = 0
        while p < 5:
            p += 1
            headp, tempp = headss, tempss
            heads = 0.5 * (headss + headt)
            temps = 0.5 * (tempss + tempt)

            self.cells[0].update_head(heads), self.cells[0].update_temperature(temps)
            self.update_boundaries()

            head, temp = self.solve(delta_t=delta_t, equation='v')

            headss = headt + head
            headss[headss > -0.13] = -0.13

            tempss = temp + tempt

            max1 = np.max(np.abs(headss - headp) / np.max(np.abs(headp)))
            max2 = np.max(np.abs(tempss - tempp) / np.max(np.abs(tempp)))

            if np.maximum(max1, max2) < epsilon:
                p, head, temp = 6, headss, tempss
                delta_t = np.minimum(dt_max, delta_t * 1.5)

            if p == 5:
                delta_t = delta_t / 2

        self.cells[0].update_head(head), self.cells[0].update_temperature(temp)

        return delta_t

    def solve(self, delta_t, equation='w'):

        storages = list(self.get_layers())

        A = lil_matrix((len(storages), len(storages)))  # create a n x n matrix
        B = np.zeros(len(storages))

        try:
            if equation == 'w':
                for storage in storages:
                    s_index = storages.index(storage)

                    A[s_index, s_index] += storage.w_capacity
                    B[s_index] += storage.w_capacity * storage.head
                self.water_flux_matrix(A, B, storages, dt=delta_t)

            elif equation == 'h':
                for storage in storages:
                    s_index = storages.index(storage)

                    A[s_index, s_index] += storage.cs
                    B[s_index] += storage.cs * storage.T
                self.heat_flux_matrix(A, B, storages, dt=delta_t)

            elif equation == 'v':

                return self.vapor_solve(storages, dt=delta_t)

            else:
                raise NotImplementedError

            a = A.tocsr()
            b = B

        except ValueError:
            raise NotImplementedError

        return spsolve(a, b)

    def water_flux_matrix(self, mat_A, mat_B, storages, dt):

        try:
            c_evap = self.cells[0].connection_evap
            # check which storage node is connected to boundary
            if isinstance(c_evap.left_node, storage.soil_layer):
                s = c_evap.left_node
            elif isinstance(c_evap.right_node, storage.soil_layer):
                s = c_evap.right_node
            else:
                raise ValueError("Boundary connection should have one node connected to soil_layer")

            index_s = storages.index(s)

            mat_B[index_s] -= c_evap.q_evap * dt / s.thickness

            for c in self.cells[0].water_connections:

                if isinstance(c, fluxes.water_flux):

                    index_l = storages.index(c.left_node)  # index of left storage node
                    index_r = storages.index(c.right_node)  # index of right storage node

                    r = self.r(c, dt)
                    hy_cond_pot, q_tmp_water = c.hy_conduct_potential() * r, c.q_tmp_water() * r

                    mat_A[index_l, index_r] -= hy_cond_pot
                    mat_A[index_r, index_l] -= hy_cond_pot
                    mat_A[index_l, index_l] += hy_cond_pot
                    mat_A[index_r, index_r] += hy_cond_pot

                    mat_B[index_l] += q_tmp_water
                    mat_B[index_r] -= q_tmp_water

                else:
                    raise NotImplementedError

            mat_A[0, 0], mat_B[0] = 1, self.cells[0].layers[0].head

        except ValueError:
            raise NotImplementedError

    def heat_flux_matrix(self, mat_A, mat_B, storages, dt):

        try:

            for c in self.cells[0].heat_connections:

                index_l = storages.index(c.left_node)  # index of left storage node
                index_r = storages.index(c.right_node)  # index of right storage node

                r = self.r(c, dt)
                th_con, th_con_pot = c.th_conductivity() * r, c.th_hy_conductivity_T0() * r

                mat_A[index_l, index_r] -= th_con
                mat_A[index_r, index_l] -= th_con
                mat_A[index_l, index_l] += th_con
                mat_A[index_r, index_r] += th_con

                if c.right_node.head > c.left_node.head:
                    mat_A[index_l, index_r] -= c.th_hy_conductivity_correction() * r
                    mat_A[index_l, index_r] += c.th_hy_conductivity_correction() * r
                else:
                    mat_A[index_r, index_l] += c.th_hy_conductivity_correction() * r
                    mat_A[index_l, index_r] -= c.th_hy_conductivity_correction() * r

                mat_B[index_l] += th_con_pot
                mat_B[index_r] -= th_con_pot

            l_left, l_right = self.cells[0].layers[0], self.cells[0].layers[-1]

            mat_A[0, 0], mat_B[0] = 1, l_left.T
            mat_A[-1, -1], mat_B[-1] = 1, l_right.T

        except ValueError:
            raise NotImplementedError

    def vapor_solve(self, storages, dt):

        try:

            q_vapor = [self.r(c, dt) * c.q_vapor() for c in self.cells[0].vapor_connections]
            qv_latent = [s.latent for s in storages]
            h, T = [], []

            # left boundary
            qv_right, qv_left = q_vapor[0], 0

            if qv_right > 0:
                qv_latent_right = qv_latent[1] * qv_right
            elif qv_right < 0:
                qv_latent_right = qv_latent[0] * qv_right
            else:
                qv_latent_right = 0

            qv_latent_left = 0

            B = [qv_right - qv_left,
                 qv_latent_right - qv_latent_left]
            A = np.array([[storages[0].H1, storages[0].H2],
                          [storages[0].T2, storages[0].T1]])

            mat_A = csr_matrix(A)

            kk = spsolve(mat_A, B)
            h.append(kk[0]), T.append(kk[1])

            for i, storage in enumerate(storages[1:-1]):

                qv_right, qv_left = q_vapor[i + 1], q_vapor[i]

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

                B = [qv_right - qv_left, qv_latent_right - qv_latent_left]
                A = np.array([[storage.H1, storage.H2], [storage.T2, storage.T1]])

                mat_A = csr_matrix(A)
                kk = spsolve(mat_A, B)
                h.append(kk[0]), T.append(kk[1])

            # right boundary

            qv_right, qv_left = 0, q_vapor[-1]

            qv_latent_right = 0
            if qv_left > 0:
                qv_latent_left = qv_latent[-1] * qv_right
            elif qv_left < 0:
                qv_latent_left = qv_latent[-2] * qv_right
            else:
                qv_latent_left = 0

            B = [qv_right - qv_left, qv_latent_right - qv_latent_left]
            A = np.array([[storages[0].H1, storages[0].H2], [storages[0].T2, storages[0].T1]])

            mat_A = csr_matrix(A)
            kk = spsolve(mat_A, B)
            h.append(kk[0]), T.append(kk[1])

        except ValueError:
            raise NotImplementedError

        return h, T

    def update_boundaries(self, head_boundary="True", temp_boundary="True", evaporation="True"):

        try:
            c = self.cells[0]
            if evaporation:
                c.update_evaporation()
            if head_boundary:
                c.update_boundary_head()
            if temp_boundary:
                c.update_boundary_T()

        except ValueError:
            raise NotImplementedError

    def r(self, connection, dt):

        dist = connection.left_node.distance(connection.right_node)
        r = dt / dist / dist

        return r

