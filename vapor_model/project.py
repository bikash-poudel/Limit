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

    def run(self, dt=1, dt_max=60, epsilon=1e-4, model='proposed',
            water_flux=True, heat_flux=True, vapor_flux=True):

        if model == 'proposed':

            if water_flux:
                self.water_equation(delta_t=dt, epsilon=epsilon)

            if heat_flux:
                dt = self.heat_equation(delta_t=dt, epsilon=epsilon)

            if vapor_flux:
                dt = self.vapor_equation(delta_t=dt, dt_max=dt_max, epsilon=epsilon)

        elif model == "coupled":

            dt = self.coupled_equation(delta_t=dt, dt_max=dt_max, epsilon=epsilon)

        else:
            raise NotImplementedError

        return dt

    def water_equation(self, delta_t, epsilon=1e-4):

        ph = self.get_layers()[0].phi

        self.cells[0].update_boundary_head()
        self.cells[0].update_boundary_T()

        headt = self.solve(delta_t, equation='w')
        headt[headt > ph] = ph

        water_fail = 1
        p = 0
        while p < 5:
            p += 1
            head = headt
            self.cells[0].update_head(head)
            self.cells[0].update_boundary_head()
            self.cells[0].update_boundary_T()

            headt = self.solve(delta_t=delta_t, equation='w')
            headt[headt > ph] = ph

            self.cells[0].update_head(headt)

            max1, max2 = np.max(np.abs(headt - head) / np.max(np.abs(head))), 0.00

            if np.maximum(max1, max2) < epsilon:
                p, head, water_fail = 6, headt, 0

            if p == 5:
                delta_t, water_fail = delta_t / 2, 1

    def heat_equation(self, delta_t, epsilon=1e-4):

        #self.cells[0].update_boundary_T()
        tempt = self.solve(delta_t, equation='h')

        heat_fail = 1
        p = 0
        while p < 5:
            p += 1
            temp = tempt
            self.cells[0].update_temperature(temp)
            #self.cells[0].update_boundary_T()

            tempt = self.solve(delta_t=delta_t, equation='h')
            self.cells[0].update_temperature(tempt)

            max1, max2 = 0.0, np.max(np.abs(tempt - temp) / np.max(np.abs(temp)))

            if np.maximum(max1, max2) < epsilon:
                p, temp, heat_fail = 6, tempt, 0

            if p == 5:
                delta_t, heat_fail = delta_t / 2, 1

            return delta_t

    def vapor_equation(self, delta_t, dt_max, epsilon=1e-4):

        ph = self.get_layers()[0].phi  # param phi_e

        headt = np.array([l.head for l in self.cells[0].layers])
        tempt = np.array([l.T for l in self.cells[0].layers])

        # self.cells[0].update_boundary_head()
        # self.cells[0].update_boundary_T()
        head, T = self.solve(delta_t, equation='v')

        headss = np.clip(headt + np.array(head), None, ph)  # Clip headss to ph
        tempss = tempt + np.array(T)

        p = 0
        while p < 5:
            p += 1
            headp, tempp = headss, tempss
            heads = 0.5 * (headss + headt)
            temps = 0.5 * (tempss + tempt)

            self.cells[0].update_head(heads), self.cells[0].update_temperature(temps)
            #self.cells[0].update_boundary_head()
            #self.cells[0].update_boundary_T()

            head, temp = self.solve(delta_t=delta_t, equation='v')

            headss = np.clip(headt + np.array(head), None, ph)  # Clip headss to ph
            tempss = tempt + np.array(temp)

            self.cells[0].update_head(headss), self.cells[0].update_temperature(tempss)

            max1 = np.max(np.abs(headss - headp) / np.max(np.abs(headp)))
            max2 = np.max(np.abs(tempss - tempp) / np.max(np.abs(tempp)))

            if np.maximum(max1, max2) < epsilon:
                p, head, temp = 6, headss, tempss
                delta_t = np.minimum(dt_max, delta_t * 1.5)

            if p == 5:
                delta_t = delta_t / 2

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

            a = A.tocsr()  # .astype(np.float64)
            b = B  # .astype(np.float64)

        except ValueError:
            raise NotImplementedError

        return spsolve(a, b)

    def water_flux_matrix(self, mat_A, mat_B, storages, dt):

        s_indices = {node: idx for idx, node in enumerate(storages)}  # Precompute indices

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

            ql = []
            for c in self.cells[0].water_connections:

                if isinstance(c, fluxes.water_flux):

                    index_l = s_indices[c.left_node]  # index of left storage node
                    index_r = s_indices[c.right_node]  # index of right storage node

                    r = self.r(c, dt)
                    hy_cond_pot, q_tmp_water = c.hy_conduct_potential() * r, c.q_tmp_water() * r

                    mat_A[index_l, index_r] -= hy_cond_pot
                    mat_A[index_l, index_l] += hy_cond_pot
                    mat_A[index_r, index_r] += hy_cond_pot

                    if not c.right_node == storages[-1]:
                        mat_A[index_r, index_l] -= hy_cond_pot

                    if not c.left_node == storages[0]:
                        mat_B[index_l] -= q_tmp_water

                    if not c.right_node == storages[-1]:
                        mat_B[index_r] += q_tmp_water

                    dz = c.left_node.distance(c.right_node)
                    ql.append(c.q_water() / dz)

                else:
                    raise NotImplementedError

            self.cells[0].liquid_fluxes = ql

            mat_A[-1, -1], mat_B[-1] = 1, self.cells[0].layers[-1].head

        except ValueError:
            raise NotImplementedError

    def heat_flux_matrix(self, mat_A, mat_B, storages, dt):

        s_indices = {node: idx for idx, node in enumerate(storages)}  # Precompute indices

        try:

            for c in self.cells[0].heat_connections:

                index_l = s_indices[c.left_node]  # Index of left storage node
                index_r = s_indices[c.right_node]  # Index of right storage node

                r = self.r(c, dt)
                th_con, th_con_pot = c.th_conductivity() * r, c.th_hy_conductivity_T0() * r

                if not index_l == 0:  # if left node is the boundary  storages(layer)
                    mat_A[index_l, index_r] -= th_con

                if not c.right_node == storages[-1]:  # if right node is not the boundary layer (storage)
                    mat_A[index_r, index_l] -= th_con

                mat_A[index_l, index_l] += th_con
                mat_A[index_r, index_r] += th_con

                if c.right_node.head > c.left_node.head:
                    if not index_l == 0:
                        mat_A[index_l, index_r] -= c.th_hy_conductivity_correction() * r
                    mat_A[index_r, index_r] += c.th_hy_conductivity_correction() * r
                else:
                    if not c.right_node == storages[-1]:
                        mat_A[index_r, index_l] += c.th_hy_conductivity_correction() * r
                    mat_A[index_l, index_l] -= c.th_hy_conductivity_correction() * r

                mat_B[index_l] -= th_con_pot
                mat_B[index_r] += th_con_pot

            l_left, l_right = self.cells[0].layers[0], self.cells[0].layers[-1]

            mat_A[0, 0], mat_B[0] = 1, l_left.T
            mat_A[-1, -1], mat_B[-1] = 1, l_right.T

        except ValueError:
            raise NotImplementedError

    def vapor_solve(self, storages, dt):

        try:

            vapor_connections = self.cells[0].vapor_connections
            q_vapor = np.array([self.r(c, dt) * c.q_vapor() for c in vapor_connections])
            qv_latent = np.array([s.latent for s in storages])

            distances = np.array([c.left_node.distance(c.right_node) for c in vapor_connections])
            qv = q_vapor / distances

            self.cells[0].vapor_fluxes = qv

            h, T = [], []
            # left boundary
            qv_right, qv_left = q_vapor[0],0

            if qv_right > 0:
                qv_latent_right = qv_latent[1] * qv_right
            elif qv_right < 0:
                qv_latent_right = qv_latent[0] * qv_right
            else:
                qv_latent_right = 0

            qv_latent_left = qv_latent[0] * qv_left

            B = np.array([qv_right - qv_left, qv_latent_right - qv_latent_left])
            A = np.array([[storages[0].H1, storages[0].H2], [storages[0].T2, storages[0].T1]])
            mat_A = csr_matrix(A)

            kk = spsolve(mat_A, B)
            h.append(kk[0])
            T.append(kk[1])

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

                B = np.array([qv_right - qv_left, qv_latent_right - qv_latent_left])
                A = np.array([[storage.H1, storage.H2], [storage.T2, storage.T1]])
                mat_A = csr_matrix(A)

                kk = spsolve(mat_A, B)
                h.append(kk[0])
                T.append(kk[1])

            # right boundary

            qv_right, qv_left = 0, q_vapor[-1]

            qv_latent_right = 0
            if qv_left > 0:
                qv_latent_left = qv_latent[-1] * qv_right
            elif qv_left < 0:
                qv_latent_left = qv_latent[-2] * qv_right
            else:
                qv_latent_left = 0

            B = np.array([qv_right - qv_left, qv_latent_right - qv_latent_left])
            A = np.array([[storages[-1].H1, storages[-1].H2], [storages[-1].T2, storages[-1].T1]])
            mat_A = csr_matrix(A)

            kk = spsolve(mat_A, B)
            h.append(kk[0])
            T.append(kk[1])

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

    def coupled_equation(self, delta_t, dt_max, epsilon=1e-4):

        ph = self.get_layers()[0].phi  # param phi_e

        storages = list(self.get_layers())

        headt = np.array([l.head for l in self.cells[0].layers])
        tempt = np.array([l.T for l in self.cells[0].layers])

        # self.cells[0].update_boundary_head()
        self.cells[0].update_boundary_T()
        head, T = self.coupled_solver(storages, delta_t)

        headss = headt + np.array(head)
        tempss = tempt + np.array(T)

        headss[-1] = head[-1]
        headss[headss > ph] = ph
        tempss[0], tempss[-1] = T[0], T[-1]

        p = 0
        while p < 5:

            p += 1
            headpp, tempp = headss, tempss
            heads = 0.5 * (headss + headt)
            temps = 0.5 * (tempss + tempt)

            self.cells[0].update_head(heads), self.cells[0].update_temperature(temps)
            # self.cells[0].update_boundary_head()
            self.cells[0].update_boundary_T()

            head, temp = self.coupled_solver(storages, delta_t)

            headss = headt + np.array(head)
            tempss = tempt + np.array(temp)

            headss[-1] = head[-1]
            headss[headss > ph] = ph
            tempss[0], tempss[-1] = T[0], T[-1]

            max1 = np.max(np.abs(headss - headpp) / np.max(np.abs(headpp)))
            max2 = np.max(np.abs(tempss - tempp) / np.max(np.abs(tempp)))

            if np.maximum(max1, max2) < epsilon:
                p, head, temp = 6, headss, tempss
                delta_t = np.minimum(dt_max, delta_t * 1.5)

            if p == 5:
                delta_t = delta_t / 2

        self.cells[0].update_head(headss), self.cells[0].update_temperature(tempss)

        self.cells[0].liquid_fluxes = [c.q_water() / c.left_node.distance(c.right_node)
                                       for c in self.cells[0].water_connections]
        self.cells[0].vapor_fluxes = [c.q_vapor() / c.left_node.distance(c.right_node)
                                      for c in self.cells[0].vapor_connections]

        return delta_t

    def coupled_solver(self, storages, dt):

        h, T = [], []

        ql_latent = [s.latent_liquid for s in storages]
        qv_latent = [s.latent for s in storages]

        # upper Boundary
        s0 = storages[0]
        left_c, right_c = s0.connections_to_left, s0.connections_to_right
        q_evap = s0.connections_to_boundaries[0].q_evap

        ql_left, ql_right = 0, right_c[0].q_water() * self.r(right_c[0], dt)
        qv_left, qv_right = dt / s0.thickness * q_evap, right_c[2].q_vapor() * self.r(right_c[2], dt)

        # solving
        A = [[s0.H1, 0], [0, 1]]
        B = [ql_right + qv_right - ql_left - qv_left, s0.T]

        mat_A = csr_matrix(A)
        kk = spsolve(mat_A, B)
        h.append(kk[0]), T.append(kk[1])

        for i, s in enumerate(storages[1:-1]):

            left_c, right_c = s.connections_to_left, s.connections_to_right

            ql_left = left_c[0].q_water() * self.r(left_c[0], dt)
            ql_right = right_c[0].q_water() * self.r(right_c[0], dt)

            qv_left = left_c[2].q_vapor() * self.r(left_c[2], dt)
            qv_right = right_c[2].q_vapor() * self.r(right_c[2], dt)

            q_cond_left = left_c[1].q_cond() * self.r(left_c[1], dt)
            q_cond_right = right_c[1].q_cond() * self.r(right_c[1], dt)

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

            # solving
            A = [[s.H1, s.H2],
                 [s.T2, s.T1]]
            B = [ql_right + qv_right - ql_left - qv_left,
                 q_cond_right - q_cond_left + ql_latent_right - ql_latent_left + qv_latent_right - qv_latent_left]

            mat_A = csr_matrix(A)
            kk = spsolve(mat_A, B)
            h.append(kk[0]), T.append(kk[1])

        # Lower Boundary
        sn = storages[-1]
        left_c, right_c = sn.connections_to_left, sn.connections_to_right

        ql_left, ql_right = left_c[0].q_water() * self.r(left_c[0], dt), 0
        qv_left, qv_right = left_c[2].q_vapor() * self.r(left_c[2], dt), 0

        A = [[sn.H1, 0], [0, 1]]
        B = [ql_right + qv_right - ql_left - qv_left, sn.T]

        mat_A = csr_matrix(A)
        kk = spsolve(mat_A, B)
        h.append(kk[0]), T.append(kk[1])

        return h, T

    def heat_solver(self, storages, dt):

        ql_latent = [s.latent_liquid for s in storages]
        qv_latent = [s.latent for s in storages]

        for i, s in enumerate(storages[1:-1]):

            left_c, right_c = s.connections_to_left, s.connections_to_right

            ql_left = left_c[0].q_water() * self.r(left_c[0], dt)
            ql_right = right_c[0].q_water() * self.r(right_c[0], dt)

            qv_left = left_c[2].q_vapor() * self.r(left_c[2], dt)
            qv_right = right_c[2].q_vapor() * self.r(right_c[2], dt)

            q_cond_left = left_c[1].q_cond() * self.r(left_c[1], dt)
            q_cond_right = right_c[1].q_cond() * self.r(right_c[1], dt)

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

            # solving
            A = [[s.H1, s.H2],
                 [s.T2, s.T1]]
            B = [ql_right + qv_right - ql_left - qv_left,
                 q_cond_right - q_cond_left + ql_latent_right - ql_latent_left + qv_latent_right - qv_latent_left]

            mat_A = csr_matrix(A)
            kk = spsolve(mat_A, B)

            T = s.T + (q_cond_right - q_cond_left
                       + ql_latent_right - ql_latent_left
                       + qv_latent_right - qv_latent_left) / s.T1

    def r(self, connection, dt):

        dist = connection.left_node.distance(connection.right_node)
        r = dt / dist / dist

        return r



