'''
Created on 20.06.2024

@author: poudel-b
'''
# -*- coding: utf-8 -*-
from __future__ import division

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

import iso_cell
import iso_fluxes
import iso_storages


def storage_fluxes(storage, Isotopologue, **kwargs):
    try:
        assert isinstance(storage, iso_storages.iso_storage)

        q = 0  # sum of fluxes
        for c in storage.connections:

            flux = c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

            if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)) and c.right_node is storage:

                q += flux

            else:
                q -= flux

    except Exception as err:
        raise NotImplementedError("An unexpected error occurred.") from err

    return q


def right_fluxes(storage, Isotopologue, **kwargs):
    try:
        assert isinstance(storage, iso_storages.iso_storage)
        s_connections = storage.connections_to_iso_storages

        q = 0  # sum of fluxes
        for c in s_connections:

            if c.left_node is storage:
                if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):
                    q -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                else:
                    q += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

    except Exception as err:
        raise NotImplementedError("An unexpected error occurred.") from err

    return q


def left_fluxes(storage, Isotopologue, **kwargs):
    try:
        assert isinstance(storage, iso_storages.iso_storage)
        s_connections = storage.connections_to_iso_storages

        q = 0  # sum of fluxes
        for c in s_connections:

            if c.right_node is storage:
                q += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

    except Exception as err:
        raise NotImplementedError("An unexpected error occurred.") from err

    return q


def i_fluxes(storage, Isotopologue, **kwargs):
    try:
        assert isinstance(storage, iso_storages.iso_storage)
        s_connections = storage.connections

        qi = 0  # sum of all iso_fluxes
        for c in s_connections:

            flux_i = c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

            if isinstance(c, iso_fluxes.boundary_connection) or c.left_node is storage:

                qi += flux_i

            else:
                qi -= flux_i

    except Exception as err:
        raise NotImplementedError("An unexpected error occurred.") from err

    return qi


class iso_project(object):
    """
    Layer storing isotopes

    TODO: move all connections and fluxes to the layer (Evaporation, Transpiration, Extraction, Seepage) as e.g. Neumann boundary condition

    Description
    ===========
      The
    """

    def __init__(self):
        """
        Constructor of iso_project.

        Projects are the owner of cells and cells are the owner of storages and storages are the owner of fluxconnections.
        """
        self.__cells = []  # list holding all cells associated with this project

        self.__connection_matrix = []

    def get_flux_nodes(self):
        """
        Returns all flux nodes (storages, soil layers, atmospheres,...) associated with this project.
        """
        return iso_storages.flux_node.Dict_of_all_flux_nodes.values()

    def get_iso_storages(self):
        """
        Returns all isotope storages (storages,pond, soil layers, ,...) associated with this project.
        """
        return iso_storages.iso_storage.Dict_of_all_iso_storages.values()

    def get_flux_connections(self):
        """
        Returns all flux connections associated with this project.
        """
        # iso_fluxes.flux_connection.Dict_of_all_flux_connections.values()
        return self.__cells[0].storage_connections + self.__cells[0].boundary_connections

    def get_cells(self):
        """
        Returns all cells associated with this project
        """
        return self.__cells

    def new_cell(self, atmosphere, area, x, y, z):
        """
        Creates a new cell and registers it in the project

        Parameters
        ----------
        atmosphere
        area
        x
        y
        z

        """
        i_point = iso_storages.Point(x, y, z)
        i_cell = iso_cell.iso_cell(location=i_point, atmosphere=atmosphere, area=area)
        self.__cells.append(i_cell)

    def _run(self, Isotopologue, delta_time, delta_cv=[], deltaS_liq=[], **kwargs):

        """
        Runs SLI for one time step. All state variables (Temperatures) and flux parameters (q_l, q_v) need to be updated before running SLI.
        After one run cycle the states of the sytem will be updated

        @param delta_time:
        @type delta_time:
        """
        storages = list(self.get_iso_storages())

        matrix_A = lil_matrix((len(storages), len(storages)))  # create a n x n matrix
        matrix_B = np.zeros(len(storages))  # current storage

        top_boundary = storages[0]
        bot_boundary = storages[-1]
        for i_storage in storages:

            s_index = storages.index(i_storage)

            matrix_A[s_index, s_index] = - i_storage.get_eff_liquid_volume(Isotopologue=Isotopologue,
                                                                           **kwargs) / delta_time \
                                         + storage_fluxes(i_storage, Isotopologue=Isotopologue, **kwargs)

            if i_storage is not top_boundary:
                matrix_A[s_index, s_index - 1] = left_fluxes(storage=i_storage,
                                                             Isotopologue=Isotopologue,
                                                             **kwargs)
            if i_storage is not bot_boundary:
                matrix_A[s_index, s_index + 1] = right_fluxes(storage=i_storage,
                                                              Isotopologue=Isotopologue,
                                                              **kwargs)

            matrix_B[s_index] = i_storage.get_storage_i(Isotopologue=Isotopologue,
                                                        delta_cv=delta_cv[s_index],
                                                        deltaS_liq=deltaS_liq[s_index],
                                                        **kwargs) / delta_time \
                                + i_fluxes(storage=i_storage, Isotopologue=Isotopologue, **kwargs)

        a = matrix_A.tocsr()
        b = matrix_B

        delta_c = spsolve(a, b)

        return delta_c

    def run(self, Isotopologue, delta_time, delta_cv=[], deltaS_liq=[], **kwargs):
        """
        Runs SLI for one time step. All state variables (Temperatures) and flux parameters (q_l, q_v) need to be updated before running SLI.
        After one run cycle the states of the sytem will be updated

        @param delta_time:
        @type delta_time:
        """
        matrix_A, matrix_B = self.coeff_matrix(Isotopologue=Isotopologue,
                                               delta_time=delta_time,
                                               delta_cv=delta_cv,
                                               deltaS_liq=deltaS_liq,
                                               **kwargs)

        a = matrix_A.tocsr()
        b = matrix_B
        delta_c = spsolve(a, b)

        return delta_c

    def coeff_matrix(self, Isotopologue, delta_time, delta_cv=[], deltaS_liq=[], **kwargs):

        A, B = self.flux_matrix(Isotopologue=Isotopologue, **kwargs)
        storages = list(self.get_iso_storages())

        for i_storage in storages:
            s_index = storages.index(i_storage)

            A[s_index, s_index] -= - i_storage.get_eff_liquid_volume(Isotopologue=Isotopologue,
                                                                            **kwargs) / delta_time
            B[s_index] += i_storage.get_storage_i(Isotopologue=Isotopologue,
                                                         delta_cv=delta_cv[s_index],
                                                         deltaS_liq=deltaS_liq[s_index],
                                                         **kwargs) / delta_time
        return A, B

    def flux_matrix(self, Isotopologue, **kwargs):

        storages = list(self.get_iso_storages())
        s_connections = self.get_flux_connections()

        A = lil_matrix((len(storages), len(storages)))  # create a n x n matrix
        B = np.zeros(len(storages))  # current storage

        for c in s_connections:

            if isinstance(c, iso_fluxes.boundary_connection):

                if isinstance(c.left_node, iso_storages.iso_storage):
                    s = c.left_node
                else:
                    s = c.right_node

                index_s = storages.index(s)
                A[index_s, index_s] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                B[index_s] += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)
            else:
                index_l = storages.index(c.left_node)  # index of left storage node
                index_r = storages.index(c.right_node)  # index of right storage node

                if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):
                    A[index_r, index_r] += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                    A[index_l, index_r] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                else:
                    A[index_r, index_r] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                    A[index_l, index_r] += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

                A[index_r, index_l] += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                A[index_l, index_l] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

                B[index_l] += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)
                B[index_r] -= c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

        return A, B

    def mass_balance(self, Isotopologue, dt, **kwargs):
        """"Returns the total mass [kg] of Isotopologue in the system"""

        m = self.storage_mass(Isotopologue=Isotopologue) \
            + self.mass_flux(Isotopologue=Isotopologue, **kwargs) * dt
        return m

    def storage_mass(self, Isotopologue):
        """"Returns the total mass in [kg] of given Isotopologue in the cell"""

        storage_mass = [l.get_conc_iso_liquid(Isotopologue=Isotopologue) * l.thickness * l.theta
                        for l in self.__cells[0].layers]  # kg m-3 * m2 * m * theta
        total_mass = sum(np.array(storage_mass))  # kg

        return total_mass

    def mass_flux(self, Isotopologue, **kwargs):

        """"Returns the net mass flux [kg s-1] leaving the cell"""

        connections = self.__cells[0].boundary_connections
        q = 0
        for c in connections:
            q += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)  # kg s-1

        return q
