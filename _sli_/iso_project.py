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
        return self.__cells[0].boundary_connections + self.__cells[0].storage_connections

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

    def run(self, Isotopologue, delta_time, **kwargs):
        """
        Runs SLI for one time step. All state variables (Temperatures) and flux parameters (q_l, q_v) need to be updated before running SLI.
        After one run cycle the states of the sytem will be updated

        @param delta_time:
        @type delta_time:
        """
        try:
            matrix_A, matrix_B = self.coeff_matrix(Isotopologue=Isotopologue,
                                                   delta_time=delta_time,
                                                   **kwargs)

            a = matrix_A.tocsr()
            b = matrix_B
            delta_c = spsolve(a, b)

        except ValueError as err:
            print(err)
            raise NotImplementedError

        return delta_c

    def coeff_matrix(self, Isotopologue, delta_time, **kwargs):

        try:
            A, B = self.flux_matrix(Isotopologue=Isotopologue, **kwargs)
            storages = list(self.get_iso_storages())

            for i_storage in storages:
                s_index = storages.index(i_storage)

                A[s_index, s_index] -= i_storage.get_eff_liquid_volume(Isotopologue=Isotopologue,
                                                                       **kwargs) / delta_time
                B[s_index] += i_storage.get_storage_i(Isotopologue=Isotopologue,
                                                      **kwargs) / delta_time

        except ValueError as err:
            print(err)
            raise NotImplementedError

        return A, B

    def flux_matrix(self, Isotopologue, **kwargs):

        try:
            storages = list(self.get_iso_storages())
            connections = self.get_flux_connections()

            A = lil_matrix((len(storages), len(storages)))  # create a n x n matrix
            B = np.zeros(len(storages))  # current storage

            for c in connections:

                if isinstance(c, iso_fluxes.boundary_connection):

                    # check which storage node is connected to boundary
                    if isinstance(c.left_node, iso_storages.iso_storage):
                        s = c.left_node
                    elif isinstance(c.right_node, iso_storages.iso_storage):
                        s = c.right_node
                    else:
                        raise ValueError("Boundary connection should have one node connected to iso_storage")

                    index_s = storages.index(s)
                    A[index_s, index_s] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                    B[index_s] += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

                else:
                    index_l = storages.index(c.left_node)  # index of left storage node
                    index_r = storages.index(c.right_node)  # index of right storage node

                    if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):
                        A[index_r, index_r] += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                        A[index_l, index_r] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

                    elif isinstance(c, (iso_fluxes.liquid_diffusion, iso_fluxes.vapor_diffusion)):
                        A[index_r, index_r] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                        A[index_l, index_r] += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                    else:
                        raise ValueError("storage connections should be either advection or diffusion")

                    A[index_r, index_l] += c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                    A[index_l, index_l] -= c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)

                    B[index_l] += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)
                    B[index_r] -= c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

        except ValueError as err:
            print(err)
            raise NotImplementedError

        return A, B

    def mass_balance(self, Isotopologue, dt, delta_cv=[], deltaS_liq=[], **kwargs):
        """"Returns the total mass [kg] of Isotopologue in the system"""

        m = self.storage_mass(Isotopologue=Isotopologue, delta_cv=delta_cv, deltaS_liq=deltaS_liq, **kwargs) \
            - self.mass_flux(Isotopologue=Isotopologue, **kwargs) * dt

        return m

    def storage_mass(self, Isotopologue, delta_cv=[], deltaS_liq=[], **kwargs):
        """"Returns the total mass in [kg] of given Isotopologue in the cell"""

        mass = 0
        for i, l in enumerate(self.__cells[0].layers):
            mass += l.get_conc_iso_liquid(Isotopologue=Isotopologue) * \
                    (l.theta + l.beta(Isotopologue=Isotopologue, **kwargs) * l.cv *
                     (l.theta_sat - l.theta)) * l.thickness
                    #* self.eff_saturation(storage=l, Isotopologue=Isotopologue, **kwargs) * l.thickness

        return mass

    def mass_flux(self, Isotopologue, **kwargs):

        """"Returns the net mass flux [kg s-1] leaving the cell"""

        connections = self.__cells[0].boundary_connections
        q = 0
        for c in connections:

            """"
            if isinstance(c, iso_fluxes.evaporation):
                q += c.q_evaporation_out * c.delta_evapout_c_iso(Isotopologue=Isotopologue, **kwargs) \
                     * c.right_node.get_conc_iso_liquid(Isotopologue=Isotopologue)
            """

            q += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)  # kg s-1

        return q

    def eff_saturation(self, storage, Isotopologue, **kwargs):

        #thre = storage.theta_sat - storage.theta_0
        #eff_sat = storage.S * thre + storage.theta_0 \
         #         + storage.beta(Isotopologue=Isotopologue, **kwargs) * storage.cv * thre * (1 - storage.S)

        #eff_sat =
        #return eff_sat

        pass

