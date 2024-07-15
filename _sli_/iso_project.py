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
        return iso_fluxes.flux_connection.Dict_of_all_flux_connections.values()

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

    def run(self, Isotopologue, delta_time, delta_cv=[], deltaS_liq=[], **kwargs):

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

            storage_index = storages.index(i_storage)

            storage_volume = i_storage.get_eff_liquid_volume(Isotopologue=Isotopologue, **kwargs)
            i_storage_volume = i_storage.get_storage_i(Isotopologue=Isotopologue,
                                                       delta_cv=delta_cv[storage_index],
                                                       deltaS_liq=deltaS_liq[storage_index],
                                                       **kwargs)

            left_c_diff_fluxes_to_left, left_c_diff_fluxes_to_right, \
                left_c_adv_fluxes_to_left, left_c_adv_fluxes_to_right = 0, 0, 0, 0
            right_c_diff_fluxes_to_left, right_c_diff_fluxes_to_right, \
                right_c_adv_fluxes_to_left, right_c_adv_fluxes_to_right = 0, 0, 0, 0
            left_fluxes_i, right_fluxes_i = 0, 0
            boundary_fluxes, boundary_fluxes_i = 0, 0

            # Boundary Connections to each storage
            for boundary_connection in i_storage.connections_to_boundaries:
                if not isinstance(boundary_connection, iso_fluxes.precipitation):
                    # sli.solve.f90: L2521
                    boundary_fluxes += boundary_connection.calc_flux(Isotopologue=Isotopologue, **kwargs)
                boundary_fluxes_i += boundary_connection.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

            # Two-way fluxes in left connections
            for l_connection in i_storage.connections_to_left:
                left_fluxes_i += l_connection.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

                try:
                    # For diffusive fluxes in left connection
                    if isinstance(l_connection, iso_fluxes.vapor_diffusion) or \
                            isinstance(l_connection, iso_fluxes.liquid_diffusion):

                        left_c_diff_fluxes_to_left += l_connection.flux_to_left(Isotopologue=Isotopologue, **kwargs)
                        left_c_diff_fluxes_to_right += l_connection.flux_to_right(Isotopologue=Isotopologue, **kwargs)

                    # For advective fluxes in left connection
                    else:
                        left_c_adv_fluxes_to_left += l_connection.flux_to_left(Isotopologue=Isotopologue, **kwargs)
                        left_c_adv_fluxes_to_right += l_connection.flux_to_right(Isotopologue=Isotopologue, **kwargs)

                except Exception as err:
                    raise NotImplementedError("An unexpected error occurred.") from err

            # Two-way fluxes in right connections
            for r_connection in i_storage.connections_to_right:
                right_fluxes_i += r_connection.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

                try:
                    # For diffusive fluxes in right connection
                    if isinstance(r_connection, iso_fluxes.vapor_diffusion) or \
                            isinstance(r_connection, iso_fluxes.liquid_diffusion):

                        right_c_diff_fluxes_to_left += r_connection.flux_to_left(Isotopologue=Isotopologue, **kwargs)
                        right_c_diff_fluxes_to_right += r_connection.flux_to_right(Isotopologue=Isotopologue, **kwargs)
                    # For advective fluxes in right connection
                    else:
                        right_c_adv_fluxes_to_left += r_connection.flux_to_left(Isotopologue=Isotopologue, **kwargs)
                        right_c_adv_fluxes_to_right += r_connection.flux_to_right(Isotopologue=Isotopologue, **kwargs)

                except Exception as err:
                    raise NotImplementedError("An unexpected error occurred.") from err

            matrix_A[storage_index, storage_index] = - storage_volume / delta_time - boundary_fluxes \
                                                     + left_c_adv_fluxes_to_right - right_c_adv_fluxes_to_left \
                                                     - left_c_diff_fluxes_to_right - right_c_diff_fluxes_to_left

            if i_storage is not top_boundary:
                matrix_A[storage_index, storage_index - 1] = left_c_adv_fluxes_to_left + left_c_diff_fluxes_to_left

            if i_storage is not bot_boundary:
                matrix_A[storage_index, storage_index + 1] = - right_c_adv_fluxes_to_right \
                                                             + right_c_diff_fluxes_to_right

            matrix_B[storage_index] = i_storage_volume / delta_time + boundary_fluxes_i - left_fluxes_i + right_fluxes_i

        a = matrix_A.tocsr()
        b = matrix_B

        delta_c = spsolve(a, b)

        return delta_c

    """
    def mass_balance(self, Isotopologue, delta_time, delta_c, **kwargs):

        storages = list(self.get_iso_storages())

        i_storage = storages[0]
        storage_index = storages.index(storage)
        delta_storage = storage.thetasat * storage.thickness / delta_time * \
                        (delta_c[storage_index] * storage.eff_saturation(Isotopologue, **kwargs))

        # Boundary Connections to each storage
        for boundary_connection in storage.connections_to_boundaries:
            if not isinstance(boundary_connection, iso_fluxes.precipitation):
                # sli.solve.f90: L2521
                boundary_fluxes += boundary_connection.calc_flux(Isotopologue=Isotopologue, **kwargs)
            boundary_fluxes_i += boundary_connection.calc_flux_i(Isotopologue=Isotopologue, **kwargs)
        #delta_fluxes =
        
    """
