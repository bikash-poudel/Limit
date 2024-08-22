'''
Created on 20.06.2024

@author: poudel-b
'''
# -*- coding: utf-8 -*-

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

from . import iso_cell
from . import iso_fluxes
from . import iso_storages


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

    def run(self, Isotopologue, delta_time, error_tol=1e-16, **kwargs):
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

            # Check for mass balance
            self.mass_balance(delta_c=delta_c, dt=delta_time, Isotopologue=Isotopologue, tolerance=error_tol, **kwargs)

            return delta_c

        except ValueError:
            raise NotImplementedError

    def coeff_matrix(self, Isotopologue, delta_time, **kwargs):

        try:
            storages = list(self.get_iso_storages())

            A = lil_matrix((len(storages), len(storages)))  # create a n x n matrix
            B = np.zeros(len(storages))

            for i_storage in storages:
                s_index = storages.index(i_storage)

                A[s_index, s_index] -= i_storage.get_eff_liquid_volume(Isotopologue=Isotopologue,
                                                                       **kwargs) / delta_time
                B[s_index] += i_storage.get_storage_i(Isotopologue=Isotopologue,
                                                      **kwargs) / delta_time
        except ValueError as err:
            print(err)
            raise NotImplementedError

        return self.flux_matrix(A, B, storages, Isotopologue=Isotopologue, **kwargs)

    def flux_matrix(self, A, B, storages, Isotopologue, **kwargs):

        try:
            for c in self.get_flux_connections():

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

    def mass_balance(self, delta_c, dt, Isotopologue, tolerance=1e-16, **kwargs):

        try:
            LHS = self.storage_mass(delta_c=delta_c, dt=dt, Isotopologue=Isotopologue, **kwargs)
            RHS = self.flux_mass(delta_c=delta_c, Isotopologue=Isotopologue, **kwargs)

            #print('error: ', max(abs(np.array(LHS) - np.array(RHS))))

            if tolerance is not None:
                if max(abs(np.array(LHS) - np.array(RHS))) > tolerance:
                    raise ValueError("Mass balance exceeded the tolerance")
            else:
                return None

        except ValueError:
            raise NotImplementedError

    def storage_mass(self, delta_c, dt, Isotopologue, **kwargs):

        try:
            storages = list(self.get_iso_storages())  # current storage

            if len(storages) != len(delta_c):
                raise ValueError("mismatch length of storages and concentration")
            else:
                mass = []
                for i, s in enumerate(storages):

                    m = (s.get_storage_i(Isotopologue=Isotopologue, **kwargs) +
                         s.get_eff_liquid_volume(Isotopologue=Isotopologue, **kwargs) * delta_c[i]) / dt

                    mass.append(m)

        except ValueError:
            raise NotImplementedError

        return mass

    def flux_mass(self, delta_c, Isotopologue, **kwargs):

        try:
            storages = list(self.get_iso_storages())  # current storage

            mass = np.zeros(len(storages))
            for c in self.get_flux_connections():

                if isinstance(c, iso_fluxes.boundary_connection):

                    # check which storage node is connected to boundary
                    if isinstance(c.left_node, iso_storages.iso_storage):
                        s = c.left_node
                    elif isinstance(c.right_node, iso_storages.iso_storage):
                        s = c.right_node
                    else:
                        raise ValueError("Boundary connection should have one node connected to iso_storage")

                    index_s = storages.index(s)

                    mass[index_s] -= c.calc_flux_i(Isotopologue=Isotopologue, **kwargs) + \
                                     c.calc_flux_liquid(Isotopologue=Isotopologue) * delta_c[index_s]

                else:
                    index_l = storages.index(c.left_node)  # index of left storage node
                    index_r = storages.index(c.right_node)  # index of right storage node

                    if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):

                        mass[index_l] -= c.calc_flux_i(Isotopologue=Isotopologue, **kwargs) + \
                                         c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs) \
                                         * (delta_c[index_l] + delta_c[index_r])

                        mass[index_r] += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs) + \
                                         c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs) * (
                                                 delta_c[index_l] + delta_c[index_r])

                    elif isinstance(c, (iso_fluxes.liquid_diffusion, iso_fluxes.vapor_diffusion)):

                        mass[index_l] -= c.calc_flux_i(Isotopologue=Isotopologue, **kwargs) + \
                                         c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs) \
                                         * (delta_c[index_l] - delta_c[index_r])

                        mass[index_r] += c.calc_flux_i(Isotopologue=Isotopologue, **kwargs) + \
                                         c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs) \
                                         * (delta_c[index_l] - delta_c[index_r])

                    else:
                        raise ValueError("storage connections should be either advection or diffusion")

        except ValueError:
            raise NotImplementedError

        return mass
