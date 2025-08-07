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


def estimate_dt(current_dt, tolerance_ratio, sf=0.9, alpha=0.5):
    """
    Computes new step size (adaptive time-stepping)

    @param sf: safety factor
    @type sf:

    @param alpha: safety factor
    @type alpha:

    """

    return sf * current_dt * tolerance_ratio ** alpha


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

    def update_concentration(self, Isotopologue, delta_c):

        """update iso concentrations to current time step """

        c = self.get_cells()[0]

        c_t = list(np.array(c.get_conc_layers(Isotopologue=Isotopologue)) + np.array(delta_c))
        c.update_c_layers(conc_iso=c_t, Isotopologue=Isotopologue)

    def run(self, Isotopologue, delta_time, error_tolerance=1e-16, safety_factor=0.9, dt_min=1e-8, **kwargs):

        """
        @param delta_time:
        @type delta_time:

        @param error_tolerance:
        @type error_tolerance:

        @param safety_factor: safety factor to avoid overshooting
        @type factor:

        @param dt_min: safety factor to avoid overshooting
        @type dt_min:

        """

        # solve for the current time step
        dc = self.solve(Isotopologue, delta_time, **kwargs)

        # Check for mass balance
        abs_error = self.mass_balance(delta_c=dc,
                                      dt=delta_time,
                                      Isotopologue=Isotopologue,
                                      tolerance=error_tolerance,
                                      **kwargs)

        print('abs_error: ', abs_error)

        # update the storage to current time step when error is within desirable limits
        self.update_concentration(Isotopologue=Isotopologue, delta_c=dc)

        """
        alpha = 0.5  # exponent for implicit Euler (Backward differentiation)
        cum_dt = 0
        current_dt = delta_time

        # estimate rate of change of storage (states) for given time steps

        delta_storage = self.rate_storage_change(dt=delta_time)
        while cum_dt < delta_time - dt_min:  # minor tolerance to avoid floating point issues

            try:
                # solve for the current time step
                dc = self.solve(Isotopologue, current_dt, **kwargs)

                # Check for mass balance
                abs_error = self.mass_balance(delta_c=dc,
                                              dt=current_dt,
                                              Isotopologue=Isotopologue,
                                              tolerance=error_tolerance,
                                              **kwargs)

                print('abs_error: ', abs_error)

                ################# check error tolerance and adjust time steps #####################
                tolerance_ratio = error_tolerance / abs_error
                if tolerance_ratio >= 1:

                    # update the storage to current time step when error is within desirable limits
                    self.update_concentration(Isotopologue=Isotopologue, delta_c=dc)

                    ######################## estimate next time step ###############################
                    cum_dt += current_dt
                    dt_new = estimate_dt(current_dt, tolerance_ratio, safety_factor, alpha)

                    # ensure new time step within the remaining time limits
                    current_dt = max(dt_min, min(dt_new, delta_time - cum_dt))  # Prevents overshoot

                    # update the interpolated storage corresponding the new time step
                    self.update_storage(delta_storages=delta_storage, dt=current_dt)

                else:

                    ######################## Reject and reduce step ################################
                    print(f"[REJECT] dt={current_dt:.2e} error={abs_error:.2e}")
                    current_dt = estimate_dt(current_dt, tolerance_ratio, safety_factor, alpha)

                    if current_dt < dt_min:
                        raise RuntimeError(f"Time step dropped below min_dt={dt_min}, stopping.")

                    ####################### update states to current time steps ######################
                    # linear interpolation of states over time (checks whether these storage changes are accepetd
                    # and updated later if they are accepted.)
                    self.check_update_storages(delta_storages=delta_storage, dt=current_dt)

            except ValueError:
                raise NotImplementedError("Matrix solver or state updates failed.")
         
         """

        return dc, abs_error

    def solve(self, Isotopologue, delta_time, **kwargs):

        """
        Runs  for one time step. All state variables and flux parameters (q_l, q_v) need to be updated before running.
        After one run cycle the states of the system needs to be updated.

        @param delta_time:
        @type delta_time:
        """

        try:

            matrix_A, matrix_B = self.coeff_matrix(Isotopologue=Isotopologue,
                                                   delta_time=delta_time,
                                                   **kwargs)
            a = matrix_A.tocsr()
            b = matrix_B

            return spsolve(a, b)

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

                flux_liquid = c.calc_flux_liquid(Isotopologue=Isotopologue, **kwargs)
                flux_iso = c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

                if isinstance(c, iso_fluxes.boundary_connection):

                    # check which storage node is connected to boundary
                    if isinstance(c.left_node, iso_storages.iso_storage):
                        s = c.left_node
                    elif isinstance(c.right_node, iso_storages.iso_storage):
                        s = c.right_node
                    else:
                        raise ValueError("Boundary connection should have one node connected to iso_storage")

                    index_s = storages.index(s)
                    A[index_s, index_s] -= flux_liquid
                    B[index_s] += flux_iso

                else:
                    index_l = storages.index(c.left_node)  # index of left storage node
                    index_r = storages.index(c.right_node)  # index of right storage node
                    if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):
                        A[index_r, index_r] += flux_liquid
                        A[index_l, index_r] -= flux_liquid

                    elif isinstance(c, (iso_fluxes.liquid_diffusion, iso_fluxes.vapor_diffusion)):
                        A[index_r, index_r] -= flux_liquid
                        A[index_l, index_r] += flux_liquid
                    else:
                        raise ValueError("storage connections should be either advection or diffusion")

                    A[index_r, index_l] += flux_liquid
                    A[index_l, index_l] -= flux_liquid

                    B[index_l] += flux_iso
                    B[index_r] -= flux_iso

        except ValueError as err:
            print(err)
            raise NotImplementedError

        return A, B

    def mass_balance(self, delta_c, dt, Isotopologue, **kwargs):

        try:

            LHS = self.storage_mass(delta_c=delta_c, dt=dt, Isotopologue=Isotopologue, **kwargs)
            RHS = self.flux_mass(delta_c=delta_c, Isotopologue=Isotopologue, **kwargs)

            # return maximum error
            return max(abs(LHS) - abs(RHS))

        except ValueError:
            raise NotImplementedError

    def storage_mass(self, delta_c, dt, Isotopologue, **kwargs):

        try:
            storages = list(self.get_iso_storages())  # current storage

            if len(storages) != len(delta_c):
                raise ValueError("Mismatch length of storages and concentration")

            liquid_storage = np.array([
                s.get_eff_liquid_volume(Isotopologue=Isotopologue, **kwargs)
                for s in storages])

            iso_storage = np.array([
                s.get_storage_i(Isotopologue=Isotopologue, **kwargs)
                for s in storages])

            delta_c = np.asarray(delta_c)

            storage_mass = (iso_storage + liquid_storage * delta_c) / dt

        except ValueError:
            raise NotImplementedError

        return storage_mass

    def flux_mass(self, delta_c, Isotopologue, **kwargs):

        try:
            storages = list(self.get_iso_storages())  # current storage
            delta_c = np.asarray(delta_c)

            storage_index = {s: i for i, s in enumerate(storages)}
            mass = np.zeros(len(storages))

            for c in self.get_flux_connections():

                flux_l = c.calc_flux_liquid(Isotopologue=Isotopologue)
                flux_i = c.calc_flux_i(Isotopologue=Isotopologue, **kwargs)

                if isinstance(c, iso_fluxes.boundary_connection):

                    # check which storage node is connected to boundary
                    if isinstance(c.left_node, iso_storages.iso_storage):
                        s = c.left_node
                    elif isinstance(c.right_node, iso_storages.iso_storage):
                        s = c.right_node
                    else:
                        raise ValueError("Boundary connection should have one node connected to iso_storage")

                    index_s = storage_index[s]
                    mass[index_s] -= flux_i + flux_l * delta_c[index_s]

                else:
                    index_l = storages.index(c.left_node)  # index of left storage node
                    index_r = storages.index(c.right_node)  # index of right storage node

                    if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):

                        flux = flux_i + flux_l * (delta_c[index_l] + delta_c[index_r])

                        mass[index_l] -= flux

                        mass[index_r] += flux

                    elif isinstance(c, (iso_fluxes.liquid_diffusion, iso_fluxes.vapor_diffusion)):

                        flux = flux_i + flux_l * (delta_c[index_l] - delta_c[index_r])

                        mass[index_l] -= flux

                        mass[index_r] += flux

                    else:
                        raise ValueError("storage connections should be either advection or diffusion")

        except ValueError:
            raise NotImplementedError

        return mass

    def rate_storage_change(self, dt):

        """
        rate of change of state variables

        @param dt: initial time step
        """

        c = self.get_cells()[0]

        delta_storages = {"delta_theta": np.array([l.theta - l.theta_t0 for l in c.layers]) / dt,
                          "delta_h": np.array([l.psi - l.psi_0 for l in c.layers]) / dt,
                          "delta_T": np.array([l.dT for l in c.layers]) / dt}

        return delta_storages

    def check_update_storages(self, delta_storages, dt):

        """
        update states to current time step

        @param delta_storages: dictionary of time derivative of storages

        @param dt: current time step
        """

        c = self.get_cells()[0]

        # update interpolated states to current time steps
        d_theta = delta_storages["delta_theta"]
        d_psi = delta_storages["delta_h"]
        d_T = delta_storages["delta_T"]

        # update storage with interpolated states for current time
        # without changing the states for previous time
        for l, dth, dh, dT in zip(c.layers, d_theta, d_psi, d_T):

            l.theta = l.theta + dth * dt
            l.psi = l.psi + dh * dt
            l.T = l.T + dT * dt

    def update_storage(self, delta_storages, dt):

        """
        update states to current time step

        @param delta_storages: dictionary of time derivative of storages

        @param dt: current time step
        """

        c = self.get_cells()[0]

        # current states
        _theta = np.array([l.theta for l in c.layers])
        _psi = np.array([l.psi for l in c.layers])
        _T = np.array([l.T for l in c.layers])

        # update storage with interpolated states for current time
        # and replace the states for previous time
        c.update_layers(theta=delta_storages["delta_theta"] * dt + _theta,
                        T=delta_storages["delta_T"] * dt + _T,
                        rH=None,
                        psi=delta_storages["delta_h"] * dt + _psi)

