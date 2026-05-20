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


class StorageMapping(object):

    def __init__(self):

        super().__init__()

        self._cells = []

        self._storage_dic = {}  # {matrix_index: (cell_obj, layer_obj)}

        self._storage_indexes = {}  # {(cell_id, layer_id): matrix_index}

        self.n_storages = None

        self._mapped = None

    def build_mapping(self):
        """
        Build mapping between matrix indices and storage IDs
        """
        matrix_index = 0

        for cell in sorted(self._cells, key=lambda x: x.get_ID()):

            for layer in sorted(cell.layers, key=lambda x: x.get_ID()):
                # Both forward and reverse mapping
                # creating global index for each storage and storage for each index

                self._storage_dic[matrix_index] = (cell, layer)
                self._storage_indexes[(cell, layer)] = matrix_index

                matrix_index += 1

        self.n_storages = matrix_index

        self._mapped = True

    def get_storage_indexes(self):
        """
        Returns the dictionary of corresponding index associated with each storage (cell, layer) = matrix_index
        """
        return dict(self._storage_indexes)

    def get_storage_index(self, layer):

        """
        Returns the corresponding matrix index for a given storage(cell, layer)
        """
        return self._storage_indexes.get((layer.cell, layer))

    def get_storage_from_index(self, matrix_index):

        """
        Refer to the corresponding storage (cell, layer) from assigned matrix index.
        """

        try:
            return self._storage_dic.get(matrix_index)
        except KeyError:
            raise IndexError(f"Unknown storage index: {matrix_index}")


class cell_connection(object):

    def __init__(self,
                 left_node,
                 right_node,
                 flow_width):

        self._flow_width = flow_width
        try:
            assert isinstance(left_node, iso_cell.cell_node) and isinstance(right_node, iso_cell.cell_node), \
                "Given nodes must be an instance of iso_cell.cell_node"

            self.__left_node = left_node
            self.__right_node = right_node

            # register connection
            left_node.RegisterConnection(newConnection=self)
            right_node.RegisterConnection(newConnection=self)

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    @property
    def flow_width(self):
        return self._flow_width

    def get_left_node(self):
        return self.__left_node

    left_node = property(get_left_node, None, None, "Returns the left node of the flux connection")

    def get_right_node(self):
        return self.__right_node

    right_node = property(get_right_node, None, None, "Returns the right node of the flux connection")


class CellConnector(object):

    def __init__(self):

        super().__init__()
        self.cell_node_connections = []
        self._all_cell_connections = []
        self._dic_cell_connections = {}

    def get_all_cell_connections(self):
        """"Returns every cell connections within the project"""
        return self._all_cell_connections

    def get_cell_connections(self, left_cell, right_cell):

        """Returns the connections between the two  cells"""

        key1 = (left_cell, right_cell)
        key2 = (right_cell, left_cell)

        if key1 in self._dic_cell_connections.keys():
            return self._dic_cell_connections[key1]
        elif key2 in self._dic_cell_connections.keys():
            return self._dic_cell_connections[key2]
        else:
            raise ValueError(f"Cells {left_cell} and {right_cell} are not connected")

    def install_cell_connection(self, flow_width: list,
                                liquid_diffusion=True, vapor_diffusion=True,
                                liquid_advection=True, vapor_advection=True,
                                ):

        """Creats connections between all the cells"""

        cells = self._cells
        for left_cell, right_cell, flow in zip(cells[:-1], cells[1:], flow_width):
            self.connect_cells(left_cell=left_cell,
                               right_cell=right_cell,
                               flow_width=flow,
                               liquid_diffusion=liquid_diffusion, vapor_diffusion=vapor_diffusion,
                               liquid_advection=liquid_advection, vapor_advection=vapor_advection
                               )

    def connect_cells(self, left_cell, right_cell, flow_width,
                      liquid_diffusion=True, vapor_diffusion=True,
                      liquid_advection=True, vapor_advection=True
                      ):

        """Connect layers between two adjacent cells"""

        # create cell_node connections
        conn = cell_connection(left_node=left_cell, right_node=right_cell, flow_width=flow_width)
        self.cell_node_connections.append(conn)

        self._dic_cell_connections[(left_cell, right_cell)] = []
        try:
            # For overlapping layers between two cell_1, cell_2
            for left_node in left_cell.layers:
                for right_node in right_cell.layers:

                    # Check which layers have lateral connection
                    if left_node.connects_to(right_node):
                        # print('connection found')
                        # If connected, define connections and assign to each layers

                        interface = conn.flow_width * self.overlap(l1=left_node, l2=right_node)

                        if liquid_diffusion:
                            ld = iso_fluxes.liquid_diffusion(left_node=left_node,
                                                             right_node=right_node,
                                                             interface_area=interface
                                                             )
                            self._all_cell_connections.append(ld)
                            self._dic_cell_connections[(left_cell, right_cell)].append(ld)

                        if vapor_diffusion:
                            vd = iso_fluxes.vapor_diffusion(left_node=left_node,
                                                            right_node=right_node,
                                                            interface_area=interface)
                            self._all_cell_connections.append(vd)
                            self._dic_cell_connections[(left_cell, right_cell)].append(vd)

                        if liquid_advection:
                            l_ad = iso_fluxes.liquid_advection(left_node=left_node,
                                                               right_node=right_node,
                                                               climate=self._climate)
                            self._all_cell_connections.append(l_ad)
                            self._dic_cell_connections[(left_cell, right_cell)].append(l_ad)

                        if vapor_advection:
                            v_ad = iso_fluxes.vapor_advection(left_node=left_node,
                                                              right_node=right_node,
                                                              climate=self._climate)
                            self._all_cell_connections.append(v_ad)
                            self._dic_cell_connections[(left_cell, right_cell)].append(v_ad)

        except Exception as err:
            raise RuntimeError(f"Error connecting {left_cell} <-> {right_cell}: {err}") from err

    def overlap(self, l1, l2):

        """
        Compute vertical overlap thickness between two soil layers.

        All depths in meters, z=0 at surface, positive downward.
        """

        overlap_top = max(abs(l1.upper_boundary), abs(l2.upper_boundary))
        overlap_bot = min(abs(l1.lower_boundary), abs(l2.lower_boundary))

        overlap = overlap_bot - overlap_top

        if overlap > 0.0:
            return overlap
        else:
            return 0.0


class iso_project(CellConnector, StorageMapping):
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

        self._cells = []  # list holding all cells associated with this project

        self.Isotopologue = None

        self.dc = None

        self._climate = 'humid'  # 0:arid, 1: Humid

        super().__init__()

    @property
    def climate(self):
        return self._climate

    @climate.setter
    def climate(self, value: str):
        try:
            if value not in ['humid', 'arid']:
                raise ValueError("Please enter either 'arid' or 'humid'")

            self._climate = value

        except:
            raise NotImplementedError

    def get_cells(self):
        """
        Returns all cells associated with this project
        """
        return self._cells

    def get_flux_nodes(self):
        """
        Returns all flux nodes (storages, soil layers, atmospheres,...) associated with this project.
        """
        return iso_storages.flux_node.Dict_of_all_flux_nodes.values()

    def get_iso_storages(self):
        """
        Returns all isotope storages (storages,pond, soil layers, ,...) associated with this project.
        """
        # return iso_storages.iso_storage.Dict_of_all_iso_storages.values()
        return [layer for c in self._cells for layer in c.layers]

    def get_flux_connections(self):
        """
        Returns all flux connections associated with this project.

        # return iso_fluxes.flux_connection.Dict_of_all_flux_connections.values()
        connections = self.cell_connections
        for c in self._cells:
            connections.extend(c.boundary_connections + c.storage_connections)
        """

        if not self._cells:
            return []

        connections = []
        for c in self._cells:
            connections.extend(getattr(c, 'boundary_connections', []))
            connections.extend(getattr(c, 'storage_connections', []))

        return connections + self.get_all_cell_connections()

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
        i_cell.project = self
        self._cells.append(i_cell)

    def remove_cell(self, cell):
        """
        Removes the cell from project
        """
        self._cells.remove(cell)

    def update_concentration(self, Isotopologue, delta_c):

        """update iso concentrations to current time step """

        try:
            for c in self._cells:
                # get all the indices for the layers in this cell
                idxs = np.array([self.get_storage_index(l) for l in c.layers])

                # extracts solution values for this cell
                cell_solution_vector = delta_c[idxs]

                c_t = np.array(c.get_conc_layers(Isotopologue=Isotopologue)) + cell_solution_vector
                c.update_c_layers(conc_iso=c_t.tolist(), Isotopologue=Isotopologue)

        except ValueError:
            raise NotImplementedError

    def run(self, Isotopologue, delta_time, error_tolerance=1e-16, check_mass_balance=False, **kwargs):

        """
        @param delta_time:
        @type delta_time:

        @param error_tolerance:
        @type error_tolerance:

        @param safety_factor: safety factor to avoid overshooting
        @type factor:

        @param dt_min: safety factor to avoid overshooting
        @type dt_min:

        @param check_mass_balance:

        """

        self.Isotopologue = Isotopologue

        if not self._mapped and len(self._cells) > 1:
            raise RuntimeError \
                ('Please build storage map for matrix indexing with "indexing()" method before project run')

        dc = self.solve(Isotopologue, delta_time, **kwargs)

        # Check for mass balance
        if check_mass_balance:

            abs_error = self.mass_balance(
                delta_c=dc,
                dt=delta_time,
                Isotopologue=Isotopologue,
                tolerance=error_tolerance,
                **kwargs
            )
            if abs_error > error_tolerance:
                raise f"Mass balance {abs_error} exceeded tolerance {error_tolerance}"

        # print("error:", abs_error)
        self.dc = dc
        self.update_concentration(Isotopologue=Isotopologue, delta_c=self.dc)

        return dc, 0

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
            # storages = list(self.get_iso_storages())
            storages = self.get_iso_storages()

            A = lil_matrix((self.n_storages, self.n_storages), dtype=float)  # create a n x n matrix
            B = np.zeros(self.n_storages, dtype=float)

            for i_storage in storages:
                s_index = self.get_storage_index(i_storage)  # storages.index(i_storage)

                A[s_index, s_index] -= i_storage.get_eff_liquid_volume(Isotopologue=Isotopologue,
                                                                       **kwargs) / delta_time
                B[s_index] += i_storage.get_storage_i(Isotopologue=Isotopologue,
                                                      **kwargs) / delta_time
        except ValueError as err:
            print(err)
            raise NotImplementedError

        return self.flux_matrix(A, B, Isotopologue=Isotopologue, **kwargs)

    def flux_matrix(self, A, B, Isotopologue, **kwargs):

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

                    index_s = self.get_storage_index(s)  # storages.index(s)
                    A[index_s, index_s] -= flux_liquid
                    B[index_s] += flux_iso
                else:
                    index_l = self.get_storage_index(
                        c.left_node)  # storages.index(c.left_node)  # index of left storage node
                    index_r = self.get_storage_index(
                        c.right_node)  # storages.index(c.right_node)  # index of right storage node

                    da, db = self.flux_derivatives(flux_liquid)
                    if isinstance(c, (iso_fluxes.liquid_advection, iso_fluxes.vapor_advection)):

                        A[index_r, index_r] += flux_liquid * db
                        A[index_l, index_r] -= flux_liquid * db

                        A[index_r, index_l] += flux_liquid * da
                        A[index_l, index_l] -= flux_liquid * da

                    elif isinstance(c, (iso_fluxes.liquid_diffusion, iso_fluxes.vapor_diffusion)):
                        A[index_r, index_r] -= flux_liquid
                        A[index_l, index_r] += flux_liquid

                        A[index_r, index_l] += flux_liquid
                        A[index_l, index_l] -= flux_liquid
                    else:
                        raise ValueError("storage connections should be either advection or diffusion")

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

                    index_s = self.get_storage_index(s)  # storage_index[s]
                    mass[index_s] -= flux_i + flux_l * delta_c[index_s]

                else:
                    index_l = self.get_storage_index(
                        c.left_node)  # storages.index(c.left_node)  # index of left storage node
                    index_r = self.get_storage_index(
                        c.right_node)  # storages.index(c.right_node)  # index of right storage node

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

    def flux_derivatives(self, flux):

        if self._climate == 'humid':

            # sli_solve:L2432, dcqldca & dcqldcb
            if flux > 0.0:

                return 1.0, 0.0

            else:
                return 0.0, 1.0

        elif self._climate == 'arid':
            return 0.5, 0.5  # sli_solve:L2432, dcqldca & dcqldcb = 0.5

        else:
            raise ValueError("Wrong climate input, accepted ones are 'arid' or 'humid'")

    def update_lateral_fluxes(slf, ql: dict | None = None):

        """
        Update lateral liquid advection fluxes.

        Parameters
        ----------
        ql : dict
            {(left_cell, right_cell): flux_value}
        """

        if ql is None:
            return

        for nodes, flux in ql.items():

            if nodes not in self._dic_cell_connections:
                raise ValueError(f"Lateral connection between {nodes} does not exist")

            connections = self._dic_cell_connections[nodes]

            updated = False
            for c in connections:
                if isinstance(c, iso_fluxes.liquid_advection):
                    c.q_l = flux
                    updated = True

            if not updated:
                raise ValueError(
                    f"Missing liquid advection connection between {nodes}"
                )

    def delta_S(self, storage):

        l = storage

        S = l.S
        S0 = l.calc_S(l.theta_sat, l.theta_0, l.theta_t0)

        return S - S0

    def delta_mass(self):

        m = []
        for s in self.get_iso_storages():
            idx = self.get_storage_index(s)

            m.append(self.dc[idx] * s.cell.area * s.thickness)

        return m

    def delta_flux_mass(self, delta_t):

        m = []

        for s in self.get_iso_storages():
            m.append(s.flux_mass(Isotopologue=self.Isotopologue,
                                 delta_t=delta_t
                                 )
                     )
        return m
