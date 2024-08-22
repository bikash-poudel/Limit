'''
Created on 23.05.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

from . import iso_storages
from . import iso_fluxes
from . import vapor_flux
from . import iso_delta


class iso_cell(object):
    """
    Description
    ===========
    This class is the basic landscape object. It is the owner of water storages, the upper and lower boundary conditions
    of the system (atmosphere, deep groundwater) and the vegetation.

    Calculates the movement of isotopes in the gaseous and liquide phase of cells
    and stores the isotope concentration in each layer.

    The isotope species transport is calculated once the heat and water transport model have been solved for (t+1)
    to provide the soil temperature and water content at a future state.

    Based on SLI by Haverd and Cuntz in cable_sli_solve.f90;  SUBROUTINE isotope_vap ;lines 2149
    """

    Dict_of_all_cells = {}
    current_ID = 1

    def __init__(self,
                 location,
                 atmosphere,  # instance of iso_atmosphere
                 area=1.0,
                 # area in m2 --> fluxes are all in m3/d and only by multiplying them with delta time you get a volume
                 ):
        """
        @param location: Instance of iso_storages.Point
        @type location: tuple of size 3 (x, y, z)
        """

        self.__ID = self.__class__.current_ID
        self.__class__.current_ID += 1  # increment the current ID for the next flux node to be created
        # self.__class__.Dict_of_all_cells[self.__ID] = self  # add Instance of flux node to Dict_of_all_flux_nodes

        assert (isinstance(location, iso_storages.Point)), "The surface center of a cell must be an instance of " \
                                                           "iso_Storages.Point"
        assert (isinstance(atmosphere, iso_storages.iso_atmosphere)), "Atmosphere must be an instance of " \
                                                                      "iso_storages.iso_atmosphere"
        # Cell Properties
        self.location = location
        self.area = area

        # Upper Boundaries
        self.__atmosphere = atmosphere  # iso>storage.iso_atmosphere

        self.__Ts = 283.15  # Surface temperature [K]
        self.__q_evap = 0.0  # evaporation flux
        self.__ql_surface = 0.0  # liquid fluxes into the surface
        self.__qv_surface = 0.0  # vapor fluxes into the surface
        self.__q_runoff = 0.0  # surface runoff
        self.__q_precipitation = 0.0  # precipitation flux
        self.__c_precipitation = {"2H": 1.0, "18O": 1.0}  # iso conc in precipitation
        self.__q_neuman = 0.0  # neuman flux
        self.__c_neuman = {"2H": 1.0, "18O": 1.0}  # iso conc in neuman flux
        self.__c_dirichlet = {"2H": 1.0, "18O": 1.0}  # iso conc as dirichlet boundary

        self.__pond = None

        # Soil Layers
        self.__layers = []
        self.__liquid_fluxes = []  # liquid fluxes within the layers
        self.__vapor_fluxes = []  # vapor fluxes within the soil
        self.__transpiration = []  # transpiration fluxes in the soil layers
        self.__conc_layers = []  # list of isotopic concentration in soil layers

        # Lower boundaries
        self.__aquifer = None
        # Storage Connections
        self.__connection_l_diff = []
        self.__connection_v_diff = []
        self.__connection_l_adv = []
        self.__connection_v_adv = []
        # Boundary Connections
        self.__connection_ev = None
        self.__connections_trans = None
        self.__connection_prec = None
        self.__connection_runoff = None
        self.__connection_to_aquifer = None
        self.__connection_neuman = None
        self.__connection_dirichlet = None

    """Properties"""

    @property
    def Ts(self):
        return self.__Ts

    @property
    def flux_nodes(self):
        """
        @return: Returns a list of flux nodes belonging to this cell (soil_layers and atmosphere).
        """
        return [self.atmosphere] + self.layers + [self.aquifer]

    @property
    def atmosphere(self):
        """
        @return: Returns the atmosphere of the cell.
        """
        return self.__atmosphere

    @property
    def layers(self):
        """
        @return: Returns a list of layers belonging to this cell.
        """
        return self.__layers

    @property
    def conc_2H(self):
        """"Returns the concentration of 2H in the soil layers"""

        c_2H = []
        for l in self.__layers:
            c_2H.append(l.get_conc_iso_liquid('2H'))

        return c_2H

    @property
    def conc_2H_delta(self):

        delta_2H = []
        for l in self.__layers:
            c = l.get_conc_iso_liquid("2H")
            delta_2H.append(iso_delta.concentration_to_delta(c, "2H"))

        return delta_2H

    @property
    def conc_18O(self):
        """"Returns the concentration of 2H in the soil layers"""

        c_18O = []
        for l in self.__layers:
            c_18O.append(l.get_conc_iso_liquid('18O'))

        return c_18O

    @property
    def conc_18O_delta(self):

        delta_18O = []
        for l in self.__layers:
            c = l.get_conc_iso_liquid("18O")
            delta_18O.append(iso_delta.concentration_to_delta(c, "18O"))

        return delta_18O

    @property
    def iso_storages(self):
        """
        Returns all the flux nodes that are the instances of iso_storages
        -------
        """
        try:
            storages = []
            for node in self.flux_nodes:
                if isinstance(node, iso_storages.iso_storage):
                    storages.append(node)

            return storages

        except Exception as err:
            raise NotImplementedError("An unexpected error occurred.") from err

    @property
    def top_layer(self):
        """"
        Returm the top boundary layer in the cell
        """
        return self.__layers[0]

    @property
    def bottom_layer(self):
        """"
        Returm the bottom boundary layer in the cell
        """
        return self.__layers[-1]

    @property
    def pond(self):
        """
        @return: Returns the aquifer isotope storage of the cell.
        """
        return self.__pond

    @property
    def aquifer(self):
        """
        @return: Returns the aquifer isotope storage of the cell.
        """
        return self.__aquifer

    @property
    def liquid_fluxes(self):
        """
       @return: Returns the liquid fluxes within the layers cell.
       """
        return self.__liquid_fluxes

    @property
    def vapor_fluxes(self):
        """
       @return: Returns the vapor fluxes within the layers cell.
       """
        return self.__vapor_fluxes

    @property
    def q_evap(self):

        """
        @return: Returns the evaporation flux into the atmosphere.
        """
        return self.__q_evap

    @property
    def ql_surface(self):
        """
        @return: Returns the liquid flux into the cell.
        """
        return self.__ql_surface

    @property
    def qv_surface(self):
        """
        @return: Returns the vapor flux into the cell.
        """
        return self.__qv_surface

    @property
    def q_prec(self):
        """
        @return: Returns the precipitation flux into the cell.
        """
        return self.__q_precipitation

    @property
    def c_prec(self):
        """
        @return: Returns the isotope concentration in precipitation.
        """
        return self.__c_precipitation

    @property
    def q_runoff(self):
        """
        @return: Returns the surface runoff flux.

        """
        return self.__q_runoff

    @property
    def q_transpiration(self):
        """
        @return: Returns the transpiration fluxes from the soil layers.

        """
        return self.__transpiration

    @property
    def q_neuman(self):
        return self.__q_neuman

    @property
    def c_neuman(self):
        return self.__c_neuman

    @property
    def c_dirichlet(self):
        return self.__c_dirichlet

    @property
    def connections_l_diff(self):
        """
        @return: Returns the flux connection related to liquid diffusion.

        """
        return self.__connection_l_diff

    @property
    def connections_v_diff(self):
        """
        @return: Returns the flux connection related to vapor diffusion.

        """
        return self.__connection_v_diff

    @property
    def connections_l_adv(self):
        """
        @return: Returns the flux connection related to liquid advection.

        """
        return self.__connection_l_adv

    @property
    def connections_v_adv(self):
        """
        @return: Returns the flux connection related to vapor advection.

        """
        return self.__connection_v_adv

    @property
    def connection_evap(self):
        """
        @return: Returns the flux connection related to vapor advection.

        """
        return self.__connection_ev

    @property
    def connections_transpiration(self):
        """
        @return: Returns the flux connection related to transpiration.

        """
        return self.__connections_trans

    @property
    def connection_prec(self):
        """
        @return: Returns the flux connection related to transpiration.

        """
        return self.__connection_prec

    @property
    def connection_runoff(self):
        """
        @return: Returns the flux connection related to transpiration.

        """
        return self.__connection_runoff

    @property
    def connection_to_aquifer(self):
        """
        @return: Returns the flux connection related to transpiration.

        """
        return self.__connection_to_aquifer

    @property
    def connection_neuman(self):
        """
        @return: Returns the neuman boundary flux connection related.

        """
        return self.__connection_neuman

    @property
    def connection_dirichlet(self):
        """
        @return: Returns the dirichlet boundary connection related.

        """
        return self.__connection_dirichlet

    @property
    def storage_connections(self):
        return self.connections_l_adv + self.connections_v_adv + self.connections_l_diff + self.connections_v_diff

    @property
    def boundary_connections(self):

        connections = [self.connection_neuman,
                       self.connection_evap,
                       self.connection_prec,
                       self.connection_runoff,
                       self.connection_to_aquifer]

        if self.connections_transpiration is not None:
            connections.extend(self.connections_transpiration)

        # Filter out any None values from the connections list
        connections = [conn for conn in connections if conn is not None]

        return connections

    def is_upper_boundary(self, storage):

        """"
        Checks whether the storage is the upper boundary layer
        """

        try:
            return storage == self.top_layer

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def is_lower_boundary(self, storage):

        """"
        Checks whether the storage is the lower boundary layer
        """

        try:
            return storage == self.bottom_layer

        except ValueError as err:
            print(err)
            raise NotImplementedError

    """Functions"""
    def get_conc_layers(self, Isotopologue):

        """"List of Isotope concentration in layers """

        c = []
        for l in self.__layers:
            c.append(l.get_conc_iso_liquid(Isotopologue))

        return c

    def add_layer(self, new_layer):
        try:
            assert (isinstance(new_layer, iso_storages.iso_soil_layer)), \
                "new_layer must be an instance of iso_storages.iso_soil_layer"

            new_layer.cell = self
            self.__layers.append(new_layer)

            self.__layers.sort(key=lambda
                x: x.lower_boundary)  # sort them to make sure that they are in the right order from top to lowest

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def install_connections(self, liquid_diffusion=True, vapor_diffusion=True,
                            liquid_advection=True, vapor_advection=True):
        """"
        Registers all the relevant connections within the layers
        """

        try:
            for left_node, right_node in zip(self.layers[:-1], self.layers[1:]):

                # Define connections and assign to each layers
                if liquid_diffusion:
                    lD = iso_fluxes.liquid_diffusion(left_node=left_node, right_node=right_node)
                    self.__connection_l_diff.append(lD)
                if vapor_diffusion:
                    vD = iso_fluxes.vapor_diffusion(left_node=left_node, right_node=right_node)
                    self.__connection_v_diff.append(vD)
                if liquid_advection:
                    ladv = iso_fluxes.liquid_advection(left_node=left_node, right_node=right_node)
                    self.__connection_l_adv.append(ladv)
                if vapor_advection:
                    vadv = iso_fluxes.vapor_advection(left_node=left_node, right_node=right_node)
                    self.__connection_v_adv.append(vadv)

        except Exception as err:
            raise NotImplementedError("An unexpected error occurred.")

    def add_pond(self, pond):
        """
        Assign pond at top boundary ..

        Returns
        -------
        """
        try:
            assert isinstance(pond, iso_storages.iso_pond), 'Pond should be an instance of iso_storages.iso_pond'

            self.__layers[0].add_pond(pond)
            self.__pond = self.__layers[0].pond

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def add_evaporation(self, hydrodynamic_dispersivity=0.0):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            # Assign Atmospheric boundary

            ev = iso_fluxes.evaporation(atmosphere=self.atmosphere, top_layer=self.top_layer, q_Evap=self.q_evap,
                                        T_surface=self.Ts, q_l=self.ql_surface, q_v=self.qv_surface,
                                        hydrodynamic_dispersivity=hydrodynamic_dispersivity)
            self.__connection_ev = ev

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    def add_transpiration(self):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            self.__transpiration = [0.0] * len(self.__layers)
            tr = []
            for layer, ql in zip(self.layers, self.q_transpiration):
                # Assign Atmospheric boundary to each layer
                tp = iso_fluxes.transpiration(atmosphere=self.atmosphere, soil_layer=layer, ql_transpiration=ql)

                tr.append(tp)

            self.__connections_trans = tr

        except ValueError:
            raise NotImplementedError("number of transpiration fluxes must be equal to that of soil layers.")

    def add_surface_runoff(self):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            # Assign Atmospheric boundary
            rn = iso_fluxes.surface_runoff(atmosphere=self.atmosphere, top_layer=self.top_layer, q_runoff=self.q_runoff)

            self.__connection_runoff = rn

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    def add_precipitation(self):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            # Assign Atmospheric boundary
            pr = iso_fluxes.precipitation(atmosphere=self.atmosphere, top_layer=self.top_layer,
                                          q_prec=self.q_prec, c_prec=self.c_prec)

            self.__connection_prec = pr

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def add_aquifer(self, aquifer, soil_layer, ql_layer=None):
        """
       Assign boundary connection to water storage body..

       Returns
       -------
       
        # soil_layer connected to the aquifer
        """

        self.__aquifer = aquifer
        try:
            # Assign Atmospheric boundary
            if soil_layer is self.bottom_layer:
                if self.__liquid_fluxes:
                    ql_aq = self.__liquid_fluxes[-1]  # if not None
                else:
                    ql_aq = None
            else:
                ql_aq = ql_layer

            aq = iso_fluxes.aquifer_connection(soil_layer=soil_layer, aquifer=self.aquifer, ql_layer=ql_aq)
            self.__connection_to_aquifer = aq

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def add_neuman_boundary(self, soil_layer):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            # Assign Atmospheric boundary
            nm = iso_fluxes.neuman_boundary(atmosphere=self.atmosphere, soil_layer=soil_layer,
                                            q_neuman=self.q_neuman, c_neuman=self.c_neuman)

            self.__connection_neuman = nm

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def add_dirichlet_boundary(self, soil_layer):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            # Assign dirichlet  boundary

            dr = iso_fluxes.dirichlet_boundary(atmosphere=self.atmosphere, soil_layer=soil_layer,
                                               c_layer=self.c_dirichlet)

            self.__connection_dirichlet = dr

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    """"Update storages"""

    def update_atmosphere(self, c_atm={'2H': 1.0, '18O': 1.0}, T=273.0, Rh=0.2, Pa=10 ** 5, wind_speed=2,
                          hc=10, d0=0.67 * 10, z0m=0.1 * 10, LAI=1.1, extku=1.5):

        """"
        updates the atmospheric states to current time step.
        """
        try:
            atm = self.__atmosphere
            atm.T = T  # Temperature [k]
            atm.Rh = Rh  # atmospheric Relative Humidity
            atm.Pa = Pa  # Atmospheric pressure
            atm.wind_speed = wind_speed  # wind speed

            # Update isotopic concentration
            atm.set_conc_iso_vapor(c_atm['2H'], '2H')
            atm.set_conc_iso_vapor(c_atm['18O'], '18O')

            self.__atmosphere = atm
        except ValueError:
            raise NotImplementedError("A required value was not provided.")

    def update_layers(self, theta, T, rH, psi):

        """
        updates the  states to current time step.
        """
        try:
            if rH is None:
                rH = [None] * len(self.__layers)

            if psi is None:
                raise ValueError("Missing matric potential (psi) for all soil layers")

            for lr, th, T_soil, r_H, pot in zip(self.__layers, theta, T, rH, psi):
                lr.theta_t0 = lr.theta  # u[date theta for previous time
                lr.theta = th  # update theta for current time

                lr.T0 = lr.T  # update temperature for previous time
                lr.T = T_soil  # updata current temperature Temperature

                lr.rH = r_H  # update soil relative humidity
                lr.psi_0 = lr.psi  # update matric potential previosu time
                lr.psi = pot  # updat matric potential current time

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def update_c_layers(self, conc_iso, Isotopologue):
        """
        updates the isotope concentration to current time step.
        """

        try:
            if len(conc_iso) != len(self.__layers):
                raise ValueError("number concentration iso must be equal to the number of layers")
            else:
                for lr, c in zip(self.__layers, conc_iso):
                    lr.set_conc_iso_liquid(c, Isotopologue)

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    def update_pond(self, pond_height):
        """
        updates the isotope concentration to current time step.
        """

        self.__layers[0].pond.pond_height_t0 = self.__pond.pond_height
        self.__layers[0].pond.pond_height = pond_height

        self.__pond = self.__layers[0].pond

    def update_aquifer(self, c_iso={"2H": 1.0, "18O": 1.0}):
        """
        Updates the iso concentration in aquifer to current time step

        Returns
        -------
        """
        try:
            self.__aquifer.set_conc_iso_liquid(c_iso["2H"], "2H")
            self.__aquifer.set_conc_iso_liquid(c_iso["18O"], "18O")

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    """"Update Fluxes"""

    def update_liquid_fluxes(self, liquid_fluxes):

        """"
        updates the list of liquid fluxes within the soil layers for current time step
        """

        try:
            if len(liquid_fluxes) == len(self.__layers):
                self.__liquid_fluxes = liquid_fluxes
            else:
                raise ValueError("Number of liquid fluxes must be equal to the cell layers")

            # assign liquid fluxes to the connections
            if len(liquid_fluxes[:-1]) == len(self.connections_l_adv):

                for ql, c_l_adv in zip(self.liquid_fluxes, self.connections_l_adv):
                    c_l_adv.q_l = ql
            else:
                raise ValueError("Number of liquid fluxes must be equal to the flux connections")

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def update_vapor_fluxes(self, vapor_fluxes=None):

        """"
        updates the list of vapor fluxes within the soil layers for current time step
        """

        try:
            if vapor_fluxes is not None:

                if len(vapor_fluxes) == len(self.__layers):
                    self.__vapor_fluxes = vapor_fluxes
                else:
                    raise ValueError("Number of vapor fluxes must be equal to the cell layers")

                # assign vapor fluxes to the connections
                if len(vapor_fluxes[:-1]) == len(self.connections_v_adv):

                    for qv, c_v_adv in zip(self.vapor_fluxes, self.connections_v_adv):
                        c_v_adv.q_v = qv
                else:
                    raise ValueError("Number of vapor fluxes must be equal to the flux connections")

            else:
                q = vapor_flux.vapor_flux()
                q_v = []
                for i, c_v_adv in enumerate(self.connections_v_adv):
                    qv = q.q_vapor(left_node=c_v_adv.left_node, right_node=c_v_adv.right_node)
                    c_v_adv.q_v = qv

                    q_v.append(qv)

                self.__vapor_fluxes = q_v

        except ValueError:
            raise NotImplementedError

    def update_evaporation(self, q_ev=0.0, T_surface=273.0, ql_surface=0.0, qv_surface=0.0):
        """
        Updates the flux related to iso - evaporation to current time step

        Returns the evaporative connection to the cell
        -------
        """

        self.__Ts = T_surface
        self.__q_evap = q_ev
        self.__ql_surface = ql_surface
        self.__qv_surface = qv_surface

        self.connection_evap.atmosphere = self.atmosphere
        self.connection_evap.q_evap = self.q_evap
        self.connection_evap.ql = self.ql_surface
        self.connection_evap.qv = self.qv_surface
        self.connection_evap.T_surface = self.Ts

    def update_transpiration(self, q_trans=[]):
        """
        Updates the flux related to iso - transpiration from each layer to current time step

        Returns
        -------
        """
        try:
            if len(q_trans) == len(self.__layers):
                self.__transpiration = q_trans
            else:
                raise ValueError("number of transpiration fluxes must be equal to that of soil layers.")

            # assign liquid fluxes to the connections
            if len(q_trans) == len(self.connections_transpiration):

                for qtrans, c_trans in zip(self.q_transpiration, self.connections_transpiration):
                    c_trans.ql_transpiration = qtrans
            else:
                return "Number of liquid fluxes must be equal to the flux connections"

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def update_precipitation(self, q_prec=0.0, c_prec={"2H": 1.0, "18O": 1.0}):

        """
        Updates the precipitation flux to current time step

        Returns
        -------
        """
        try:
            self.__q_precipitation = q_prec
            self.__c_precipitation = c_prec

            self.connection_prec.q_prec = self.q_prec
            self.connection_prec.ci_prec = self.c_prec

            # self.connection_prec.set_conc_prec(c_prec["2H"], "2H")
            # self.connection_prec.set_conc_prec(c_prec["18O"], "18O")

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    def update_runoff(self, q_runoff=0.0):

        """
        Updates the runoff flux to current time step

        Returns
        -------
        """
        try:
            self.__q_runoff = q_runoff
            self.connection_runoff.q_runoff = self.q_runoff

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def update_connection_to_aquifer(self, ql_layer=0.0):

        """
        Updates the runoff flux to current time step

        Returns
        -------
        """

        try:
            # Assign Atmospheric boundary
            if self.__connection_to_aquifer.soil_layer is self.bottom_layer:

                if self.__liquid_fluxes:
                    ql_aq = self.__liquid_fluxes[-1]  # if ql is not None
                else:
                    ql_aq = None
            else:
                ql_aq = ql_layer

            # self.__connection_to_aquifer.ql = ql_aq
            self.connection_to_aquifer.ql = ql_aq

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.")

    def update_neuman_boundary(self, q_neuman=0.0, c_neuman={"2H": 1.0, "18O": 1.0}):

        """
        Updates the neuman flux to current time step

        Returns
        -------
        """
        try:
            self.__q_neuman = q_neuman
            self.__c_neuman = c_neuman

            self.connection_neuman.q_neuman = self.q_neuman
            self.connection_neuman.ci_neuman = self.c_neuman

            # self.connection_prec.set_conc_prec(c_prec["2H"], "2H")
            # self.connection_prec.set_conc_prec(c_prec["18O"], "18O")

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    def update_dirichlet_boundary(self, c_dirichlet={"2H": 1.0, "18O": 1.0}):

        """
        Updates the neuman flux to current time step

        Returns
        -------
        """
        try:
            self.__c_dirichlet = c_dirichlet
            self.__connection_dirichlet.ci_dirichlet = self.c_dirichlet

            self.connection_dirichlet.set_conc_dirichlet(self.c_dirichlet["2H"], "2H")
            self.connection_dirichlet.set_conc_dirichlet(self.c_dirichlet["18O"], "18O")

        except ValueError:
            raise NotImplementedError("A required value was not provided.")



