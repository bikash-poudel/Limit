'''
Created on 11.12.2024
@author: poudel-b
'''
# -*- coding: utf-8 -*-

from . import storage
from . import fluxes


class cell(object):

    def __init__(self, atmosphere):

        self.__atmosphere = atmosphere
        self.__layers = []
        self.__connection_water = []
        self.__connection_heat = []
        self.__connection_vapor = []
        self.__connection_evap = None

        self.__hleft = None
        self.__hright = None
        self.__Tleft = None
        self.__Tright = None

        self.__qevap = None
        self.__ql = None
        self.__qv = None

    @property
    def atmosphere(self):
        return self.__atmosphere

    @property
    def layers(self):

        return self.__layers

    @property
    def hleft(self):
        return self.__hleft

    @hleft.setter
    def hleft(self, h):
        self.__hleft = h

    @property
    def hright(self):
        return self.__hright

    @hright.setter
    def hright(self, h):
        self.__hright = h

    @property
    def Tleft(self):
        return self.__Tleft

    @Tleft.setter
    def Tleft(self, T):
        self.__Tleft = T

    @property
    def Tright(self):
        return self.__Tright

    @Tright.setter
    def Tright(self, T):
        self.__Tright = T

    @property
    def q_evap(self):
        return self.__qevap

    @q_evap.setter
    def q_evap(self, q):
        self.__qevap = q

    @property
    def liquid_fluxes(self):
        return self.__ql

    @liquid_fluxes.setter
    def liquid_fluxes(self, q):
        self.__ql = q

    @property
    def vapor_fluxes(self):
        return self.__qv

    @vapor_fluxes.setter
    def vapor_fluxes(self, q):
        self.__qv = q

    @property
    def water_connections(self):
        return self.__connection_water

    @property
    def heat_connections(self):
        return self.__connection_heat

    @property
    def vapor_connections(self):
        return self.__connection_vapor

    @property
    def storage_connections(self):
        return self.water_connections + self.heat_connections + self.vapor_connections

    @property
    def connection_evap(self):
        return self.__connection_evap

    @property
    def flux_connections(self):
        return self.storage_connections + self.connection_evap

    """Functions"""

    def add_layer(self, new_layer):
        try:
            assert (isinstance(new_layer, storage.soil_layer)), \
                "new_layer must be an instance of layer.iso_soil_layer"

            self.__layers.append(new_layer)

            self.__layers.sort(key=lambda
                x: x.lower_boundary)  # sort them to make sure that they are in the right order from top to lowest

        except AssertionError:
            raise NotImplementedError

    def install_connections(self, water=True, heat=True, vapor=True):

        """"
        Registers all the relevant connections within the layers
        """

        try:
            for left_node, right_node in zip(self.layers[:-1], self.layers[1:]):

                # Define connections and assign to each layers
                if water:
                    w = fluxes.water_flux(left_node=left_node, right_node=right_node)
                    self.__connection_water.append(w)
                if heat:
                    h = fluxes.heat_flux(left_node=left_node, right_node=right_node)
                    self.__connection_heat.append(h)
                if vapor:
                    v = fluxes.vapor_flux(left_node=left_node, right_node=right_node)
                    self.__connection_vapor.append(v)

        except Exception:
            raise NotImplementedError("An unexpected error occurred.")

    def add_evaporation(self, top_layer):
        """
        Assign boundary connection eq: atmosphere etc..

        Returns
        -------
        """
        try:
            # Assign Atmospheric boundary

            ev = fluxes.evaporation(atmosphere=self.atmosphere, soil_layer=top_layer)
            self.__connection_evap = ev

        except ValueError as err:
            raise NotImplementedError("A required value was not provided.") from err

    """Update_storages"""

    def update_head(self, head):

        try:

            for lr, h in zip(self.__layers, head):

                lr.head = h

        except:
            raise NotImplementedError

    def update_temperature(self, T):

        try:

            for lr, t in zip(self.__layers, T):

                lr.T = t

        except:
            raise NotImplementedError


    """update_boundaries"""

    def update_boundary_head(self):

        l_left, l_right = self.layers[0], self.layers[-1]
        #l_left.head = self.hleft
        l_right.head = self.hright

    def update_boundary_T(self, Tleft=25, Tright=25):

        l_left, l_right = self.layers[0], self.layers[-1]
        #l_left.T = self.Tleft
        l_right.T = self.Tright

    def update_evaporation(self, qevap):

        self.q_evap = qevap
        self.__connection_evap.q_evap = self.q_evap

