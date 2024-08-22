""""

"""

# -*- coding: utf-8 -*-
from math import exp


class iso_storage(object):
    """
    Description
    ===========
      Layer storing the state, location and dimension:
      -water storage
      -isotopes
      -temperature
    Each storage belong to a cell defining the horizontal shape (area in m2)
    """

    # Class attribute
    Dict_of_all_iso_storages = {}

    def __init__(self,
                 upper_boundary,  # in (m)
                 lower_boundary,  # in (m)
                 #center,  # as LSI_Storages.Point
                 area,
                 conc_iso_liquid={"2H": 1.0, "18O": 100.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 theta=0.35,  # volumetric water content m3/m3
                 theta_sat=0.35,  # #volumetric water content m3/m3 at saturation, or porosity of the soil m3/m3
                 tortuosity=0.67,  # tortuosity of the soil m/m
                 T=283.15  # temperature in Kelvin
                 ):
        """
        Constructor of iso_storage
        """
        #flux_node.__init__(self, ID, conc_iso_liquid, T)
        #self.__class__.Dict_of_all_iso_storages[ID] = self  # add Instance of flux node to Dict_of_all_flux_nodes

        self.upper_boundary = upper_boundary
        self.lower_boundary = lower_boundary
        #assert (isinstance(new_layer, Point)), "the center of a layer must be an instance of LSI_Storages.Point"
        #self.center = center
        self.__area = area
        self.__conc_iso_liquid = conc_iso_liquid
        self.thickness = lower_boundary - upper_boundary
        self.theta = theta
        self.theta_sat = theta_sat
        self.tortuosity = tortuosity
        self.flux_connections = []
        self.T = T

    # Property: layers
    def get_area(self):
        """
        @return: Returns the horizontal area of the storage (the same as the cell it belongs to)  in m2.
        """
        return self.__area

    area = property(get_area, None, None, "Horizontal area of the storage (m)")

    def get_volume(self):
        """
        @return: Returns the volume of the storage (the same as the cell it belongs to) in m3.
        """
        return self.area * self.thickness

    volume = property(get_volume, None, None, "Volume of the storage (m3)")

    def get_liquid_H2O_volume(self):
        """
        @return: Returns the volume of the storage filled with liquids  in m3.
        """
        return self.volume * self.theta

    liquid_H2O_volume = property(get_liquid_H2O_volume, None, None, "Volume of the storage (m3)")

    def get_gas_volume(self):
        """
        @return: Returns the volume of the storage filled with gas in m3 (not only water).
        """
        return self.volume * (1 - self.theta)

    gas_volume = property(get_gas_volume, None, None, "Volume of the storage (m3)")

    def get_gas_H2O_volume(self):

        """
        @return: Returns the volume of the storage filled with gas in m3 H2O.

        var%cvsat  = esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3
        var%cv     = var%rh*var%cvsat
        """
        cv_sat = self.cv_sat  # (m**3_H20/m**3)
        rH = self.soil_relative_humidity(psi=100, T=self.T)  # relative humidity in Percent

        return cv_sat * rH * self.gas_volume

    gas_H2O_volume = property(get_gas_H2O_volume, None, None, "Volume of the storage (m3)")

    def get_storage_i(self, Isotopologue):
        """
        @return: Returns the total amount of the isotopes in kg.
        """
        return self.liquid_H2O_volume * self.get_conc_iso_liquid(
            Isotopologue) + self.gas_H2O_volume * self.get_conc_iso_vapor(Isotopologue)

    def delta_theta(self, ql):
        """
        Returns the change in theta based on the given liquid flux in m3

        @param ql: Liquid flux leaving or entering the current layer in m3
        @type ql: float
        """

        raise NotImplementedError

    def soil_relative_humidity(self, psi, T):
        """
        Calculates the relative humidity in the soil air space based on the temperature and the matrix potential (Psi)

        @param psi: Matrix potential of the soil layer
        @type psi: float

        @param T: Temperature of the soil layer in kelvin
        @type T: float
        """
        Mw = 0.018016  # molecular wt: water    (kg/mol)
        gravity = 9.80  # gravity acceleration (m/s2)
        R = 8.3144621  # gas constant(J/(mol*K))

        hr = exp(Mw * gravity * psi / (R * T))

        return hr

    def get_conc_iso_vapor(self, Isotopologue):
        """
        @return: Returns the concentration of the given isotopologue in the vapor phase in kg/m3.
        """
        return self.__conc_iso_liquid[Isotopologue] * self.alpha_i(Isotopologue)

    def get_conc_iso_liquid(self, Isotopologue):
        """
        @return: Returns the concentration of the given isotopologue in the liquid phase in kg/m3.
        """
        return self.__conc_iso_liquid[Isotopologue]

    def set_conc_iso_liquid(self, new_conc_iso_liquid, Isotopologue):
        try:
            assert new_conc_iso_liquid > 0.0, "Given isotope concentration of the isotopologue must be in kg/m2 and >0.0"

            self.__conc_iso_liquid[Isotopologue] = new_conc_iso_liquid

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def get_cv_sat(self):
        """
        Returns the saturated water vapor volumetric mass (m**3_H20/m**3) .

        T_atmosphere = Temperature in Kelvin

        pv_sat = saturated water vapor volumetric mass (m**3_H20/m**3)

        TODO: Check for which temperature range it is valid

        Control status: Checked on 15.05.2013 --> Results are the same as for SLI
        """

        # SLI: sli_main.f90::L205
        cv_sat = 0.002166 * self.p_sat / (self.T + 273.16) * 1000  # m3/m3

        return cv_sat

    cv_sat = property(get_cv_sat, None, None,
                      "Returns the saturated water vapor volumetric mass (kg/m**3) under free air")

    def get_p_sat(self):
        """
        Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch and valid over the range 0 to 60 Degree celsius

        T_atmosphere = Temperature in Kelvin

        p_sat = Saturation vapour pressure, ps, in pascal

        SLI: cable_sli_utils.f90 Lines: 2513-2525 <--checked and own implementation is okay for temp. above 273.15 Kelvin
        """

        # SLI: cable_sli_utils.f90::L2513-2525

        p_sat = 610.6 * exp(
            17.27 * (self.T - 273.15) / ((self.T - 273.15) + 237.3))  # Saturation vapour pressure, ps, in pascal

        if self.T < 273.16:  # below 0 Degree celsius
            p_sat = -4.86 + 0.855 * p_sat + 0.000244 * p_sat ** 2  # ps for ice

        return p_sat

    p_sat = property(get_p_sat, None, None,
                     "Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch and valid over the range 0 to 60 Degree celsius")


    def alpha_i(self, Isotopologue, T=None, ignore_alpha_i=False):
        """
        Calculates the liquid-vapour isotopic frctionation factor at equilibrium given by Majoube (1971) as a function of temperature T (K)

        Description
        ===========

        @param Isotopologue: Identifier for the solute to be calculated (currently supported '2H' and '18O')
        @type Isotopologue: Sting

        @param T:  Temperature in Kelvin. Per default the temparature of the flux node is used
        @type T: Float

        @param ignore_alpha_i:  If set to "True" isotopic frctionation factor will be set to 1 = no fractionation
        @type ignore_alpha_i: Boolean

        Control status: Checked on 16.05.2013 --> Results are the same as for SLI
        """
        try:
            # SLI: cable_sli_solve.f90::L2275 - L2305
            if T == None:
                T = self.T

            if ignore_alpha_i == False:
                if "18O" == Isotopologue:
                    alpha_i = exp(-(24844 / T ** 2 + (-76.248) / T + 0.052612))
                elif "2H" == Isotopologue:
                    alpha_i = exp(-(1137 / T ** 2 + (-0.4156) / T - 0.0020667))
                else:
                    raise NotImplementedError
                return alpha_i
            else:
                return 1.0

        except ValueError as err:
            print(err)
            raise NotImplementedError
