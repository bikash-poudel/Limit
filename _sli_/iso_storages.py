'''
Created on 23.05.2024
@author: poudel-b
'''
import math
# -*- coding: utf-8 -*-

from math import exp


class flux_node(object):
    """
    Description
    ===========
      A flux node can be a boundary condition a Layer storing the state, location and dimension or the atmosphere.

    """

    # Class attribute
    Dict_of_all_flux_nodes = {}
    current_ID = 1

    def __init__(self, conc_iso_liquid={"2H": 1.0, "18O": 1.0}, T=283.15):
        """
        Constructor of flux_node

        @param conc_iso_liquid: Concentration of isotope species (kg/m3) in the liquid phase {"2H" : 1.0,"18O" : 1.0}
        @type conc_iso_liquid: Dict

        @param T: Temperature in Kelvin
        @type T: float
        """
        self.__ID = self.__class__.current_ID
        self.__class__.current_ID += 1  # increment the current ID for the next flux node to be created
        self.__class__.Dict_of_all_flux_nodes[self.__ID] = self  # add Instance of flux node to Dict_of_all_flux_nodes

        self.__Connections = []  # list holding all flux connections from and to this layer (e.g. advection, diffusion of the liquid or vapor phase, or any boundary condition)
        self.__Connections_to_iso_storages = []
        self.__Connections_to_boundaries = []
        self.__Connections_to_right = []
        self.__Connections_to_left = []
        self.__conc_iso_liquid = conc_iso_liquid
        self.T = T  # temperature in Kelvin

    def get_connections(self):
        return self.__Connections

    connections = property(get_connections, None, None, "List of all flux connections connected with this node")

    def get_connections_to_iso_storages(self):
        return self.__Connections_to_iso_storages

    connections_to_iso_storages = property(get_connections_to_iso_storages, None, None,
                                           "List of all flux connections between this node and iso storages")

    def get_connections_to_boundaries(self):
        return self.__Connections_to_boundaries

    connections_to_boundaries = property(get_connections_to_boundaries, None, None,
                                         "List of all flux connections between this node and boundaries (like atmopshere)")

    def get_connections_to_left(self):
        return self.__Connections_to_left

    connections_to_left = property(get_connections_to_left, None, None,
                                   "List of all flux connections from this node to left node if available")

    def get_connections_to_right(self):
        return self.__Connections_to_right

    connections_to_right = property(get_connections_to_right, None, None,
                                    "List of all flux connections from this node to right node if available")

    def RegisterConnection(self, newConnection):
        """
        Registers the given connection.
        """
        self.__Connections.append(newConnection)

        if newConnection.left_node == self:
            if isinstance(newConnection.right_node, iso_storage):
                # left connection only to the storages, to ta boundaries are connection to the boundaries
                self.__Connections_to_right.append(newConnection)
            other_node = newConnection.right_node
        elif newConnection.right_node == self:
            if isinstance(newConnection.left_node, iso_storage):
                # left connection only to the storages, to ta boundaries are connection to the boundaries
                self.__Connections_to_left.append(newConnection)
            other_node = newConnection.left_node

        else:
            raise NotImplementedError

        if isinstance(other_node, iso_storage):
            self.__Connections_to_iso_storages.append(newConnection)
        elif isinstance(other_node, iso_atmosphere) or isinstance(other_node, iso_aquifer):
            self.__Connections_to_boundaries.append(newConnection)

    def DeregisterConnection(self, oldConnection):
        """
        Deregisters the given connection.
        """
        if oldConnection in self.__Connections:
            self.__Connections.remove(oldConnection)
            if oldConnection in self.__Connections_to_iso_storages:
                self.__Connections_to_iso_storages.remove(oldConnection)
            elif oldConnection in self.__Connections_to_boundaries:
                self.__Connections_to_boundaries.remove(oldConnection)
            else:
                raise NotImplementedError

            return True
        else:
            return False

    """Functions"""

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

    def set_conc_iso_liquid(self, new_c_iso_liquid, Isotopologue):

        try:
            assert new_c_iso_liquid >= 0.0, "Given isotope concentration of the isotopologue must be in kg/m2 and >0.0"

            self.__conc_iso_liquid[Isotopologue] = new_c_iso_liquid

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def alpha_i(self, Isotopologue, T=None, **kwargs):
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
        ignore_alpha_i = kwargs.get('ignorealphai', False)
        try:
            # SLI: cable_sli_solve.f90::L2275 - L2305
            if T is None:
                T = self.T

            if not ignore_alpha_i:
                if "2H" == Isotopologue:
                    alpha_i = exp(-(24844 / T ** 2 + (-76.248) / T + 0.052612))
                elif "18O" == Isotopologue:
                    alpha_i = exp(-(1137 / T ** 2 + (-0.4156) / T - 0.0020667))
                else:
                    raise NotImplementedError
                return alpha_i
            else:
                return 1.0

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def get_cv_sat(self):
        """
        Returns the saturated water vapor volumetric mass (m**3_H20/m**3) .

        T_atmosphere = Temperature in Kelvin

        pv_sat = saturated water vapor volumetric mass (m**3_H20/m**3)

        TODO: Check for which temperature range it is valid
        TODO: recheck the result with SLI
        Control status: Checked on 03.06.2024 --> Results are the same as for SLI
        """

        # SLI: sli_utils.f90::L1416
        Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
        R = 8.3142995834350586  # universal gas constant (j / mol / k)

        cv_sat = self.p_sat * Mw / R / self.T / 1000  # m3/m3
        # cv_sat = 0.002166 * self.p_sat / self.T * 1000  # m3/m3

        return cv_sat

    cv_sat = property(get_cv_sat, None, None,
                      "Returns the saturated water vapor volumetric mass (kg/m**3) under free air")

    def get_p_sat(self):
        """
        Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch and valid over the range 0 to 60 Degree celsius

        T_atmosphere = Temperature in Kelvin

        p_sat = Saturation vapour pressure, ps, in pascal

        <--checked and own implementation is okay for temp. above 273.15 Kelvin
        """
        # SLI: cable_sli_utils.f90::L2195-2195
        Tzero = 273.16000366210938
        p_sat = 610.59999465942383 * exp(
            17.270000457763672 * (self.T - Tzero) /
            ((self.T - Tzero) + 237.30000305175781))  # Saturation vapour pressure, ps, in pascal

        if self.T < 273.16:  # below 0 Degree celsius
            p_sat = -4.86 + 0.855 * p_sat + 0.000244 * p_sat ** 2  # ps for ice

        return p_sat

    p_sat = property(get_p_sat, None, None,
                     "Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch "
                     "and valid over the range 0 to 60 Degree celsius")

    @classmethod
    def delta_to_concentration(self, delta_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}, M_w=0.018,
                               M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
        """
        Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

        delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
        solute_i = sting identifying the isotope species currently supported '2H' or '18O'
        R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]
        M_w = molar mass of water (kg)  0.018 kg for H_2 O
        M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
        mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

        SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
        """
        R_i = self.delta_to_ratio(delta_i, solute_i, R_ref)  # isotopic ratio of isotopic species i in sample
        c_solutes = R_i * density_H2O * M_i[solute_i] / M_w
        return c_solutes

    @classmethod
    def delta_to_ratio(self, delta_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}):
        """
        Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

        delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
        solute_i = sting identifying the isotope species currently supported '2H' or '18O'
        R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]

        SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
        """
        R_i = ((delta_i) / 1000 + 1) * R_ref[solute_i]  # isotopic ratio of isotopic species i in sample
        return R_i

    @classmethod
    def ratio_to_delta(self, R_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}):
        """
        Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

        delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
        solute_i = sting identifying the isotope species currently supported '2H' or '18O'
        R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]

        SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
        """
        delta_i = (R_i / R_ref[solute_i] - 1) * 1000
        return delta_i

    @classmethod
    def concentration_to_delta(self, c_solute_i, solute_i, R_ref={"2H": 0.00015576, "18O": 0.00200520}, M_w=0.018,
                               M_i={"2H": 0.019, "18O": 0.020}, density_H2O=1000.0):
        """
        Converts the given delta signature of isotope species i into the concentration of isotope species i in kg/m**3

        delta_i = delta signature in permil of the isotope species i relative to the V-SMOW standard
        solute_i = sting identifying the isotope species currently supported '2H' or '18O'
        R_ref = isotopic ratio of isotopic species i in V-SMOW standard (Gonfiantini, 1978) [V-SMOW_2H/1H = 0.00015576 and V-SMOW_18O/16O = 0.00200520]
        M_w = molar mass of water (kg)  0.018 kg for H_2 O
        M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
        mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

        SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
        """
        R_i = self.concentration_to_ratio(c_solute_i, solute_i, M_w, M_i,
                                          density_H2O)  # isotopic ratio of isotopic species i

        delta_i = self.ratio_to_delta(R_i, solute_i, R_ref)

        return delta_i

    @classmethod
    def concentration_to_ratio(self, c_solutes_i, solute_i, M_w=0.018, M_i={"2H": 0.019, "18O": 0.020},
                               density_H2O=1000.0):
        """
        Converts the given concentration of isotope species i  (kg/m**3) into the ratio of isotope species i

        solute_i = sting identifying the isotope species currently supported '2H' or '18O'
        M_w = molar mass of water (kg)  0.018 kg for H_2 O
        M_i=molar mass of isotopic species i (kg) 0.019 kg for (_^2)H_2 O and 0.020 kg for H_2 (_^18)O.
        mass_density_H2O = Density/liquid volumetric mass of water (kg/m^3 )  1000.0

        SLI: sli_init.f90 Lines: 115 -142 <--checked and own implementation is okay
        """

        R_i = c_solutes_i / (density_H2O * M_i[solute_i] / M_w)

        return R_i


class iso_atmosphere(flux_node):
    """
    Layer storing isotopes

    Description
    ===========
      The
    """

    def __init__(self,
                 conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                 conc_iso_vapor={"2H": 1.0, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 T=283.15,  # temperature in the atmosphere in Kelvin
                 Rh_atmosphere=0.20,  # relative humidity of the atmosphere (-)
                 Pa_atmosphere=10 ** 5,  # atmospheric pressure (Pa)
                 initial_wind_speed=2.0,  # wind speed at top of canopy in m/s
                 wind_speed=2.0,  # wind speed at top of canopy in m/s
                 hc=10,  # canopy height [m] (e.g. 40)
                 d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                 z0m=0.1 * 10,  # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                 LAI=1.1,  # leaf area index (e.g. 2.0)
                 extku=1.5):  # extinction coeff't for windspeed (e.g. 0.5)

        """
        Constructor of LSI_iso_cell

        @param conc_iso_vapor: Concentration of isotope species (kg/m3) in the vapor phase {"2H" : 1.0,"18O" : 1.0}
        @type conc_iso_vapor: Dict

        @param T: atmospheric temperature in (K)
        @type T: float

        @param Rh_atmosphere: relative humidity of the atmosphere (-)
        @type Rh_atmosphere: float

        @param Pa_atmosphere: atmospheric pressure (Pa)
        @type Pa_atmosphere: float

        @param wind_speed: wind speed at top of canopy (m/s)
        @type wind_speed: float

        @param hc: canopy height [m] (e.g. 40)
        @type hc: float

        @param d0: displacement height (e.g. 0.7)
        @type d0: float

        @param z0m: roughness height for momentum (e.g. 0.1)
        @type z0m: float

        @param LAI: leaf area index (e.g. 2.0)
        @type LAI: float

        @param extku: extinction coeff't for windspeed (e.g. 0.5)
        @type extku: float
        """
        flux_node.__init__(self, conc_iso_liquid, T)
        self.__conc_iso_vapor = conc_iso_vapor
        self.Rh = Rh_atmosphere  # porosity of the soil m3/m3
        self.Pa = Pa_atmosphere  # Atmospheric preassure (Pa)
        self.wind_speed = wind_speed  # wind speed at top of canopy in m/s
        self.hc = hc  # canopy height [m]
        self.d0 = d0  # displacement height
        self.z0m = z0m  # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
        self.LAI = LAI  # leaf area index (e.g. 2.0)
        self.extku = extku  # extinction coeff't for windspeed

    def get_conc_iso_vapor(self, Isotopologue):
        """
        @return: Returns the concentration of the given isotopologue in the vapor phase in kg/m3.
        """
        return self.__conc_iso_vapor[Isotopologue]

    def set_conc_iso_vapor(self, new_conc_iso_vapor, Isotopologue):
        try:
            assert new_conc_iso_vapor > 0.0, "Given isotope concentration of the isotopologue must be in kg/m2 and >0.0"

            self.__conc_iso_vapor[Isotopologue] = new_conc_iso_vapor

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def get_cv_a(self):
        """
        Converts the Rh to concentration of the vapor in the air

        T = Temperature in Kelvin
        Rh = Relative humidity

        Control status: Checked on 15.05.2013 --> Results are the same as for SLIsli_main.f90::L205
        """

        # SLI: sli_main.f90::L205
        # vmet%cva  = vmet(:)%rha * esat(vmet(:)%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet(:)%Ta+Tzero)

        cv_a = self.Rh * self.p_sat * 0.018 / (1000 * 8.314 * self.T)

        return cv_a

    cv_a = property(get_cv_a, None, None, "Returns concentration of the vapor in the air ")

    def civ_a(self, Isotopologue):
        """
        Returns the concentration of isotope i in the atmosphere with respect to the current relative humidity
        """

        # SLI: sli_main.f90::L206
        # vmet%civa = Ra * vmet%cva * rhow * Mi/Mw18  ! kg HDO (liq) m-3 (air) !JLU David: civa = concentration of isotope in the atmosphere

        civ_a = self.cv_a * self.get_conc_iso_vapor(Isotopologue)

        return civ_a


class iso_aquifer(flux_node):
    """"
    Storage of isotopes in boundary can be an instance of flux node.
    """

    def __init__(self,
                 conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                 T=283.15,  # temperature in Kelvin
                 ):
        flux_node.__init__(self, conc_iso_liquid, T)


class iso_storage(flux_node):
    """
    Description
    ===========
      Layer storing the state, location and dimension:
      -water storage
      -isotopes
      -temperature
    Each storage belong to a cell defining the horizontal shape (area in m2)
    """

    Dict_of_all_iso_storages = {}

    def __init__(self,
                 ID,
                 upper_boundary,  # in (m)
                 lower_boundary,  # in (m)
                 conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 theta_t0=0.3,  # volumetric water content m3/m3 at previous time
                 theta=0.3,  # volumetric water content m3/m3
                 theta_0=0.0,  # volumetric water content m3/m3 at high suctions
                 theta_sat=0.35,  # volumetric water content m3/m3 at saturation, or porosity of the soil m3/m3
                 T0=283.15,  # Temperature previous time [kelvin]
                 T=283.15,  # temperature current time [Kelvin]
                 rH=0,  # soil relative humidity
                 psi_0=None,  # soil matric potential at previous time
                 psi=None  # matric potential at current time
                 ):
        """
        Constructor of iso_storage
        """
        flux_node.__init__(self, conc_iso_liquid, T)
        self.__ID = ID
        self.__class__.Dict_of_all_iso_storages[self.__ID] = self  # add Instance of flux node to Dict_of_all_flux_nodes
        self.__cell = None
        self.upper_boundary = upper_boundary
        self.lower_boundary = lower_boundary
        self.thickness = lower_boundary - upper_boundary
        self.theta_t0 = theta_t0
        self.theta = theta
        self.theta_sat = theta_sat
        self.theta_0 = theta_0
        self.T0 = T0
        self.T = T
        self.rH = rH
        self.psi_0 = psi_0
        self.psi = psi
        self.pond = None

    """Properties"""

    @property
    def get_ID(self):
        return self.__ID

    @property
    def cell(self):
        return self.__cell

    @property
    def area(self):
        return self.__cell.area

    @property
    def center(self):
        """
        @return: Returns the 3D coordinates of center of iso_storage / iso_layer
        """
        x = self.cell.location.x  # x coordinate of center of storage
        y = self.cell.location.y  # y coordinate of center of storage

        center_depth = (self.upper_boundary + self.lower_boundary) / 2
        return Point(x, y, center_depth)

    @cell.setter
    def cell(self, c):

        """
        @return: sets the cell for the current layer / storage
        """
        self.__cell = c

    @property
    def dT(self):
        return self.T - self.T0

    """Functions"""

    def get_storage_i(self, Isotopologue, delta_cv=0.0, deltaS_liq=0.0, **kwargs):
        """
        @return: Returns the total amount of the isotopes
        """
        storage_v = self.get_liquid_volume()
        storage_i = self.get_conc_iso_liquid(Isotopologue=Isotopologue) \
                    * self.del_eff_saturation(Isotopologue=Isotopologue, d_cv=delta_cv, d_sliq=deltaS_liq, **kwargs)
        return storage_v * storage_i

    def get_eff_liquid_volume(self, Isotopologue, **kwargs):
        """
        @return: Returns the effective volume of the storage filled with liquids  in m3.
        """
        return self.eff_saturation(Isotopologue=Isotopologue, **kwargs) * (self.theta_sat - self.theta_0) * self.thickness

    def get_liquid_volume(self):

        """
        @return: Returns the effective volume of the storage filled with liquids  in m3.
        """
        return (self.theta_sat - self.theta_0) * self.thickness

    def get_saturation(self):

        """
        Returns the degree of saturation (effective saturation) according to (Brooks and Corey, 1964)

        Description
        ===========

        @param theta: liquid phase (m**3/m**3 )
        @type theta: Float

        @param theta_sat: theta at saturation equal to soil porosity (m**3/m**3 )
        @type theta_sat: Float

        @param theta_0: defined residual water content (water remaining at high suctions) (m**3/m**3 ) --> should be 0
        @type theta_0: Float

        """
        # SLI: sli_main.f90::L275
        S = (self.theta - self.theta_0) / (self.theta_sat - self.theta_0)

        return S

    S = property(get_saturation, None, None, "Returns effective saturation")

    def get_liq_saturation(self):

        """
        Returns the degree of saturation in liquid phase only (without vapour phase)

        Description
        ===========

        @param S: degree of saturation (m**3/m**3 )
        @type S: Float

        """
        # SLI: sli_solve.f90::L1788,
        liq_S = (self.S - self.cv) / (1 - self.cv)

        return liq_S

    @property
    def Sl(self):
        """
        Returns the liquid saturation considering the presence of a pond.
        If a pond is present, its height is factored into the liquid saturation.
        """
        if self.pond is not None:
            return self.get_liq_saturation() + self.pond.pond_height / self.thickness / (self.theta_sat - self.theta_0)

        return self.get_liq_saturation()

    def get_cv(self):
        """
        Returns the saturated water vapor volumetric mass (m**3_H20/m**3) .


        rH = relative humidity in the soil air space based on the temperature and the matrix potential (Psi)

        c_sat = saturated water vapor volumetric mass (m**3_H20/m**3)

        TODO: Check for which temperature range it is valid
        TODO: recheck the result with SLI
        Control status: Checked on 03.06.2024 --> Results are the same as for SLI
        """
        """if self.psi is None:
            cv = self.rH * self.cv_sat
        elif self.rH is None:
            cv = self.relative_humidity(psi=self.psi, T=self.T) * self.cv_sat  # sli_utils.f90:L1417
        else:
            raise ValueError("Either of matric potential or relative humidity do not match should be not None")"""

        cv = self.relative_humidity(psi=self.psi, T=self.T) * self.cv_sat  # sli_utils.f90:L1417

        return cv

    cv = property(get_cv, None, None,
                  "Returns the saturated water vapor volumetric mass (kg/m**3) under free air")

    def eff_saturation(self, Isotopologue, **kwargs):
        try:
            s_liq = self.get_liq_saturation()
            if s_liq < 1:
                beta = self.beta(Isotopologue=Isotopologue, **kwargs)

                eff_S = self.Sl + self.cv * beta - self.cv * s_liq * beta \
                        + (self.theta_0 / self.theta_sat)
            else:
                eff_S = self.Sl + self.theta_0 / self.theta_sat

            return eff_S

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def del_eff_saturation(self, Isotopologue, **kwargs):

        try:

            Tzero = 273.16000366210938
            Mw = 0.018015999346971512  # Molecular weight of water (kg / mol)
            R = 8.3142995834350586  # universal gas constant (j / mol / k)

            p_sat_0 = 610.59999465942383 * exp(
                17.270000457763672 * (self.T0 - Tzero) /
                ((self.T0 - Tzero) + 237.30000305175781))  # Saturation vapour pressure, ps, in pascal
            cv_sat_0 = p_sat_0 * Mw / R / self.T0 / 1000  # m3/m3
            cv_0 = self.relative_humidity(psi=self.psi_0, T=self.T0) * cv_sat_0  # sli_utils.f90:L1417

            S0 = (self.theta_t0 - self.theta_0) / (self.theta_sat - self.theta_0)
            liq_S_0 = (S0 - cv_0) / (1 - cv_0)

            if self.pond is not None:
                Sl_0 = liq_S_0 + self.pond.pond_height_t0 / self.thickness / (self.theta_sat - self.theta_0)
            else:
                Sl_0 = liq_S_0

            d_sliq = self.Sl - Sl_0
            d_cv = self.cv - cv_0

            dT = self.dT
            dbeta = self.d_beta(Isotopologue=Isotopologue, T=self.T, dT=dT, **kwargs)

            s_liq = self.get_liq_saturation()
            beta = self.beta(Isotopologue=Isotopologue, **kwargs)

            if s_liq < 1:

                d_eff_S = d_sliq + beta * d_cv + self.cv * dbeta \
                          - (s_liq * beta * d_cv + self.cv * beta * d_sliq + self.cv * s_liq * dbeta)
            else:
                d_eff_S = d_sliq

            return d_eff_S

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def relative_humidity(self, psi, T):
        """
        Calculates the relative humidity in the soil air space based on the temperature and the matrix potential (Psi)

        @param psi: Matrix potential of the soil layer
        @type psi: float

        @param T: Temperature of the soil layer in kelvin
        @type T: float
        """
        # SLI_utils: L1409

        Mw = 0.018015999346971512  # molecular wt: water    (kg/mol)
        gravity = 9.8000001907348633  # gravity acceleration (m/s2)
        R = 8.3142995834350586  # gas constant(J/(mol*K))
        rhmin = 5.0000000745058060E-002  # min relative humidity SLI_utils:L1409

        if self.S < 1:
            hr = max(exp(Mw * gravity * psi / (R * T)), rhmin)
        else:
            hr = 1

        return hr

    def beta(self, Isotopologue, **kwargs):
        """"

        """
        beta = self.alpha_i(Isotopologue=Isotopologue, T=self.T, **kwargs) \
               + self.d_beta(Isotopologue=Isotopologue, T=self.T, dT=self.dT, **kwargs)

        return beta

    def d_beta(self, Isotopologue, T=None, dT=None, **kwargs):
        """
        Calculates the factor relating liquid and vapour isotope concentration for species i (-)

        Description
        ===========

        @param T:  liquid-vapour isotopic frctionation factor at equilibrium given by Majoube (1971) as a function of temperature T (K) at soil depth
        @type T: Float

        @ param dT: change in temperature in current time
        @type T: Float

        @param solute: Identifier for the solute to be calculated (currently supported '2H' and '18O')
        @type solute: Sting

        @param ignore_alpha_i:  If set to "True" isotopic frctionation factor will be set to 1 = no fractionation
        @type ignore_alpha_i: Boolean

        Control status: Checked on 16.05.2013 --> Results are the same as for SLI
        """
        ignore_alpha_i = kwargs.get('ignorealphai', False)

        try:
            # SLI: cable_sli_solve.f90::L2275 - L2305
            if T is None:
                T = self.T

            if dT is None:
                dT = self.dT

            if not ignore_alpha_i:
                if "2H" == Isotopologue:
                    d_alpha = (2 * 24844 / T ** 3 + (-76.248) / T ** 2) / exp(
                        24844 / T ** 2 + (-76.248) / T + 0.052612)
                elif "18O" == Isotopologue:
                    d_alpha = (2 * 1137 / T ** 3 + (-0.4156) / T ** 2) / exp(
                        1137 / T ** 2 + (-0.4156) / T - 0.0020667)
                else:
                    raise NotImplementedError

            else:
                d_alpha = 0

            delta_beta = d_alpha * dT

            return delta_beta
        except ValueError as err:
            print(err)
            raise NotImplementedError

    def add_pond(self, pond):
        """
        Adds iso_pond instance to the iso_storage.

        @param pond: An instance of iso_pond to be added to the storage.
        """
        self.pond = pond


class iso_pond(flux_node):
    """
    Description
    ===========
      A pond like layer storing the state, location and dimension:
      -water storage
      -isotopes
      -temperature
    Each storage belong to a cell defining the horizontal shape (area in m2)
    """

    def __init__(self,
                 pond_height_t0=0.0,
                 pond_height=0.0,
                 conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                 T=283.15,  # temperature in Kelvin
                 ):
        """
        Constructor of iso_storage
        """
        flux_node.__init__(self, conc_iso_liquid, T)
        self.pond_height_t0 = pond_height_t0
        self.__pond_height = pond_height
        self.delta_h0 = self.__pond_height - pond_height_t0  # change in pond height

    @property
    def pond_height(self):
        return self.__pond_height

    @pond_height.setter
    def pond_height(self, pond_height):
        self.__pond_height = pond_height


class iso_soil_layer(iso_storage):
    """
    Description
    ===========
    Soil Layer storing the state of the:
    -water storage
    -isotopes
    -temperature

    and the physical parameters of the soil:
    -tortuosity
    -theta_0
    """

    def __init__(self,
                 ID,
                 upper_boundary,  # in (m)
                 lower_boundary,  # in (m)
                 # cell,  # as iso_cell
                 conc_iso_liquid={"2H": 1.0, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 theta_t0=0.3, # volumetric water content m3/m3 at previous time
                 theta=0.3,  # volumetric water content m3/m3 at current time
                 theta_0=0.0,  # volumetric water content m3/m3 at high suctions
                 theta_sat=0.35,  # volumetric water content m3/m3 at saturation, or porosity of the soil m3/m3
                 tortuosity=0.67,  # tortuosity of the soil m/m
                 T0 = 283.15,  # temperature at previous time step
                 T=283.15,  # temperature in Kelvin at current time
                 rH=0,  # soil relative humidity
                 psi_0=1,  # matric potential for previous time
                 psi=1,  # matric potential current time
                 ):
        iso_storage.__init__(self, ID, upper_boundary, lower_boundary, conc_iso_liquid,
                             theta_t0, theta, theta_0, theta_sat, T0, T, rH, psi_0, psi)
        flux_node.__init__(self, conc_iso_liquid, T)
        self.tortuosity = tortuosity
        self.rH = rH


class Point(object):
    """
    A point in space.

    Description
    ===========
      A point with x, y, and z-coordinates.
    """

    def __init__(self, x, y, z=None):
        """

        @param x: x coordinate of this point.
        @type x:  Numeric

        @param y: y coordinate of this point.
        @type y:  Numeric

        @param z: z coordinate of this point.
        @type z:  Numeric

        @return: Point
        @rtype: map.Map.Point
        """
        self.__x = x
        self.__y = y
        self.__z = z

    def __repr__(self):
        """
        @return: Returns a string identifying the point in space.
        @rtype: String
        """
        return_Str = "Point: x=" + str(self.x) + " y=" + str(self.y) + " z=" + str(self.z)
        return return_Str

    """Properties"""

    @property
    def x(self):
        """


        @return: Returns the x coordinate of this point.
        @rtype:  Numeric
        """
        return self.__x

    @property
    def y(self):
        """
        "The y coordinate of the Point as Numeric"

        @return: Returns the y coordinate of this point.
        @rtype:  Numeric
        """
        return self.__y

    @property
    def z(self):
        """
        @return: Returns the z coordinate of this point.
        @rtype:  Numeric
        """
        return self.__z

    @x.setter
    def x(self, new_x):
        """
        @param new_x: Sets the x coordinate of this point.
        @type new_x:  Numeric
        """
        try:
            # checks new_x
            assert (isinstance(new_x, (int, float))), "Value must be Numeric"

            self.__x = new_x

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    @y.setter
    def y(self, new_y):
        """
        @param new_y: Sets the y coordinate of this point.
        @type new_y:  Numeric
        """
        try:
            # checks new_y
            assert (isinstance(new_y, (int, float))), "Value must be Numeric"

            self.__y = new_y

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    @z.setter
    def z(self, new_z):
        """
        "The z coordinate of the Point as Numeric"

        @param new_z: Sets the z coordinate of this point.
        @type new_z:  Numeric
        """
        try:
            # checks new_z
            assert (isinstance(new_z, (int, float))), "Value must be Numeric"

            self.__z = new_z

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    """Functions"""

    def distance(self, to_Point):
        """
        Calculates the distance between two points.

        Description
        ===========
            Returns the linear (Euclidean)distance to the to_Point as Units.Distance!

        @param to_Point: Point in question.
        @type to_Point: map.Map.Point
        """
        try:
            assert (isinstance(to_Point, Point)), "Value must be an instance of LSI_storages.Point"

            if self.z is not None and to_Point.z is not None:
                point_1 = (self.x, self.y, self.z)
                point_2 = (to_Point.x, to_Point.y, to_Point.z)

                distance = math.dist(point_1, point_2)
                # distance = math.sqrt(math.pow(self.x - to_Point.x, 2) + math.pow(self.y - to_Point.y, 2)
                #                    + math.pow(self.z - to_Point.z, 2))
            else:

                point_1 = (self.x, self.y)
                point_2 = (to_Point.x, to_Point.y)

                distance = math.dist(point_1, point_2)
                # distance = math.sqrt(math.pow(self.x - to_Point.x, 2) + math.pow(self.y - to_Point.y, 2)
                #                   + math.pow(self.z - to_Point.z, 2))

            return distance

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def __add__(self, other):
        """
        Returns self + other.

        @param other: Value of other will be added to self.
        @type other: map.Map.Point

        @return: The sum of self and other.
        @rtype: map.Map.Point
        """
        if isinstance(other, Point):
            my_x = self.x + other.x
            my_y = self.y + other.y
            if self.z is not None and other.z is not None:
                my_z = self.z + other.z
            else:
                my_z = None
            return Point(x=my_x, y=my_y, z=my_z)
        else:
            raise TypeError('Point must be added with Point')

    def __iadd__(self, other):
        """
        Returns self += other.

        @param other: Value of other will be added to self.
        @type other: map.Map.Point

        @return: self
        @rtype: map.Map.Point
        """
        if isinstance(other, Point):
            self.x += other.x
            self.y += other.y
            if self.z is not None and other.z is not None:
                self.z += other.z
            else:
                self.z = None
            return self
        else:
            raise TypeError('Point must be added with Point')

    def __sub__(self, other):
        """
        Returns self - other.

        @param other: Value of other will be subtracted from self.
        @type other: map.Map.Point

        @return: The subtraction of self and other.
        @rtype: map.Map.Point
        """
        if isinstance(other, Point):
            my_x = self.x - other.x
            my_y = self.y - other.y
            if self.z is not None and other.z is not None:
                my_z = self.z - other.z
            else:
                my_z = None
            return Point(x=my_x, y=my_y, z=my_z)
        else:
            raise TypeError('Point must be subtracted with Point')

    def __isub__(self, other):
        """
        Returns self -= other.

        @param other: Value of other will be subtracted from self.
        @type other: map.Map.Point

        @return: self
        @rtype: map.Map.Point
        """
        if isinstance(other, Point):
            self.x -= other.x
            self.y -= other.y
            if self.z is not None and other.z is not None:
                self.z -= other.z
            else:
                self.z = None
            return self
        else:
            raise TypeError('Point must be subtracted with Point')

    def __div__(self, other):
        """
        Returns self / other.

        @param other: Self will be divided with other.
        @type other: Numeric

        @return: The division of self and other.
        @rtype: map.Map.Point
        """
        if isinstance(other, (int, float)):
            my_x = self.x / float(other)
            my_y = self.y / float(other)
            if self.z is not None:
                my_z = self.z / float(other)
            else:
                my_z = None
            return Point(x=my_x, y=my_y, z=my_z)
        else:
            raise TypeError('Point must be divided with a number')

    def __idiv__(self, other):
        """
        Returns self /= other.

        @param other: Self will be divided with other.
        @type other: Numeric

        @return: The division of self and other.
        @rtype: map.Map.Point
        """
        if isinstance(other, (int, float)):
            if isinstance(other, (int, float)):
                self.x /= float(other)
                self.y /= float(other)
            if self.z is not None:
                self.z /= float(other)
            else:
                self.z = None
            return self
        else:
            raise TypeError('Point must be divided with a number')

