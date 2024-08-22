'''
Created on 23.05.2024
@author: poudel-b
'''

# -*- coding: utf-8 -*-

from math import exp, log
from . import iso_storages


class flux_connection(object):
    """
    Description
    ===========
    Base class for all flux connections.

    The connections hold the processes for the calculation of fluxes [kg] between storages and model boundaries based on a given or formulated flux rate [m/d].

    A positive flux rate [m/d] is always considered to be from left node to right node a  negative from right to left.

    """

    # Class attribute
    Dict_of_all_flux_connections = {}
    current_ID = 1

    def __init__(self,
                 left_node,
                 right_node
                 ):
        """
        Constructor of flux_connection

        @param left_node: The left node of this flux connection
        @type left_node: LSI_storages.iso_storage

        @param right_node: The right node of this flux connection
        @type right_node: LSI_storages.iso_storage

        @param flux_crosssectional_area: Area shared by the left and the right node of the flux connection in m2
        @type flux_crosssectional_area: float
        """

        try:
            assert isinstance(left_node, iso_storages.flux_node) and isinstance(right_node, iso_storages.flux_node), \
                "Given nodes must be an instance of LSI_storages.flux_node"
            # assert (flux_crosssectional_area > 0.0, "Given flux_crosssectional_area must be in m2 and >0.0")

            self.__ID = self.__class__.current_ID
            self.__class__.current_ID += 1  # increment the current ID for the next flux node to be created
            # add Instance of flux node to Dict_of_all_flux_nodes
            self.__class__.Dict_of_all_flux_connections[self.__ID] = self
            self.__left_node = left_node
            self.__right_node = right_node
            # self.__flux_crosssectional_area = flux_crosssectional_area

            # register connection
            left_node.RegisterConnection(newConnection=self)
            right_node.RegisterConnection(newConnection=self)

        except AssertionError as err:
            print(err)
            raise NotImplementedError

    def get_left_node(self):
        return self.__left_node

    left_node = property(get_left_node, None, None, "Returns the left node of the flux connection")

    def get_right_node(self):
        return self.__right_node

    right_node = property(get_right_node, None, None, "Returns the right node of the flux connection")

    def calc_flux(self, Isotopologue):
        """
        Positive flux of water [m/day] means left to right
        """
        raise NotImplementedError

    def calc_flux_i(self, Isotopologue):
        """
        Positive flux of isotopes [m/day] means left to right
        flux_i = flux * isotope_concentration
        """
        raise NotImplementedError

    def calc_flux_liquid(self, Isotopologue):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux
        """
        raise NotImplementedError

    def calc_q_i(self, Isotopologue, delta_time):
        """
        Positive q_i [kg Isotope] means flux left to right

        @param delta_time: Time step
        @type delta_time: datetime.timedelta
        """
        totalseconds_per_day = 86400  # 24* 60 * 60

        return  # self.calc_flux_i(Isotopologue) * self.flux_crosssectional_area *.total_seconds() / totalseconds_per_day

    def is_node(self, node):
        """
        Returns True if the given node is a node of this flux connection and False if not.

        @param node: Node to be tested
        @type node: SLI_storages.flux_node
        """
        if node == self.__left_node or node == self.__right_node:
            return True
        else:
            return False


class boundary_connection(flux_connection):
    """
    Description
    ===========
    Base class for all flux connections.

    The connections hold the processes for the calculation of fluxes [m/day] between water storages and model boundaries.

    A positive flux rate [m/d] is always considered to be from left node to right node a  negative from right to left.


    """

    def __init__(self,
                 left_node,
                 right_node
                 ):
        """
        Constructor of boundary_connection

        """
        flux_connection.__init__(self, left_node, right_node)


class liquid_diffusion_base_class(object):
    """
    Description
    ===========
    Functions to calculate diffusive fluxes in liquids
    """

    def __init__(self, q_l=None, hydrodynamic_dispersivity=0.0):
        self.q_l = q_l
        self.hydrodynamic_dispersivity = hydrodynamic_dispersivity

    def dli_eff_sli(self, dli_eff, theta_sat, theta_r, liq_sat):

        """
        Calculates the total (effective) liquid diffusivity according to SLI for given saturation i in the soil liquide phase (m**2/s)

        Description
        ===========

        @param theta_sat: Saturated liquid phase (m**3/m**3 )
        @type theta_sat: Float

        @param theta_r: Residual liquid phase (m**3/m**3 )
        @type theta_r: Float

        @param liq_sat: Liquid Saturation (m/m)
        @type liq_sat: Float

        """
        dli_sli = dli_eff * (min(liq_sat, 1) * (theta_sat - theta_r) + theta_r)

        return dli_sli

    def dl_i_eff(self, dl_i, q_l=None, hydrodynamic_dispersivity=0.0):
        """
        Calculates the total (effective) liquid diffusivity for isotope i in the soil liquide phase (m**2/s)

        Description
        ===========

        Under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones and are frequently set to zero (Auriault and Adler, 1995)

        @param dl_i: liquid diffusivity for isotope i in the soil liquide phase (m**2/day)
        @type dl_i: Float

        @param q_l: liqiud water flux (m/s)
        @type q_l: Float

        @param hydrodynamic_dispersivity: per default set to 0 (m) because under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones
        @type hydrodynamic_dispersivity: Float
        """
        try:
            if q_l is None or hydrodynamic_dispersivity == 0.0:
                dl_i_eff = dl_i
            else:  # q_l is not None and hydrodynamic_dispersivity != 0.0:
                dl_i_eff = dl_i + hydrodynamic_dispersivity * abs(q_l)

            return dl_i_eff

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def dl_i(self, T, Isotopologue, theta, tortuosity, formulation="Cuntz", **kwargs):
        """
        Calculates the liquid diffusivity for isotope i in the soil liquide phase (m**2/day)

        Description
        ===========

        @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
        @type T: Float

        @param Isotopologue: Identifier for the Isotopologue to be calculated (currently supported '2H' and '18O')
        @type Isotopologue: Sting

        @param theta: liquid phase (m**3/m**3 )
        @type theta: Float

        @param tortuosity: tortuosity (m/m)
        @type tortuosity: Float

        @param ignore_dl_i: If set to "True" liquid diffusivity for isotope will be set to 0 = no liquid diffusivity
        @type ignore_dl_i: Boolean

        """
        ignore_dl_i = kwargs.get('ignoredli', False)
        try:
            if not ignore_dl_i:
                # SLI: cable_sli_solve.f90::L2381-2390
                dl_i = self.dl_i_self_diffusivity(T, Isotopologue, formulation=formulation) * theta * tortuosity
            else:
                dl_i = 0.0

            return dl_i

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def dl_i_self_diffusivity(self, T, Isotopologue, formulation="Cuntz"):
        """
        Calculates the liquide self diffusivity of an given isotope in pure water (m**2/day) either after Melayah et al.  (1996) or Cuntz et al. (2007)

        Description
        ===========

        @param T: temperature T (K) of the current layer (283.15 K = 10 Celsius)
        @type T: Float

        @param Isotopologue: Identifier for the Isotopologue to be calculated (currently supported '2H' and '18O')
        @type Isotopologue: Sting

        @param formulation: Identifier for the formulation to be used to calculated dl_i_self_diffusivity (currently supported 'Melayah'and 'Cuntz')
        @type formulation: Sting

        """
        try:

            if formulation == "Melayah":
                if "18O" == Isotopologue:
                    a_i = 0.9669
                elif "2H" == Isotopologue:
                    a_i = 0.9833
                else:
                    raise NotImplementedError

                dl_i = a_i * 1.10 ** -9 * exp(-(535400 / T ** 2) + (1393.3 / T) + 2.1876)

            # SLI: cable_sli_solve.f90::L2381-2390
            elif formulation == "Cuntz":
                if "18O" == Isotopologue:
                    a_i = 1 / 1.026
                elif "2H" == Isotopologue:
                    a_i = 1 / 1.013
                else:
                    raise NotImplementedError

                dl_i = a_i * 1e-7 * exp(-577 / (T - 145))

            else:
                raise NotImplementedError

            return dl_i  # * 86400  # convert m2 per seconds to m2 per day

        except ValueError as err:
            print(err)
            raise NotImplementedError


class liquid_diffusion(flux_connection, liquid_diffusion_base_class):
    """
    Description
    ===========
    Diffusive flux of isotopes [kg] between to nodes.

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 left_node,
                 right_node,
                 hydrodynamic_dispersivity=0.0,
                 q_l=None):
        """
        Constructor of boundary_connection

        @param hydrodynamic_dispersivity: per default set to 0 (m) because under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones
        @type hydrodynamic_dispersivity: Float

        @param q_l: liqiud water flux (m/s) (only necessary if the hydrodynamic_dispersivity shall be calculated)
        @type q_l: Float
        """

        flux_connection.__init__(self, left_node, right_node)
        liquid_diffusion_base_class.__init__(self, q_l, hydrodynamic_dispersivity)

    def flux_to_left(self, Isotopologue, dl_i_formulation="Cuntz", **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, dl_i_formulation=dl_i_formulation, **kwargs)

    def flux_to_right(self, Isotopologue, dl_i_formulation="Cuntz", **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, dl_i_formulation=dl_i_formulation, **kwargs)

    def calc_flux(self, Isotopologue, dl_i_formulation="Cuntz", **kwargs):
        """
        Positive flux of the given Isotopologue [m/day] left to right
        """
        # Under evaporation, convective and hydrodynamic dispersion processes are negligible as
        # compared to the diffusion ones and are frequently set to zero (Auriault and Adler, 1995)

        dl_i_left_node = self.dl_i(T=self.left_node.T,
                                   Isotopologue=Isotopologue,
                                   theta=1,  # self.left_node.theta, ## using saturation instead of theta
                                   tortuosity=self.left_node.tortuosity,
                                   formulation=dl_i_formulation,
                                   **kwargs)

        dl_i_right_node = self.dl_i(T=self.right_node.T,
                                    Isotopologue=Isotopologue,
                                    theta=1,  # self.right_node.theta, Saturation instead of theta
                                    tortuosity=self.right_node.tortuosity,
                                    formulation=dl_i_formulation,
                                    **kwargs)

        dl_i_eff_left_node = self.dl_i_eff(dl_i=dl_i_left_node, q_l=self.q_l,
                                           hydrodynamic_dispersivity=self.hydrodynamic_dispersivity)
        dl_i_eff_right_node = self.dl_i_eff(dl_i=dl_i_right_node, q_l=self.q_l,
                                            hydrodynamic_dispersivity=self.hydrodynamic_dispersivity)

        liq_sat_l = self.left_node.get_liq_saturation()
        liq_sat_R = self.right_node.get_liq_saturation()

        dli_eff_sli_leftnode = self.dli_eff_sli(dli_eff=dl_i_eff_left_node,
                                                theta_sat=self.left_node.theta_sat,
                                                theta_r=self.left_node.theta_0,
                                                liq_sat=liq_sat_l)
        dli_eff_sli_rightnode = self.dli_eff_sli(dli_eff=dl_i_eff_right_node,
                                                 theta_sat=self.right_node.theta_sat,
                                                 theta_r=self.right_node.theta_0,
                                                 liq_sat=liq_sat_R)

        thickness_left_node = self.left_node.thickness  # thickness in m
        thickness_right_node = self.right_node.thickness  # thickness in m

        delta_z = self.left_node.center.distance(
            to_Point=self.right_node.center)  # distance between left and right node in cm
        # delta_z = self.right_node.center - self.left_node.center

        dl_mean = ((dli_eff_sli_leftnode * thickness_left_node + dli_eff_sli_rightnode * thickness_right_node) / (
                thickness_right_node + thickness_left_node)) / delta_z

        return dl_mean

    def calc_flux_liquid(self, Isotopologue, dl_i_formulation="Cuntz", **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, dl_i_formulation=dl_i_formulation, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):
        """
        Positive flux of the given Isotopologue [m/day] left to right

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_i = self.left_node.get_conc_iso_liquid(Isotopologue) * flux
        else:
            flux_i = self.right_node.get_conc_iso_liquid(Isotopologue) * flux

        """
        flux = self.calc_flux(Isotopologue=Isotopologue, **kwargs)
        flux_i = flux * (self.left_node.get_conc_iso_liquid(Isotopologue)
                         - self.right_node.get_conc_iso_liquid(Isotopologue))

        return flux_i


class vapor_diffusion_base_class(object):
    """
    Description
    ===========
    Functions to caclculate diffusive fluxes
    """

    def __init__(self):
        pass

    def dv_i(self, dv, Isotopologue, **kwargs):
        """
        Calculates the vapour diffusivity of a Isotopologue in free air (m**2/s)

        Description
        ===========

        @param dv: vapour diffusivity of water in free air (m**2/s)
        @type dv: Float

        @param Isotopologue: Identifier for the Isotopologue to be calculated (currently supported '2H' and '18O')
        @type Isotopologue: Sting

        @param ignore_dv_i:  If set to "True" vapour diffusivity of a Isotopologue in free air will be set to dv_i = dv --> the same as normal water
        @type ignore_dv_i: Boolean

        SLI: sli_init.f90 Lines: 2366 - L2367 <--checked and own implementation is okay
        """
        ignore_dv_i = kwargs.get('ignoredvi', False)
        try:
            # SLI: cable_sli_solve.f90::L2366 - L2367
            if "2H" == Isotopologue:
                b_i = 1.0251  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
            elif "18O" == Isotopologue:
                b_i = 1.0285  # HDO diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
            else:
                raise NotImplementedError

            if ignore_dv_i:
                return dv
            else:
                # SLI: cable_sli_solve.f90::L2373 - L2374
                dv_i = dv / b_i  # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff

                return dv_i

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def dv_soil_air(self, T, theta, theta_sat, tortuosity):
        """
        Calculates the vapour diffusivity of water in soil air space (m**2/s)

        Description
        ===========

        @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
        @type T: Float

        @param dv_free_air: diffusivity of water in free air (m**2/s)
        @type dv_free_air: Float

        @param theta: liquid phase (m**3/m**3 )
        @type theta: Float

        @param theta_sat: theta at saturation equal to soil porosity (m**3/m**3 )
        @type theta_sat: Float

        @param tortuosity: tortuosity (m/m)
        @type tortuosity: Float
        """
        try:
            # SLI:: cable_sli_utils.f90 L::1494: var%Dv    = Dva*parin%tortuosity*(parin%the-theta)  * ((Tsoil+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
            # TODO: Change to actual atmospheric pressure
            dv_air = self.dv_free_air(T=T, Pa=1e5)  ## for 1e5 / Pa = 1 in SLI_utils #
            dv_soil_air = dv_air * tortuosity * (theta_sat - theta)

            return dv_soil_air

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def dv_free_air(self, T=283.15, Pa=10 ** 5):
        """
        Calculates the vapour diffusivity of water in free air (m**2/s)

        Description
        ===========

        @param T: temperature T (K) (283.15 K = 10 Celsius) <-- soil or atmosphere????
        @type T: Float

        @param Pa: air pressure in the atmosphere (1bar = 1 atm = 10**5 Pa) [Pa]
        @type Pa: Float

        TODO: SLI uses 1 for Pa and not 100.000 why????

        """
        try:
            dva = 2.1699999706470408E-005  # vapour diffusivity of water in air at 0 degC [m2/s]
            dv = dva * 1e5 / Pa * (T / 273.16000366210938) ** 1.88
            # SLI: Dva       = 2.17e-5   ! vapour diffusivity of water in air at 0 degC [m2/s]
            # SLI: Dvs  = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffuxivity of water in air (m2s-1)

            return dv

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def nD(theta, theta_0=0.05):
        """
        Calculates the exponent that describes the isotopic fractioning due to diffusion (-) Melayah et al. (1996)


        Description
        ===========


        @param theta: liquid phase (m**3/m**3 )
        @type theta: Float

        @param theta_0: defined residual water content (water remaining at high suctions) (m**3/m**3 ) --> should be 0
        @type theta_0: Float


        """
        try:
            nD = 0.67 + 0.33 * exp(1 - theta / theta_0)

            return nD

        except ZeroDivisionError as err:
            return 0.67


class vapor_diffusion(flux_connection, vapor_diffusion_base_class):
    """
    Description
    ===========
    Diffusive flux of isotopes [kg] between to nodes.

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 left_node,
                 right_node
                 ):
        """
        Constructor of boundary_connection

        """
        # assert isinstance(left_node and right_node, iso_storages.iso_soil_layer)
        flux_connection.__init__(self, left_node, right_node)
        vapor_diffusion_base_class.__init__(self)

    def flux_to_left(self, Isotopologue, **kwargs):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_liquid = self.left_node.alpha_i(Isotopologue) * flux
        else:
            flux_liquid = self.right_node.alpha_i(Isotopologue) * flux
        """
        # TODO: Check beta in sli_solve:L2515
        flux = self.calc_flux(Isotopologue, **kwargs)
        beta = self.left_node.beta(Isotopologue=Isotopologue, **kwargs)

        flux_v = flux * beta

        return flux_v

    def flux_to_right(self, Isotopologue, **kwargs):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_liquid = self.left_node.alpha_i(Isotopologue) * flux
        else:
            flux_liquid = self.right_node.alpha_i(Isotopologue) * flux
        """
        # TODO: Check beta in sli_solve:L2515
        flux = self.calc_flux(Isotopologue, **kwargs)
        beta = self.right_node.beta(Isotopologue=Isotopologue, **kwargs)

        flux_v = flux * beta

        return flux_v

    def calc_flux(self, Isotopologue, **kwargs):
        """
        Positive flux of the given Isotopologue [m/day] left to right
        """

        left_node_dv_soil_air = self.dv_soil_air(T=self.left_node.T,
                                                 theta=self.left_node.theta,
                                                 theta_sat=self.left_node.theta_sat,
                                                 tortuosity=self.left_node.tortuosity)

        right_node_dv_soil_air = self.dv_soil_air(T=self.right_node.T,
                                                  theta=self.right_node.theta,
                                                  theta_sat=self.right_node.theta_sat,
                                                  tortuosity=self.right_node.tortuosity)

        left_node_dv_i = self.dv_i(dv=left_node_dv_soil_air, Isotopologue=Isotopologue, **kwargs)
        right_node_dv_i = self.dv_i(dv=right_node_dv_soil_air, Isotopologue=Isotopologue, **kwargs)

        leftnode_dvi_sli, rightnode_dvi_sli = left_node_dv_i * self.left_node.cv, right_node_dv_i * self.right_node.cv

        delta_z = self.left_node.center.distance(to_Point=self.right_node.center)
        # delta_z = self.right_node.center - self.left_node.center  # distance between left and right node in m

        dv_i_mean = (leftnode_dvi_sli * self.left_node.thickness + rightnode_dvi_sli * self.right_node.thickness) / (
                self.left_node.thickness + self.right_node.thickness) / delta_z

        return dv_i_mean

    def calc_flux_liquid(self, Isotopologue, **kwargs):

        T_mean = (self.left_node.T + self.right_node.T) / 2  # average tempereture between nodes
        T_mean_0 = (self.left_node.T0 + self.right_node.T0) / 2  # average tempereture between nodes previous time

        dT_mean = T_mean - T_mean_0

        alpha_i_mean = self.left_node.alpha_i(Isotopologue=Isotopologue, T=T_mean, **kwargs)
        d_beta_mean = self.left_node.d_beta(Isotopologue=Isotopologue, T=T_mean, dT=dT_mean, **kwargs)

        beta_mean = alpha_i_mean + d_beta_mean

        flux = self.calc_flux(Isotopologue, **kwargs)
        flux_v = flux * beta_mean

        return flux_v

    def calc_flux_liquid_2(self, Isotopologue, **kwargs):
        beta_mean = (self.left_node.beta(Isotopologue=Isotopologue, **kwargs)
                     + self.right_node.beta(Isotopologue=Isotopologue, **kwargs)) / 2

        flux = self.calc_flux(Isotopologue, **kwargs)
        flux_v = flux * beta_mean

        return flux_v

    def calc_flux_i(self, Isotopologue, **kwargs):
        """
        Positive flux of the given Isotopologue [m/day] left to right

        if flux >= 0.0:
            flux_i = self.left_node.get_conc_iso_vapor(Isotopologue) * flux
        else:
            flux_i = self.right_node.get_conc_iso_vapor(Isotopologue) * flux

        """

        T_mean = (self.left_node.T + self.right_node.T) / 2  # average tempereture between nodes
        T_mean_0 = (self.left_node.T0 + self.right_node.T0) / 2  # average tempereture between nodes previous time

        dT_mean = T_mean - T_mean_0

        alpha_i_mean = self.left_node.alpha_i(Isotopologue=Isotopologue, T=T_mean, **kwargs)
        d_beta_mean = self.left_node.d_beta(Isotopologue=Isotopologue, T=T_mean, dT=dT_mean, **kwargs)

        beta_mean = alpha_i_mean + d_beta_mean

        flux = self.calc_flux(Isotopologue, **kwargs)

        # beta_left = self.left_node.beta(Isotopologue=Isotopologue, **kwargs)
        # beta_right = self.right_node.beta(Isotopologue=Isotopologue, **kwargs)

        flux_i = flux * (self.left_node.get_conc_iso_liquid(Isotopologue)
                         - self.right_node.get_conc_iso_liquid(Isotopologue)) * beta_mean

        return flux_i


class liquid_advection(flux_connection):
    """
    Description
    ===========
    Advective flux of isotopes [kg] between to nodes based on a given liquid flux [m/day]

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 left_node,
                 right_node,
                 q_l=0.0  # Liquide flux in m/day positive form left to right : negative from right to left
                 ):
        """
        Constructor of boundary_connection

        @param q_l: Liquide flux in m3/(m2*day). Positive form left to right : negative from right to left
        @type q_l: float
        """
        flux_connection.__init__(self, left_node, right_node)
        self.q_l = q_l

    def flux_to_left(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def flux_to_right(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def get_flux(self):
        """
        Positive flux of the given Isotopologue [m/day] left to right
        """
        flux = self.q_l

        return flux

    def calc_flux(self, Isotopologue, **kwargs):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux
        """
        return self.get_flux() * 0.5  # sli_solve:L2432, dcqldca & dcqldcb = 0.5

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):
        """
        Positive flux of the given Isotopologue [m/day] left to right

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_i = left_node.get_conc_iso_liquid(Isotopologue) * flux
        else:
            flux_i = right_node.get_conc_iso_liquid(Isotopologue) * flux

        """

        flux = self.get_flux()
        cql = 0.5 * (self.left_node.get_conc_iso_liquid(Isotopologue)
                     + self.right_node.get_conc_iso_liquid(Isotopologue))  # sli_solve:L2431
        flux_i = flux * cql
        return flux_i


class vapor_advection(flux_connection, vapor_diffusion_base_class):
    """
    Description
    ===========
    Advective flux of isotopes [kg] between to nodes based on a given vapor flux [m/day]

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 left_node,
                 right_node,
                 q_v=0.0,
                 # vapor flux in m/day positive form left to right : negative from right to left. If None q_v  will be calculated based on equation A.12 Haverd&Cuntz (2005) knowing futurue Temp and matrix potential
                 ):
        """
        Constructor of boundary_connection

        @param q_v: vapor flux in (m/day). Positive form left to right : negative from right to left
        @type q_v: float
        """
        vapor_diffusion_base_class.__init__(self)
        flux_connection.__init__(self, left_node, right_node)
        self.q_v = q_v

    def flux_to_left(self, Isotopologue, **kwargs):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_liquid = left_node.alpha_i(Isotopologue) * flux
        else:
            flux_liquid = right_node.alpha_i(Isotopologue) * flux

        """
        flux = self.get_flux()
        beta_leftnode = self.left_node.beta(Isotopologue, **kwargs)

        # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff
        betaq = self.dv_i(dv=1, Isotopologue=Isotopologue, **kwargs)

        flux_liquid = flux * betaq * beta_leftnode * 0.5  # sli_solve:L2439

        return flux_liquid

    def flux_to_right(self, Isotopologue, **kwargs):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_liquid = left_node.alpha_i(Isotopologue) * flux
        else:
            flux_liquid = right_node.alpha_i(Isotopologue) * flux

        """
        flux = self.get_flux()
        beta_rightnode = self.right_node.beta(Isotopologue=Isotopologue, **kwargs)
        # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff
        betaq = self.dv_i(dv=1, Isotopologue=Isotopologue, **kwargs)

        flux_liquid = flux * betaq * beta_rightnode * 0.5  # sli_solve:L2439

        return flux_liquid

    def get_flux(self):
        """
        Positive flux of the given Isotopologue [m/day] left to right
        """

        flux = self.q_v

        return flux

    def calc_flux_liquid(self, Isotopologue, **kwargs):

        T_mean = (self.left_node.T + self.right_node.T) / 2  # average tempereture between nodes
        T_mean_0 = (self.left_node.T0 + self.right_node.T0) / 2  # average tempereture between nodes previous time

        dT_mean = T_mean - T_mean_0

        alpha_i_mean = self.left_node.alpha_i(Isotopologue=Isotopologue, T=T_mean, **kwargs)
        d_beta_mean = self.left_node.d_beta(Isotopologue=Isotopologue, T=T_mean, dT=dT_mean, **kwargs)

        beta_mean = alpha_i_mean + d_beta_mean

        flux = self.get_flux()
        # beta_rightnode = self.right_node.beta(Isotopologue=Isotopologue, **kwargs)
        # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff
        betaq = self.dv_i(dv=1, Isotopologue=Isotopologue, **kwargs)

        flux_v = flux * betaq * beta_mean * 0.5  # sli_solve:L2439

        return flux_v

    def calc_flux_liquid_2(self, Isotopologue, **kwargs):
        flux = self.get_flux()

        beta_mean = (self.left_node.beta(Isotopologue=Isotopologue, **kwargs)
                     + self.right_node.beta(Isotopologue=Isotopologue, **kwargs)) / 2
        # beta_rightnode = self.right_node.beta(Isotopologue=Isotopologue, **kwargs)
        # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff
        betaq = self.dv_i(dv=1, Isotopologue=Isotopologue, **kwargs)

        flux_v = flux * betaq * beta_mean * 0.5  # sli_solve:L2439

        return flux_v

    def calc_flux_i(self, Isotopologue, **kwargs):
        """
        Positive flux of the given Isotopologue [m/day] left to right

        flux = self.calc_flux(Isotopologue)
        if flux >= 0.0:
            flux_i = left_node.get_conc_iso_vapor(Isotopologue) * flux
        else:
            flux_i = right_node.get_conc_iso_vapor(Isotopologue) * flux

        """

        T_mean = (self.left_node.T + self.right_node.T) / 2  # average tempereture between nodes
        T_mean_0 = (self.left_node.T0 + self.right_node.T0) / 2  # average tempereture between nodes previous time

        dT_mean = T_mean - T_mean_0

        alpha_i_mean = self.left_node.alpha_i(Isotopologue=Isotopologue, T=T_mean, **kwargs)
        d_beta_mean = self.left_node.d_beta(Isotopologue=Isotopologue, T=T_mean, dT=dT_mean, **kwargs)

        beta_mean = alpha_i_mean + d_beta_mean

        flux = self.get_flux()

        # beta_leftnode = self.left_node.beta(Isotopologue=Isotopologue, **kwargs)
        # beta_rightnode = self.right_node.beta(Isotopologue=Isotopologue, **kwargs)

        cqv = 0.5 * (self.left_node.get_conc_iso_liquid(Isotopologue)
                     + self.right_node.get_conc_iso_liquid(Isotopologue)) * beta_mean

        # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff
        betaq = self.dv_i(dv=1, Isotopologue=Isotopologue, **kwargs)

        flux_i = flux * betaq * cqv  # sli_solve:L2438

        return flux_i


class evaporation(boundary_connection, vapor_diffusion_base_class, liquid_diffusion_base_class):
    """
    Description
    ===========
    Flux of isotopes [kg] between two nodes (atmosphere and upper soil layer) due to evaporation.

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 atmosphere,
                 top_layer,
                 q_Evap=0.0,
                 # evaporative flux from soil or litter surface to atmosphere (kg/(m2*day)) -> should be a negative number flowing from UpperSoilLayer to Atmosphere
                 T_surface=273.0,  # soil surface (interface atmosphere / soil) temparature in Kelvin
                 q_l=0.0,
                 # liquid flux between the topsoil layer and the atmosphere (m/day) - including precipitation ->in case of percepitation it should be a positive number flowing from Atmosphere to UpperSoilLayer
                 q_v=0.0,  # vapour flux into the soil
                 hydrodynamic_dispersivity=0.0):  # hydrodynamic_dispersivity for liquide diffusion will be ignored when set to 0.0

        """
        Constructor of evaporation

        @param q_Evaporation: Flux of evaporation in kg/(m2*s)
        @type q_Evaporation: float

        @param T_surface: surface temperature at the interface atmosphere / soil) K
        @type q_Evaporation: float

        """
        if isinstance(atmosphere, iso_storages.iso_atmosphere) and isinstance(top_layer, iso_storages.iso_soil_layer):

            boundary_connection.__init__(self, left_node=atmosphere, right_node=top_layer)

            self.atmosphere = atmosphere
            self.top_layer = top_layer
        else:
            raise NotImplementedError

        vapor_diffusion_base_class.__init__(self)
        liquid_diffusion_base_class.__init__(self, q_l, hydrodynamic_dispersivity)
        self.q_evap = q_Evap
        self.T_surface = T_surface
        self.ql = q_l
        self.qv = q_v

    def get_flux(self):
        return self.ql

    def calc_flux(self, Isotopologue, **kwargs):

        delta_c_evapout = self.delta_evapout_c_iso(Isotopologue, **kwargs)
        flux = self.q_evaporation_out * delta_c_evapout
        return flux

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, formulation_alpha_i_k="MathieuBariac", formulation_dl_i="Cuntz", **kwargs):

        ciso_surface = self.c_iso_liq_surface(Isotopologue=Isotopologue, formulation_dl_i=formulation_dl_i, **kwargs)
        nk = self.nk_sli_solve(thetasat_surface=self.top_layer.theta_sat, Sl=self.top_layer.Sl)

        c_evap_out = self.get_conc_evapout(Isotopologue=Isotopologue, c_iso_surface=ciso_surface, nk=nk, **kwargs)
        c_evap_in = self.get_conc_evapin(Isotopologue=Isotopologue, nk=nk, **kwargs)

        evap_flux_in = self.q_evaporation_in * c_evap_in
        evap_flux_out = self.q_evaporation_out * c_evap_out

        evap_flux = evap_flux_out - evap_flux_in

        return evap_flux

    def get_conc_evapin(self, Isotopologue, nk, **kwargs):

        """"
        Isotopic concentration in the evaporation flux entering soil
        """
        dv_surface = self.dv_free_air(T=self.T_surface, Pa=self.atmosphere.Pa)
        div_surface = self.dv_i(dv=dv_surface, Isotopologue=Isotopologue, **kwargs)

        alpha_i_k = 1 / self.alpha_i_k(dv=dv_surface, dv_i=div_surface, nK_MathieuBariac=nk,
                                       formulation="MathieuBariac", **kwargs)

        c_in = self.atmosphere.civ_a(Isotopologue=Isotopologue) / self.atmosphere.cv_a * alpha_i_k

        return c_in

    def get_conc_evapout(self, Isotopologue, c_iso_surface, nk, **kwargs):

        """"
        Isotopic concentration in the evaporation flux leaving soil
        """
        dv_surface = self.dv_free_air(T=self.T_surface, Pa=self.atmosphere.Pa)
        div_surface = self.dv_i(dv=dv_surface, Isotopologue=Isotopologue, **kwargs)

        alpha_i_k = 1 / self.alpha_i_k(dv=dv_surface, dv_i=div_surface, nK_MathieuBariac=nk,
                                       formulation="MathieuBariac", **kwargs)

        alphai_surface = self.top_layer.alpha_i(Isotopologue=Isotopologue, T=self.T_surface, **kwargs)
        c_out = alphai_surface * alpha_i_k * c_iso_surface

        return c_out

    def get_evap_influx(self):

        """Evaporation flux entering soil"""

        # sli_main.f90: L359 / rbh in sli_solve.f90:L2391

        ram = self.ram(self.atmosphere.wind_speed, self.atmosphere.hc, self.atmosphere.d0, self.atmosphere.z0m)
        rbh = self.rbh(self.atmosphere.wind_speed, self.atmosphere.extku, self.atmosphere.LAI)

        q_in = self.atmosphere.cv_a / (ram + rbh)

        return q_in

    q_evaporation_in = property(get_evap_influx, None, None, "Evaporation flux entering soil")

    def get_evap_outflux(self):

        """Evaporation flux leaving soil"""

        # sli_main.f90: L359 / rbh in sli_solve.f90:L2392

        ram = self.ram(self.atmosphere.wind_speed, self.atmosphere.hc, self.atmosphere.d0, self.atmosphere.z0m)
        rbh = self.rbh(self.atmosphere.wind_speed, self.atmosphere.extku, self.atmosphere.LAI)

        q_out = self.cv_surface / (ram + rbh)

        return q_out

    q_evaporation_out = property(get_evap_outflux, None, None, "Evaporation flux exiting soil")

    def get_cv_surface(self):

        """" Returns concentration of water vapour at soil/air interface (m3 (H2O liq)/ m3 (air))"""

        # sli_main.f90: L359 / rbh in sli_solve.f90:L2323

        ram = self.ram(self.atmosphere.wind_speed, self.atmosphere.hc, self.atmosphere.d0, self.atmosphere.z0m)
        rbh = self.rbh(self.atmosphere.wind_speed, self.atmosphere.extku, self.atmosphere.LAI)

        cv_surface = self.atmosphere.cv_a + self.q_evap * (ram + rbh)

        return cv_surface

    cv_surface = property(get_cv_surface, None, None, "concentration of water vapour at soil/air interface")

    def c_iso_liq_surface(self, Isotopologue, formulation_alpha_i_k="MathieuBariac", formulation_dl_i="Cuntz",
                          BA83=False, **kwargs):

        """
        Calculates the concentration of minor isotopologue in liquid water at the surface (kg/m3_H2O)

        sli_solve.f90:L2323:2389
        """
        dv_surface = self.dv_free_air(T=self.T_surface, Pa=self.atmosphere.Pa)
        div_surface = self.dv_i(dv=dv_surface, Isotopologue=Isotopologue, **kwargs)

        if BA83:  # analytical solution for Barnes and Allison (1983)  saturated case
            nk = 1
        else:
            nk = self.nk_sli_solve(thetasat_surface=self.top_layer.theta_sat, Sl=self.top_layer.Sl)

        alpha_i_k = 1 / self.alpha_i_k(dv=dv_surface, dv_i=div_surface, nK_MathieuBariac=nk,
                                       formulation=formulation_alpha_i_k, **kwargs)

        dv_top_layer = self.dv_soil_air(self.top_layer.T, self.top_layer.theta,
                                        self.top_layer.theta_sat, self.top_layer.tortuosity)
        dvi_top_layer = self.dv_i(dv=dv_top_layer, Isotopologue=Isotopologue, **kwargs)

        dli = self.dl_i(self.top_layer.T, Isotopologue=Isotopologue, theta=1,
                        tortuosity=self.top_layer.tortuosity, **kwargs)
        dl_eff = self.dl_i_eff(dli, q_l=self.q_l)
        dl = self.dli_eff_sli(dl_eff, theta_sat=self.top_layer.theta_sat, theta_r=self.top_layer.theta_0,
                              liq_sat=self.top_layer.Sl)

        alphai_top_layer = self.top_layer.alpha_i(Isotopologue=Isotopologue, T=self.top_layer.T, **kwargs)
        alphai_surface = self.top_layer.alpha_i(Isotopologue=Isotopologue, T=self.T_surface, **kwargs)

        ram = self.ram(self.atmosphere.wind_speed, self.atmosphere.hc, self.atmosphere.d0, self.atmosphere.z0m)
        rbh = self.rbh(self.atmosphere.wind_speed, self.atmosphere.extku, self.atmosphere.LAI)

        cvs = self.cv_surface  # concentration of water vapour at soil/air interface (m3 (H2O liq)/ m3 (air))
        c_iso_atm = self.atmosphere.civ_a(
            Isotopologue=Isotopologue)  # atmospheric isotopic concentration w.r.t relative humidity / vapour concentration

        # upper part okay
        if self.ql > 0.0:  # liquide water flux (m/day) between the atmosphere and the top soil layer (negative value leaving the cell positive entering the cell)?
            w1 = 0  # in case that water enters the soil
        else:
            w1 = 1  # in case that water leaves the cell

        if dv_top_layer != 0.0:
            cv_toplayer = - self.qv * (0.5 * self.top_layer.thickness) / dv_top_layer + cvs
        else:
            cv_toplayer = cvs

        num = dvi_top_layer / (0.5 * self.top_layer.thickness) * cv_toplayer \
              * alphai_top_layer * self.top_layer.get_conc_iso_liquid(Isotopologue) \
              - self.ql * self.top_layer.get_conc_iso_liquid(Isotopologue) * w1 \
              + c_iso_atm * alpha_i_k / (ram + rbh) \
              + dl * self.top_layer.get_conc_iso_liquid(Isotopologue) / (0.5 * self.top_layer.thickness)

        den = alpha_i_k * cvs * alphai_surface / (ram + rbh) + self.ql * (1 - w1) \
              + alphai_top_layer * dvi_top_layer * cvs / (0.5 * self.top_layer.thickness) \
              + dl / (0.5 * self.top_layer.thickness)

        ci_surface = num / den

        return ci_surface

    def delta_evapout_c_iso(self, Isotopologue, **kwargs):

        # sli_solve.f90:
        dv_surface = self.dv_free_air(T=self.T_surface, Pa=self.atmosphere.Pa)
        div_surface = self.dv_i(dv=dv_surface, Isotopologue=Isotopologue, **kwargs)

        nk = self.nk_sli_solve(thetasat_surface=self.top_layer.theta_sat, Sl=self.top_layer.Sl)
        alpha_i_k = 1 / self.alpha_i_k(dv=dv_surface, dv_i=div_surface, nK_MathieuBariac=nk,
                                       formulation="MathieuBariac", **kwargs)

        dv_top_layer = self.dv_soil_air(self.top_layer.T, self.top_layer.theta,
                                        self.top_layer.theta_sat, self.top_layer.tortuosity)
        dvi_top_layer = self.dv_i(dv=dv_top_layer, Isotopologue=Isotopologue, **kwargs)

        # theta=1 sli.solve:L2285  ??
        dli = self.dl_i(self.top_layer.T, Isotopologue=Isotopologue, theta=1,
                        tortuosity=self.top_layer.tortuosity, **kwargs)

        dl_eff = self.dl_i_eff(dli, q_l=self.q_l)
        dl = self.dli_eff_sli(dl_eff, theta_sat=self.top_layer.theta_sat, theta_r=self.top_layer.theta_0,
                              liq_sat=self.top_layer.Sl)

        alphai_top_layer = self.top_layer.alpha_i(Isotopologue=Isotopologue, T=self.top_layer.T,
                                                  **kwargs)
        alphai_surface = self.top_layer.alpha_i(Isotopologue=Isotopologue, T=self.T_surface,
                                                **kwargs)

        ram = self.ram(self.atmosphere.wind_speed, self.atmosphere.hc, self.atmosphere.d0, self.atmosphere.z0m)
        rbh = self.rbh(self.atmosphere.wind_speed, self.atmosphere.extku, self.atmosphere.LAI)

        # upper part okay
        if self.ql > 0.0:  # liquide water flux (m/day) between the atmosphere and the top soil layer (negative value leaving the cell positive entering the cell)?
            w1 = 0  # in case that water enters the soil
        else:
            w1 = 1  # in case that water leaves the cell

        cvs = self.cv_surface

        den = alpha_i_k * cvs * alphai_surface / (ram + rbh) + self.ql * (1 - w1) \
              + alphai_top_layer * dvi_top_layer * cvs / (0.5 * self.top_layer.thickness) \
              + dl / (0.5 * self.top_layer.thickness)

        delta_ciso = alpha_i_k * alphai_surface
        d_evapout_ciso = delta_ciso * (dvi_top_layer / (0.5 * self.top_layer.thickness)
                                       * self.top_layer.cv * alphai_top_layer
                                       - self.ql + dl / (0.5 * self.top_layer.thickness)) / den

        return d_evapout_ciso

    def alpha_i_k(self, dv=None, dv_i=None, nK_MathieuBariac=None, nK_Brutsaert=None, rbh=None, ram=None,
                  formulation="Melayah", **kwargs):
        """
        Calculates isotopic kinetic fractionation factor under non-saturated and non-isothermal conditions at the soil surface

        Description
        ===========

        Follogwing formulations are available:

        BarnesAllsion: Based on Barnes and Allsion (1983) --> alpha_i_k = dv /dv_i
        MathieuBariac: Based on Mathieu and Bariac (1996) --> alpha_i_k = (dv /dv_i)**nK_MathieuBariac
        Melayah: Based on Melayah et al. (1996) --> alpha_i_k = 1
        Brutsaert: Based on Brutsaert (1982) -->alpha_i_k =  1 + nK_Brutsaert*(dv/dv_i-1)*(r_am/r_a)  <-- with aerodynamic resistance (r_a) = turbulent resistance (r_aT) + molecular resistance (r_am); note nK_Brutsaert != nk MathieuBariac

        @param dv: vapour diffusivity of water in air (m**2/s)
        @type dv: Float

        @param dv_i: vapour diffusivity of isotope i in air (m**2/s)
        @type dv_i: Float

        @param nK_MathieuBariac: exponent that describes the isotopic kinetic fractioning at the soil surface used by the formulation of 'MathieuBariac'
        @type nK_MathieuBariac: float

        @param nK_Brutsaert: exponent that describes the isotopic kinetic fractioning at the soil surface used by the formulation of 'Brutsaert' depending on friction velcity, roughness length and kinematic velocity
        @type nK_Brutsaert: float
        nK_Brutsaert

        @param rbh: rbh =  ra + rs resistance to water vapour transfer beween surface and lowest atmospheric layer (s*m**-1)  used by the formulation of 'MathieuBariac'
        @type rbh: float

        @param ram: molecular resistance between level  hc the soil surface  z0m (s*m**-1)  used by the formulation of 'MathieuBariac'
        @type ram: float

        @param formulation: Identifier for the formulation to be used to calculated alpha_i_k (currently supported 'Melayah', 'BarnesAllsion', 'MathieuBariac' and 'Brutsaert')
        @type formulation: Sting

        @param ignore_alpha_i_k:  If set to "True" kinetic isotopic frctionation factor will be set to 1 = no kinetic fractionation
        @type ignore_alpha_i_k: Boolean
        """
        ignore_alpha_i_k = kwargs.get('ignorealphaik', False)
        try:
            if not ignore_alpha_i_k:
                # SLI: cable_sli_solve.f90::L2275 - L2305
                if formulation == "Melayah":
                    alpha_i_k = 1  # SLI: cable_sli_solve.f90::L2379
                elif formulation == "BarnesAllsion":
                    if dv is not None and dv_i is not None:
                        alpha_i_k = dv / dv_i
                    else:
                        print("The formualtion of alpha_i_k based on Barnes and Allsion (1983) requires dv  & dv_i")
                        raise NotImplementedError
                elif formulation == "MathieuBariac":
                    if nK_MathieuBariac is not None and dv is not None and dv_i is not None:
                        alpha_i_k = (dv / dv_i) ** nK_MathieuBariac  ##SLI: cable_sli_solve.f90::L2377
                    else:
                        print \
                            ("The formualtion of alpha_i_k based on Mathieu and Bariac "
                             "(1996) requires nK_MathieuBariac, dv  & dv_i")
                        raise NotImplementedError
                elif formulation == "Brutsaert":
                    if nK_Brutsaert is not None and dv != None and dv_i is not None \
                            and ram is not None and rbh is not None:
                        alpha_i_k = 1 + nK_Brutsaert * (dv / dv_i - 1) * (ram / rbh)  ##SLI: cable_sli_solve.f90::L2378
                    else:
                        print("The formualtion of alpha_i_k based on Mathieu "
                              "and Bariac (1996) requires nK, dv, dv_i, r_a & r_am")
                        raise NotImplementedError
                else:
                    raise NotImplementedError
            else:
                alpha_i_k = 1.0

            return alpha_i_k

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def nK_MathieuBariac(self, theta_surface, theta_surface_sat, theta_surface_0=0.05):
        """
        Calculates the exponent that describes the isotopic kinetic fractioning at the soil surface (nK = 1 --> dry soils ; nK = 0.5 --> saturated soils) based on Mathieu and Bariac (1996)

        Description
        ===========

        @param theta_surface: liquid phase at soil surface (m**3/m**3 )
        @type theta_surface: Float

        @param theta_surface_sat: water content at saturation (water remaining at high suctions) (m**3/m**3 )
        @type theta_surface_sat: Float

        @param theta_surface_0: defined residual water content (water remaining at high suctions) (m**3/m**3 ) --> should be 0
        @type theta_surface_0: Float

        SLI: cable_sli_solve.f90::L2375 <--checked and own implementation is the same
        """
        try:
            # SLI: cable_sli_solve.f90::L2375
            n_a = 0.5
            n_s = 1.0

            nK = ((theta_surface - theta_surface_0) * n_a + (theta_surface_sat - theta_surface) * n_s) / (
                    theta_surface_sat - theta_surface_0)

            return nK

        except ValueError as err:
            return 0.67

    def nk_sli_solve(self, thetasat_surface, Sl):
        ## Different formulation in SLI: Appendix : B.7

        """
        Calculates the exponent that describes the isotopic kinetic fractioning at the soil surface (nK = 1 --> dry soils ; nK = 0.5 --> saturated soils) based on Mathieu and Bariac (1996)

        Description
        ===========

        @param theta_surface: liquid phase at soil surface (m**3/m**3 )
        @type theta_surface: Float
        """

        nk = ((thetasat_surface * min(Sl, 1) - 0) * 0.5 + (thetasat_surface * (1 - min(Sl, 1)))) \
             / (thetasat_surface - 0)

        # nk = 1 if testcase == 7 or 8: sli_solve

        return nk

    def rbh(self, wind_speed, extku, LAI):
        """
        Returns the soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).

        Description
        ===========
        @param wind_speed: wind speed at top of canopy [m/s]
        @type wind_speed: float

        @param LAI: leaf area index (e.g. 2.0)
        @type LAI: float

        @param extku: extinction coeff't for windspeed (e.g. 0.5)
        @type extku: float

        @return:  Returns the soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).

        Control status: Checked on 15.05.2013 --> Results are the same as for SLI
        vmet%rs in sli_main.f90:L359 / rbh in sli_solve.f90:L2323
        """
        us = wind_speed * exp(-extku * LAI)  # us = windspeed below canopy [m/s]
        rbh = 1. / (0.004 + 0.012 * us)

        return rbh

    def ram(self, wind_speed, hc, d0, z0m):
        """
        Returns the aerodyamic resistance (from z0m to hc).

        Description
        ===========
        @param wind_speed: wind speed at top of canopy [m/s]
        @type wind_speed: float

        @param hc: canopy height [m] (e.g. 40)
        @type hc: float

        @param d0: displacement height (e.g. 0.7)
        @type d0: float

        @param z0m: roughness height for momentum (e.g. 0.1)
        @type z0m: float

        @return: Returns the aerodyamic resistance (from z0m to hc).

        Control status: Checked on 15.05.2013 --> Results are the same as for SLI
        vmet%ra in sli_main.f90:L357 / ram in sli_solve.f90:L2323
        """
        ram = log((hc - d0) / z0m) ** 2 / 0.41 ** 2 / wind_speed

        return ram

    def rbw(self, ram, rbh):
        """
        Returns the boundary layer resistance. Resistance to water vapour transfer beween surface and lowest atmospheric layer.

        Description
        ===========
        @param ram: Aerodyamic resistance (from z0m to hc)
        @type ram: float

        @param rs: Soil laminar boundary layer resistance  (Kustas & Norman AFM 94 (1999) 13-29).
        @type rs: float

        @return: Returns the boundary layer resistance. Resistance to water vapour transfer beween surface and lowest atmospheric layer.

        Control status: Checked on 15.05.2013 --> Results are the same as for SLI
        vmet%rbw in sli_main.f90:L361
        """
        rbh = ram + rbh  # resistance to water vapour transfer beween surface and lowest atmospheric layer

        return rbh


class transpiration(boundary_connection):
    """
    Description
    ===========
    Flux of isotopes [kg] due to transpiratation (root extraction) between a soil storage and the atmosphere based on a given liquid flux [m/day]

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 atmosphere,
                 soil_layer,
                 ql_transpiration=0.0,  # negative liquid flux in (kg/(m2*day)) from node_SoilLayer to node_Atmosphere.
                 ):
        """
        Constructor of boundary_connection

        @param ql_Transpiration: negative liquid flux in (m3/(m2*day) or m/day) from node_SoilLayer to node_Atmosphere.
        @type ql_Transpiration: float
        """

        if isinstance(atmosphere, iso_storages.iso_atmosphere) and isinstance(soil_layer, iso_storages.iso_soil_layer):

            boundary_connection.__init__(self, left_node=atmosphere, right_node=soil_layer)
            self.soil_layer = soil_layer
        else:
            raise NotImplementedError

        self.ql_transpiration = ql_transpiration

    def calc_flux(self, Isotopologue, **kwargs):

        return self.ql_transpiration

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):

        """
        Positive flux of the given Isotopologue [m/day] left to right
        """

        return self.calc_flux(Isotopologue=Isotopologue) * self.soil_layer.get_conc_iso_liquid(
            Isotopologue=Isotopologue)


class surface_runoff(boundary_connection):
    """
    Description
    ===========
    Flux of isotopes [kg] between two nodes (atmosphere and upper soil layer) due to evaporation.

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 atmosphere,
                 top_layer,
                 q_runoff=0.0
                 ):

        """
        Constructor of evaporation

        """
        if isinstance(atmosphere, iso_storages.iso_atmosphere) and isinstance(top_layer, iso_storages.iso_soil_layer):

            boundary_connection.__init__(self, left_node=atmosphere, right_node=top_layer)
            self.top_layer = top_layer
        else:
            raise NotImplementedError

        self.q_runoff = q_runoff

    def calc_flux(self, Isotopologue, **kwargs):

        return self.q_runoff

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):

        return self.calc_flux(Isotopologue=Isotopologue) * self.top_layer.get_conc_iso_liquid(Isotopologue=Isotopologue)


class precipitation(boundary_connection):
    """
    Description
    ===========
    Flux of isotopes [kg] between two nodes (atmosphere and upper soil layer) due to evaporation.

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 atmosphere,
                 top_layer,
                 q_prec=0.0,
                 c_prec={"2H": 1.0, "18O": 1.0}
                 ):

        """
        Constructor of evaporation

        """
        if isinstance(atmosphere, iso_storages.iso_atmosphere) and isinstance(top_layer, iso_storages.iso_soil_layer):

            boundary_connection.__init__(self, left_node=atmosphere, right_node=top_layer)

        else:
            raise NotImplementedError

        self.q_prec = q_prec
        self.ci_prec = c_prec

    def calc_flux(self, Isotopologue, **kwargs):
        return - self.q_prec

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):

        return self.calc_flux(Isotopologue=Isotopologue) * self.ci_prec[Isotopologue]

    def set_conc_prec(self, c_prec, Isotopologue):

        self.ci_prec[Isotopologue] = c_prec


class aquifer_connection(boundary_connection):
    """
    Description
    ===========
    Flux of isotopes [kg] due to transpiratation (root extraction) between a soil storage and the atmosphere based on a given liquid flux [m/day]

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 soil_layer,
                 aquifer,
                 ql_layer=0.0,  # positive liquid flux in (kg/(m2*day)) from node_SoilLayer to node_aquifer.
                 ):
        """
        Constructor of boundary_connection

        """

        if isinstance(aquifer, iso_storages.iso_aquifer) and isinstance(soil_layer, iso_storages.iso_soil_layer):

            boundary_connection.__init__(self, left_node=soil_layer, right_node=aquifer)
            self.soil_layer = soil_layer
            self.aquifer = aquifer
        else:
            raise NotImplementedError

        self.ql = ql_layer

    def get_flux(self):
        return self.ql

    def calc_flux(self, Isotopologue, **kwargs):

        if self.aquifer.get_conc_iso_liquid(Isotopologue=Isotopologue) > 0.0:
            flux = self.get_flux() * 1 - self.get_flux() * 1  # sli.solve.f90:L2594
        else:
            flux = self.get_flux() * 1

        return flux

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):

        if self.aquifer.get_conc_iso_liquid(Isotopologue=Isotopologue) > 0.0:
            # sli.solve.f90:L2594
            flux_i = self.get_flux() * self.aquifer.get_conc_iso_liquid(Isotopologue=Isotopologue)
        else:
            flux_i = self.get_flux() * self.soil_layer.get_conc_iso_liquid(Isotopologue=Isotopologue)

        return flux_i


class neuman_boundary(boundary_connection):
    """
    Description
    ===========
    Flux of isotopes [kg] between two nodes (atmosphere and usoil layer).

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 atmosphere,
                 soil_layer,
                 q_neuman=0.0,
                 c_neuman={"2H": 1.0, "18O": 1.0}
                 ):

        """
        Constructor of evaporation

        """
        if isinstance(atmosphere, iso_storages.iso_atmosphere) and isinstance(soil_layer, iso_storages.iso_soil_layer):

            boundary_connection.__init__(self, left_node=atmosphere, right_node=soil_layer)

        else:
            raise NotImplementedError

        self.q_neuman = q_neuman
        self.ci_neuman = c_neuman

    def calc_flux(self, Isotopologue, **kwargs):
        return self.q_neuman

    def calc_flux_liquid(self, Isotopologue, **kwargs):
        return self.calc_flux(Isotopologue=Isotopologue, **kwargs)

    def calc_flux_i(self, Isotopologue, **kwargs):

        return self.calc_flux(Isotopologue=Isotopologue) * self.ci_neuman[Isotopologue]

    def set_conc_neuman(self, c_neuman, Isotopologue):

        self.ci_neuman[Isotopologue] = c_neuman
