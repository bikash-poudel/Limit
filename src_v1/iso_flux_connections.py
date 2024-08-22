""""

"""

# -*- coding: utf-8 -*-

# from iso_Storages import Storage
import math

class linear_connection(object):
    """
    Description
    ===========
    Creates the flux connection [advective and diffusive] in both phases between the isotope storages.

    """

    def __init__(self,
                 storage1,  # Instance of iso-storage
                 storage2,  # Instance of iso-storage
                 ql,
                 qv,
                 Isotopologue,
                 time_steps
                 ):

        self.Isotopologue = Isotopologue
        self.q_l = ql
        self.q_v = qv
        self.storage1 = storage1
        self.storage2 = storage2
        self.__time_steps = time_steps

        self.liquid_advection = liquid_advection(self.q_l, self.storage1)
        self.liquid_diffusion = liquid_diffusion(self.storage1)
        self.vapor_advection = vapor_advection(self.q_v, self.storage1)
        self.vapor_diffusion = vapor_diffusion(self.storage1)

    def run(self):

        c1 = []  # liquid isotope concentration for storage1
        c2 = []  # liquid isotope concentration for storage2

        v1 = []
        v2 = []

        c1.append(self.storage1.get_conc_iso_liquid(self.Isotopologue))
        c2.append(self.storage2.get_conc_iso_liquid(self.Isotopologue))
        v1.append(self.storage1.get_conc_iso_vapor(self.Isotopologue))
        v2.append(self.storage2.get_conc_iso_vapor(self.Isotopologue))

        for time in range(self.__time_steps):

            self.liquid_advection = liquid_advection(self.q_l, self.storage1)
            self.liquid_diffusion = liquid_diffusion(self.storage1)
            self.vapor_advection = vapor_advection(self.q_v, self.storage1)
            self.vapor_diffusion = vapor_diffusion(self.storage1)

            # Total outgoing flux is the sum of advective flux and diffusive flux
            liquid_iso_outflux = self.liquid_advection.calc_flux_i(self.Isotopologue) + \
                                            self.liquid_diffusion.calc_flux_i(self.Isotopologue)
            # TODO: vapor advective / diffusive flux due to the concentration gradient / including equilibrium fractionation ??
            vapor_iso_outflux = self.vapor_advection.calc_flux_i(self.Isotopologue) + \
                                           self.vapor_diffusion.calc_flux_i(self.Isotopologue)

            liquid_iso_influx = liquid_iso_outflux
            vapor_iso_influx = vapor_iso_outflux

            # Calculate new concentration in the storage after the incoming or outgoing flux for each time steps
            new_liquid_conc_storage1 = (self.storage1.get_conc_iso_liquid(self.Isotopologue) * self.storage1.volume -
                                        liquid_iso_outflux * self.storage1.area) / self.storage1.volume
            new_liquid_conc_storage2 = (self.storage2.get_conc_iso_liquid(self.Isotopologue) * self.storage2.volume +
                                        liquid_iso_influx * self.storage2.area) / self.storage2.volume

            new_vapor_conc_storage1 = (self.storage1.get_conc_iso_vapor(self.Isotopologue) * self.storage1.volume -
                                       vapor_iso_outflux * self.storage1.area) / self.storage1.volume
            new_vapor_conc_storage2 = (self.storage2.get_conc_iso_vapor(self.Isotopologue) * self.storage2.volume +
                                       vapor_iso_influx * self.storage2.area) / self.storage2.volume

            # Update each storage with new liquid and vapor concentration
            self.storage1.set_conc_iso_liquid(new_liquid_conc_storage1, self.Isotopologue)
            self.storage2.set_conc_iso_liquid(new_liquid_conc_storage2, self.Isotopologue)

            # Store isotope concentration for each time steps
            c1.append(new_liquid_conc_storage1)
            c2.append(new_liquid_conc_storage2)
            v1.append(new_vapor_conc_storage1)
            v2.append(new_vapor_conc_storage2)

        return [[c1, c2], [v1,v2]]


class liquid_advection(object):
    """
    Description
    ===========
    Advective flux of isotopes [kg] between to nodes based on a given liquid flux [m/day]

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 q_l,  # Liquide flux in m/day positive form left to right : negative from right to left
                 storage # instance of iso_storage
                 ):
        """
        Constructor of boundary_connection

        @param q_l: Liquide flux in m3/(m2*day). Positive form left to right : negative from right to left
        @type q_l: float
        """
        # flux_connection.__init__(self ,left_node ,right_node ,flux_crosssectional_area)
        self.q_l = q_l
        self.storage = storage

    def calc_flux(self, Isotopologue):
        """
        Positive flux of the liquid flux  [m/day] left to right
        """
        flux = self.q_l

        return flux

    def calc_flux_i(self, Isotopologue):
        """
        Positive flux of the given Isotopologue [kg/m2/day] in liquid phase left to right
        """
        flux = self.calc_flux(Isotopologue)

        flux_i = self.storage.get_conc_iso_liquid(Isotopologue) * flux

        return flux_i

    def calc_flux_liquid(self, Isotopologue):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux
        """
        return self.calc_flux(Isotopologue)


class liquid_diffusion(object):
    """
    Description
    ===========
    Functions to calculate diffusive fluxes in liquids
    """

    def __init__(self, storage, q_l=None, hydrodynamic_dispersivity=0.0):

        self.storage = storage
        self.q_l = q_l
        self.hydrodynamic_dispersivity = hydrodynamic_dispersivity

    def dl_i_eff(self, dl_i, q_l=None, hydrodynamic_dispersivity=0.0):
        """
        Calculates the total (effective) liquid diffusivity for isotope i in the soil liquide phase (m**2/day)

        Description
        ===========

        Under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones and are frequently set to zero (Auriault and Adler, 1995)

        @param theta: liquid phase (m**3/m**3 )
        @type theta: Float

        @param tortuosity: tortuosity (m/m)
        @type tortuosity: Float

        @param q_l: liqiud water flux (m/s)
        @type q_l: Float

        @param hydrodynamic_dispersivity: per default set to 0 (m) because under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones
        @type hydrodynamic_dispersivity: Float
        """
        try:
            if q_l == None or hydrodynamic_dispersivity == 0.0:
                dl_i_eff = dl_i
            elif q_l != None and hydrodynamic_dispersivity != 0.0:
                dl_i_eff = dl_i + hydrodynamic_dispersivity * abs(q_l)

            return dl_i_eff

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def dl_i(self, T, Isotopologue, tortuosity, theta, formulation="Cuntz", ignore_dl_i=False):
        """
        Calculates the liquid diffusivity for isotope i in the soil liquid phase (m**2/day)

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
        try:
            if ignore_dl_i == False:
                # SLI: cable_sli_solve.f90::L2381-2390
                dl_i = self.dl_i_self_diffusivity(T, Isotopologue, formulation="Cuntz") * theta * tortuosity
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

                dl_i = a_i * 1.10 ** -9 * math.exp(-(535400 / T ** 2) + (1393.3 / T) + 2.1876)

            # SLI: cable_sli_solve.f90::L2381-2390
            elif formulation == "Cuntz":
                if "18O" == Isotopologue:
                    a_i = 1 / 1.026
                elif "2H" == Isotopologue:
                    a_i = 1 / 1.013
                else:
                    raise NotImplementedError

                dl_i = a_i * 1e-7 * math.exp(-577 / (T - 145))

            else:
                raise NotImplementedError

            return dl_i * 86400  # convert m2 per seconds to m2 per day

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def calc_flux(self, Isotopologue, dl_i_formulation="Cuntz", ignore_dl_i=False):
        """
        Positive flux of the given Isotopologue [m/day] left to right
        Average diffusivity over the thickness
        """
        # Under evaporation, convective and hydrodynamic dispersion processes are negligible as compared to the diffusion ones and are frequently set to zero (Auriault and Adler, 1995)
        # left_node = self.left_node()
        # right_node = self.rigth_node()

        dl_i = self.dl_i(T=self.storage.T,  # TODO: where is the temperature imported from??
                         Isotopologue=Isotopologue,
                         tortuosity=self.storage.tortuosity,
                         theta=self.storage.theta,
                         formulation=dl_i_formulation,
                         ignore_dl_i=ignore_dl_i)

        dl_i_eff = self.dl_i_eff(dl_i=dl_i, q_l=self.q_l,
                                 hydrodynamic_dispersivity=self.hydrodynamic_dispersivity)

        delta_z = self.storage.thickness

        dl_mean = dl_i_eff * delta_z

        return dl_mean

    def calc_flux_i(self, Isotopologue):
        """
        Positive flux of the given Isotopologue [kg/m2/day] left to right
        """
        flux = self.calc_flux(Isotopologue)
        flux_i = self.storage.get_conc_iso_liquid(Isotopologue) * flux

        return flux_i


class vapor_advection(object):
    """
    Description
    ===========
    Advective flux of isotopes [kg] between to nodes based on a given vapor flux [m/day]

    A positive flux rate [kg] is always considered to be from left node to right node a  negative from right to left.

    """

    def __init__(self,
                 q_v, # vapor flux in m/day positive form left to right : negative from right to left. If None q_v  will be calculated based on equation A.12 Haverd&Cuntz (2005) knowing futurue Temp and matrix potential
                 storage  #  Instance of iso storage
                 ):
        """
        Constructor of boundary_connection

        @param q_v: vapor flux in (m/day). Positive form left to right : negative from right to left
        @type q_v: float
        """
        # flux_connection.__init__(self ,left_node ,right_node ,flux_crosssectional_area)
        self.q_v = q_v
        self.storage = storage

    def calc_flux(self, Isotopologue):
        """
        Positive flux of vapor  [m/day] left to right
        """
        flux = self.q_v

        return flux

    def calc_flux_i(self, Isotopologue):
        """
        Positive flux of the given Isotopologue [kg/m2/day] in vapor left to right
        """
        flux = self.calc_flux(Isotopologue)

        flux_i = self.storage.get_conc_iso_vapor(Isotopologue) * flux

        return flux_i

    def calc_flux_vapour(self, Isotopologue):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux
        """
        flux = self.calc_flux(Isotopologue)

        # TODO: Is this correct way to calculate vapour flux. (not isotope)
        flux_vapour = self.storage.alpha_i(Isotopologue) * flux

        return flux_vapour


class vapor_diffusion(object):
    """
    Description
    ===========
    Functions to caclculate diffusive fluxes
    """

    def __init__(self, storage, q_v=None, hydrodynamic_dispersivity=0.0):

        self.storage = storage
        self.q_v = q_v
        self.hydrodynamic_dispersivity = hydrodynamic_dispersivity

    @classmethod
    def dv_i(self, dv, Isotopologue, ignore_dv_i=False):
        """
        Calculates the vapour diffusivity of Isotopologue in free air (m**2/s)

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
        try:
            # SLI: cable_sli_solve.f90::L2366 - L2367
            if "18O" == Isotopologue:
                b_i = 1.0251  # H218O diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
            elif "2H" == Isotopologue:
                b_i = 1.0285  # HDO diffusivity in air (Merlivat 1978) // SLI: 1/b_i = alphak_vdiff
            else:
                raise NotImplementedError

            if ignore_dv_i == True:
                return dv
            else:
                # SLI: cable_sli_solve.f90::L2373 - L2374
                dv_i = dv / b_i  # SLI: dv = Dvs ; dv_i = Divs ; 1/b_i = alphak_vdiff

                return dv_i

        except ValueError as err:
            print(err)
            raise NotImplementedError

    @classmethod
    def dv_soil_air(self, dv_free_air, theta, theta_sat, tortuosity):
        """
        Calculates the vapour diffusivity of water in soil air space (m**2/s)

        Description
        ===========

        @param dv_free_air: vapour diffusivity of water in free air (m**2/s)
        @type dv_free_air: Float

        @param theta: air pressure in the atmosphere (1bar = 1 atm = 10**5 Pa) [Pa]
        @type theta: Float

        @param theta: liquid phase (m**3/m**3 )
        @type theta: Float

        @param theta_sat: theta at saturation equal to soil porosity (m**3/m**3 )
        @type theta_sat: Float

        @param tortuosity: tortuosity (m/m)
        @type tortuosity: Float
        """
        try:
            # SLI:: cable_sli_utils.f90 L::1494: var%Dv    = Dva*parin%tortuosity*(parin%the-theta)  * ((Tsoil+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
            dv_soil_air = dv_free_air * tortuosity * (theta_sat - theta)

            return dv_soil_air

        except ValueError as err:
            print(err)
            raise NotImplementedError

    @classmethod
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
            dva = 2.17e-5  # vapour diffusivity of water in air at 0 degC [m2/s]
            dv = dva * 1e5 / Pa * (T / 273.16) ** 1.88
            # SLI: Dva       = 2.17e-5   ! vapour diffusivity of water in air at 0 degC [m2/s]
            # SLI: Dvs  = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffuxivity of water in air (m2s-1)

            return dv

        except ValueError as err:
            print(err)
            raise NotImplementedError

    def calc_flux(self, Isotopologue, ignore_dv_i=False):
        """
        Positive flux of the given Isotopologue [m/day] left to right
        """
        # left_node = self.left_node()
        # right_node = self.rigth_node()

        # if isinstance(left_node, LSI_storages.iso_soil_layer) and isinstance(right_node, LSI_storages.iso_soil_layer):

        dv_free_air = self.dv_free_air(T=self.storage.T,
                                       Pa=1e5)  # TODO: Change to actual atmospheric preassure)

        dv_soil_air = self.dv_soil_air(dv_free_air=dv_free_air,
                                       theta=self.storage.theta,
                                       theta_sat=self.storage.theta_sat,
                                       tortuosity=self.storage.tortuosity)

        dv_i = self.dv_i(dv=dv_soil_air,
                         Isotopologue=Isotopologue,
                         ignore_dv_i=ignore_dv_i)

        delta_z = self.storage.thickness  # sqrt of area in m

        dv_i_mean = dv_i * delta_z * 86400

        return dv_i_mean

    def calc_flux_i(self, Isotopologue):
        """
        Positive flux of the given Isotopologue [kg/m2/day] left to right
        """
        flux = self.calc_flux(Isotopologue)
        flux_i = self.storage.get_conc_iso_vapor(Isotopologue) * flux

        return flux_i

    def calc_flux_liquid(self, Isotopologue):
        """
        Returns the flux relative to the liquid isotope concentration.
        For vapor flux --> vapor flux * alpha(i)
        For liquid flux --> liquid flux
        """
        flux = self.calc_flux(Isotopologue)
        flux_vapour_liquid = self.storage.alpha_i(Isotopologue) * flux

        return flux_vapour_liquid


