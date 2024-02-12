'''
Created on 22.08.2023

@author: poudel-b
'''
# -*- coding: utf-8 -*-

from math import exp


class BoundaryCondition(object):
    def __init__(self):

        self.upper_boundary_type = None
        self.upper_boundary_content = None
        self.lower_boundary_type = None
        self.lower_boundary_content = None

    def upper_boundary(self, boundary_type, content):
        self.upper_boundary_type = boundary_type
        self.upper_boundary_content = content

    def lower_boundary(self, boundary_type, content):
        self.lower_boundary_type = boundary_type
        self.lower_boundary_content = content

    def apply_boundary(self):

        pass


class iso_atmosphere(BoundaryCondition):
    """
    Layer storing isotopes

    Description
    ===========
      The
    """

    def __init__(self,
                 initial_c_atmosphere={"2H": 1.0, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 initial_T_atmosphere=283.15,  # temperature in the atmosphere in Kelvin
                 initial_Rh_atmosphere=0.20,  # relative humidity of the atmosphere (-)
                 initial_Pa=10 ** 5,  # atmospheric pressure (Pa)
                 initial_wind_speed=2.0,  # wind speed at top of canopy in m/s
                 initial_Precipitation=2.0,  # Precipitation in m3/d
                 initial_c_Precipitation={"2H": 1.0, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 ):
        """
        Constructor of LSI_iso_cell

        @param cmf_source_cell: Pattern for the LSI_iso_cell
        @type cmf_source_cell: cmf.cell
        """
        self.c_atmosphere = initial_c_atmosphere  # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
        self.T_atmosphere = initial_T_atmosphere  # temperature in the atmosphere in Kelvin
        self.Rh_atmosphere = initial_Rh_atmosphere  # porosity of the soil m3/m3
        self.Pa = initial_Pa  # Atmospheric preassure (Pa)
        self.wind_speed = initial_wind_speed  # wind speed at top of canopy in m/s
        self.Precipitation = initial_Precipitation  # Precipitation in m3/d
        self.c_Precipitation = initial_c_Precipitation  # Dict with solutes and initial concentrations in kg/m**3 of precipitation (!!NO delta signature!!) (currently supported "2H" and/or "18O")

    @classmethod
    def soil_surface_Rh(self, Rh_atmosphere, T_soil, T_atmosphere):
        """
        Returns the relative humidity at the soil surface

        Rh_atmosphere = relative humidity in the atmosphere
        T_soil = Temperature in Kelvin at the first soil layer ???
        pv_sat_soil = saturated water vapor volumetric mass of the soil (kg/m**3)
        T_atmosphere = Temperature in Kelvin in the atmosphere
        pv_sat_atmosphere = saturated water vapor volumetric mass of the atmosphere (kg/m**3)
        """

        soil_surface_Rh = Rh_atmosphere * self.pv_sat(T=T_atmosphere) / self.pv_sat(T=T_soil)

        return soil_surface_Rh

    @classmethod
    def cv_a(self, Rh, T):
        """
        Converts the Rh to concentration of the vapor in the air

        T = Temperature in Kelvin
        Rh = Relative humidity

        Control status: Checked on 15.05.2013 --> Results are the same as for SLIsli_main.f90::L205
        """

        # SLI: sli_main.f90::L205
        # vmet%cva  = vmet(:)%rha * esat(vmet(:)%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet(:)%Ta+Tzero)

        cv_a = Rh * self.p_sat(T) * 0.018 / (1000 * 8.314 * T)

        return cv_a

    @classmethod    # bikash
    def civ_a(self, solute_i):
        """
        Returns the concentration of isotope i in the atmosphere with respect to the current relative humidity
        """

        # SLI: sli_main.f90::L206
        # vmet%civa = Ra * vmet%cva * rhow * Mi/Mw18  ! kg HDO (liq) m-3 (air) !JLU David: civa = concentration of isotope in the atmosphere

        # c_atmosphere = float(self.c_atmosphere[solute_i])  # bikash
        civ_a = self.cv_a(Rh=self.Rh_atmosphere, T=self.T_atmosphere) * self.c_atmosphere[solute_i]

        return civ_a

    @classmethod
    def pv_sat(self, T):
        """
        Returns the saturated water vapor volumetric mass (kg/m**3)

        T_atmosphere = Temperature in Kelvin

        pv_sat = saturated water vapor volumetric mass (kg/m**3)

        TODO: Check for which temperature range it is valid

        Control status: Checked on 15.05.2013 --> Results are the same as for SLI
        """

        # SLI: sli_main.f90::L205
        pv_sat = 0.002166 * self.p_sat(T) / (T + 273.16)  # kg/m3

        return pv_sat

    @classmethod
    def p_sat(self, T):
        """
        Returns the Saturation vapour pressure, ps, in pascal based on Goff and Gratch and valid over the range 0 to 60 Degree celsius

        T_atmosphere = Temperature in Kelvin

        p_sat = saturated water vapor volumetric mass (kg/m**3)

        SLI: cable_sli_utils.f90 Lines: 2513-2525 <--checked and own implementation is okay for temp. above 273.15 Kelvin
        """

        # SLI: cable_sli_utils.f90::L2513-2525

        p_sat = 610.6 * exp(17.27 * (T - 273.15) / ((T - 273.15) + 237.3))  # Saturation vapour pressure, ps, in pascal

        if T < 273.16:  # below 0 Degree celsius
            p_sat = -4.86 + 0.855 * p_sat + 0.000244 * p_sat ** 2  # ps for ice

        return p_sat
