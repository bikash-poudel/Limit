
'''
Created on 20.08.2023
@author: poudel-b
'''

# -*- coding: utf-8 -*-

from math import exp

class iso_layer(object):
    """
    Layer storing isotopes

    TODO: move all connections and fluxes to the layer (Evaporation, Transpiration, Extraction, Seepage) as e.g. Neumann boundary condition

    Description
    ===========
      The
    """

    def __init__(self,
                 upper_boundary,
                 lower_boundary,
                 initial_c_solutes={"2H": 1.0, "18O": 1.0},
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 initial_theta=0.30,  # volumetric water content m3/m3
                 theta_0=0.01,  # volumetric water content m3/m3 at high suctions
                 theta_sat=0.35,  # volumetric water content m3/m3 at saturation
                 initial_T=283.15,  # soil temperature in Kelvin
                 porosity=0.35,  # porosity of the soil m3/m3
                 tortuosity=0.67,  # tortuosity of the soil m/m
                 psi=None  # soil metric potential in [m]
                 ):
        """
        Constructor of LSI_iso_cell

        @param cmf_source_cell: Pattern for the LSI_iso_cell
        @type cmf_source_cell: cmf.cell
        """
        self.lower_boundary = lower_boundary
        self.upper_boundary = upper_boundary
        self.area = 1  # None  # area of the layer in m2
        self.center = upper_boundary + (lower_boundary - upper_boundary) / 2
        self.thickness = lower_boundary - upper_boundary
        self.c_solutes = initial_c_solutes
        self.theta = initial_theta
        self.theta_0 = theta_0
        self.theta_sat = theta_sat
        self.T = initial_T
        self.porosity = porosity
        self.tortuosity = tortuosity
        self.psi = psi

    @classmethod
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



    def delta_theta(self, ql):
        """
        Returns the change in theta based on the given flux in m3

        @param ql: Liquide flux leavin or entering the current layer in m3
        @type ql: float
        """

        delta_theta = ql / (self.area * self.thickness)

        return delta_theta

    def delta_theta_bikash(self, ql_in, ql_out):

        delta_theta = (ql_in - ql_out) * self.porosity

        return delta_theta

