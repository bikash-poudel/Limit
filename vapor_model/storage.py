'''
Created on 11.12.2024
@author: poudel-b
'''

# -*- coding: utf-8 -*-

from math import exp, sqrt


class flux_node(object):
    def __init__(self):

        self.__Connections = []
        self.__Connections_to_storages = []
        self.__Connections_to_boundaries = []
        self.__Connections_to_right = []
        self.__Connections_to_left = []

    def get_connections(self):
        return self.__Connections

    connections = property(get_connections, None, None, "List of all flux connections connected with this node")

    def get_connections_to_boundaries(self):
        return self.__Connections_to_boundaries

    connections_to_boundaries = property(get_connections_to_boundaries, None, None,
                                         "List of all flux connections between this node and boundaries (like atmopshere)")

    def get_connections_to_storages(self):
        return self.__Connections_to_storages

    connections_to_storages = property(get_connections_to_storages, None, None,
                                           "List of all flux connections between this node and iso storages")

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
            if isinstance(newConnection.right_node, soil_layer):
                # left connection only to the storages, to ta boundaries are connection to the boundaries
                self.__Connections_to_right.append(newConnection)
            other_node = newConnection.right_node
        elif newConnection.right_node == self:
            if isinstance(newConnection.left_node, soil_layer):
                # left connection only to the storages, to ta boundaries are connection to the boundaries
                self.__Connections_to_left.append(newConnection)
            other_node = newConnection.left_node

        else:
            raise NotImplementedError

        if isinstance(other_node, soil_layer):
            self.__Connections_to_storages.append(newConnection)
        elif isinstance(other_node, atmosphere):
            self.__Connections_to_boundaries.append(newConnection)


class atmosphere(flux_node):

    def __init__(self,
                 # Dict with solutes and initial concentrations in kg/m**3 (!!NO delta signature!!) (currently supported "2H" and/or "18O")
                 T=283.15,  # temperature in the atmosphere in Kelvin
                 Rh_atmosphere=0.20,  # relative humidity of the atmosphere (-)
                 Pa_atmosphere=10 ** 5,  # atmospheric pressure (Pa)
                 initial_wind_speed=2.0,  # wind speed at top of canopy in m/s
                 R_net=0.0,  # net radiation absorbed
                 wind_speed=2.0,  # wind speed at top of canopy in m/s
                 hc=10,  # canopy height [m] (e.g. 40)
                 d0=0.67 * 10,  # displacement height (e.g. 0.7 * hc)
                 z0m=0.1 * 10,  # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
                 LAI=1.1,  # leaf area index (e.g. 2.0)
                 extku=1.5):  # extinction coeff't for windspeed (e.g. 0.5)

        flux_node.__init__(self)
        self.T = T
        self.Rh = Rh_atmosphere  # porosity of the soil m3/m3
        self.Pa = Pa_atmosphere  # Atmospheric preassure (Pa)
        self.R_net = R_net  # Net radiation absorbed
        self.wind_speed = wind_speed  # wind speed at top of canopy in m/s
        self.hc = hc  # canopy height [m]
        self.d0 = d0  # displacement height
        self.z0m = z0m  # roughness height for momentum (e.g. 0.1 * hc) need to be 0.0 if hc = 0.0
        self.LAI = LAI  # leaf area index (e.g. 2.0)
        self.extku = extku  # extinction coeff't for windspeed


class soil_layer(flux_node):
    """
    Description
    ===========
    Soil Layer storing the state of the:
    -water storage
    -temperature

    and the physical parameters of the soil:
    -tortuosity
    -theta_0
    """

    def __init__(self,
                 upper_boundary,
                 lower_boundary,
                 theta=0.547,
                 theta_sat=0.547,
                 tortuosity=0.67,
                 rH=0,
                 head=1.0,
                 T=25,
                 T0=25,
                 K_s=3.8e-3,  # cm /s
                 phi_e=-13.0,  # cm
                 b=6.53,
                 clay=0.249,
                 sand=0.022,
                 silt=0.729,
                 som=0.044,
                 rho_b=1.2,
                 rho_l=1.0
                 ):

        flux_node.__init__(self)
        self.upper_boundary = upper_boundary
        self.lower_boundary = lower_boundary
        self.thickness = lower_boundary - upper_boundary
        self.theta_sat = theta_sat
        self.tortuosity = tortuosity
        self.rH = rH
        self.T0 = T0

        # Soil Properties
        self.Ks = K_s
        self.phi = phi_e
        self.b = b
        self.clay = clay
        self.sand = sand
        self.silt = silt
        self.som = som
        self.rho_b = rho_b
        self.rho_l = rho_l

        self.__theta = theta
        self.__head = head
        self.__T = T  # K
        # self.__theta = 0.547

    """Properties"""

    @property
    def head(self):
        return self.__head

    @head.setter
    def head(self, h):
        self.__head = h

    @property
    def theta(self):
        return self.calc_theta(self.head)

    # @theta.setter
    #def theta(self, theta):
     #   self.__theta = theta

    @property
    def T(self):
        """in  Celcius"""
        return self.__T

    @T.setter
    def T(self, T):
        self.__T = T

    @property
    def kT0(self):

        """Calculate the hydraulic conductivity at T_0"""

        k = self.Ks * ((self.theta / self.theta_sat) ** (2 * self.b + 3))

        return k

    @property
    def kml(self):

        """The hydraulic conductivity of liquid water due to water potential (<0)
           For Kml, in unit of "cm s^-1"
           refer to Nassar and Horton 1997 and Heitman 2008"""

        T0 = self.T0

        muT = (1 / (0.021482 * (self.T - 8.435 + (self.T ** 2 - 16.87 * self.T + 8149.5492) ** 0.5) - 1.2)) / 1000
        muT_0 = (1 / (0.021482 * (T0 - 8.435 + (T0 ** 2 - 16.87 * T0 + 8149.5492) ** 0.5) - 1.2)) / 1000  # (N s/m^2)

        return self.kT0 * muT_0 / muT

    @property
    def dtl(self):

        """"The hydraulic conductivity(diffusivity) of liquid water due to temperature change
        - For KpPhi_mpT, in unit of "cm^2 s^-1 T^-1"
        - refer to Nassar and Horton 1997 and Heitman 2008
        """

        Sa = 2.4415e6  # cm-2/cm-3
        gain_factor = 4
        kk = Sa / 981 * 1.4550e-6

        return self.kml * kk * gain_factor

    @property
    def w_capacity(self):

        """The water capicity [cm-1], potention to Capa
        refer to Nassar and Horton 1997 and Heitman 2008 """

        h = min(self.head, -13.0)  # cm
        h = max(h, -3.0e10)  # cm

        c = - self.theta_sat * abs((h / self.phi) ** (-1 / self.b - 1)) / self.b / self.phi

        return c

    @property
    def lambda_bulk(self):

        """Thermal conductivity (lambda)
        - unit "w cm-1 K-1"
        """

        lambda_dry = - 0.56 * self.theta_sat + 0.51
        alpha = 0.67 * self.clay + 0.24
        beta = 1.97 * self.sand + 1.87 * self.rho_b - 1.36 * self.sand * self.rho_b - 0.95

        lam = 0.01 * (lambda_dry + exp(beta - self.theta ** (-alpha)))

        return lam

    @property
    def cs(self):

        """soil volumetric head capacity
        - unit J cm-3 K-1
        """

        # preliminary parameters for each partitions
        c_o = 1.900  # "J g^-1 K-1"
        c_m = 0.730  # "J g^-1 K-1"
        c_v = 1.864  # "J g^-1 K-1"
        c_l = 4.187  # "J g^-1 K-1"

        # weight of each partition
        phi_o = self.som
        phi_m = self.sand + self.silt + self.clay

        # potential - temperature correction
        ht = self.head
        Tt = self.T + 273.15
        theta_a = self.theta_sat - self.theta
        rho_vs = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # gm cm-3
        HR = exp(2.1238e-4 * ht / Tt)
        rho_v = rho_vs * HR  # gm cm-3

        c = self.rho_b * (c_o * phi_o + c_m * phi_m) + c_v * rho_v * theta_a + c_l * self.rho_l * self.theta

        return c

    @property
    def Dtv(self):

        """"The diffusivity of (conductivity of vapor due to temperature or potential)
        - unit "cm2 s-1 k-1"
        """

        rho_w = 1.0  # gm-3
        Tt = self.T + 273.15

        D = 0.23500 * abs((Tt / 273.15) ** 1.75)  # cm2 s-1
        omega = abs((self.theta_sat - self.theta) ** (2 / 3))
        theta_a = self.theta_sat - self.theta

        rho_s = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # g m-3
        HR = exp(2.1238e-4 * self.head / Tt)
        prho_spT = rho_s * 4975.9 / (Tt ** 2)

        pHRpT = -HR * 2.1238e-4 * self.head / (Tt ** 2)
        eta = 8 + 6 * self.theta - 7 * exp(-((1 + 2.6 / sqrt(self.clay)) * self.theta) ** 4)

        d = eta * D * omega * theta_a * (HR * prho_spT + rho_s * pHRpT) / rho_w

        return d

    @property
    def Dmv(self):

        """"The diffusivity of (conductivity of vapor due to temperature or potential)
              - unit "cm/s"
        """
        Tt = self.T + 273.15

        D = 0.23500 * abs((Tt / 273.15) ** 1.75)  # cm2 s-1
        omega = abs((self.theta_sat - self.theta) ** (2 / 3))
        theta_a = self.theta_sat - self.theta

        rho_s = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # g m-3
        HR = exp(2.1238e-4 * self.head / Tt)
        pHRph = HR * 2.1238e-4 / Tt  # cm-1

        d = D * omega * theta_a * (rho_s * pHRph) / self.rho_l

        return d

    @property
    def H1(self):

        """Hydraulic_1
        - we defined hera as H1
        - unit cm-1
        - H1 is not as same as the coefficient, butwe multiply it by water
        """

        rho_w = 1.00  # g cm-3
        theta_a = self.theta_sat - self.theta
        Tt = self.T + 273.15

        rho_vs = 1.0e-6 * exp(19.84 - 4975.9 / Tt) # g m-3
        HR = exp(2.1238e-4 * self.head / Tt)
        rho_v = rho_vs * HR
        prho_vpPhi_m = rho_v * 2.1238e-4 / Tt

        h1 = (1 - rho_v / rho_w) * self.w_capacity + (theta_a / rho_w) * prho_vpPhi_m

        return h1

    @property
    def H2(self):

        """Hydraulic_1
        - we defined hera as H1
        - unit oC-1
        - H1 is not as same as the coefficient, i.e. Hydraulic_2
        """

        rho_w = 1.00  # g cm-3
        theta_a = self.theta_sat - self.theta
        Tt = self.T + 273.15

        rho_s = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # g cm-3
        HR = exp(2.1238e-4 * self.head / Tt)
        # rho_v = rho_s * HR
        prho_SpT = rho_s * 4975.9 / (Tt ** 2)

        pHRpT = -HR * 2.1238e-4 * self.head / (Tt ** 2)
        prho_vpT = rho_s * pHRpT + prho_SpT * HR

        h2 = theta_a * prho_vpT / rho_w  # 1/K

        return h2

    @property
    def T1(self):

        """ Thermal_1
        - we defined hera as T1
        - unit J cm-3 K-1
        - T1 is not as same as the coefficient of the pTpt.
        """

        # preliminary parameters

        L0 = 2270   # J g-1
        cv = 1.864   # J g-1 K-1

        theta_a = self.theta_sat - self.theta
        Tt = self.T + 273.15

        rho_s = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # g cm-3
        HR = exp(2.1238e-4 * self.head / Tt)
        # rho_v = rho_s * HR
        prho_SpT = rho_s * 4975.9 / (Tt ** 2)

        pHRpT = -HR * 2.1238e-4 * self.head / (Tt ** 2)
        prho_vpT = rho_s * pHRpT + prho_SpT * HR

        t1 = self.cs + (L0 + cv * (self.T - self.T0)) * theta_a * prho_vpT

        return t1

    @property
    def T2(self):

        """Thermal_2
        - we defined hera as T1
        - unit J cm-3 cm-1
        """

        rho_w = 1.00  # g cm-3

        # preliminary parameters
        L0 = 2270   # J g-1
        cv = 1.864   # J g-1 K-1
        cl = 4.187   # J g-1 K-1

        w = 0.00462  # J g-1

        theta_a = self.theta_sat - self.theta
        Tt = self.T + 273.15

        rho_vs = 1.0e-6 * exp(19.84 - 4975.9 / Tt)  # g cm-3
        HR = exp(2.1238e-4 * self.head / Tt)
        prho_vpPhi_m = rho_vs * HR * 2.1238e-4 / Tt  # g cm-3 cm-1
        rho_v = rho_vs * HR

        t2 = (L0 + cv * (self.T - self.T0)) * theta_a * prho_vpPhi_m \
             + (cl * rho_w * (self.T - self.T0) - rho_w * w - cv * rho_v * (self.T - self.T0)
                - L0 * rho_v) * self.w_capacity

        return t2

    @property
    def latent(self):

        rho_l = 1  # g cm-3 density of water
        L_0 = 2270  # [J g-1]
        c_v = 1.864  # Specific Heat at Constant Temperature [J g-1 K-1]

        return rho_l * (L_0 + c_v * (self.T - self.T0))

    """Functions"""

    def calc_theta(self, head):
        h = min(head, -13)  # cm
        h = max(h, -3.0e10)  # cm

        th = abs((h / self.phi) ** (-1 / self.b)) * self.theta_sat

        return th

    def calc_head(self, theta):
        h = self.phi * abs((theta / self.theta_sat) ** (-self.b))
        h = max(h, -5.0e9)  # cm
        h = min(h, -13.0)  # cm

        return h

    def distance(self, to_point):

        if isinstance(to_point, soil_layer):

            d = (self.thickness + to_point.thickness) / 2

        else:
            raise NotImplementedError

        return d

