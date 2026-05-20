
import numpy as np


class SlI:

    def __init__(self, path):
        self.__path = path
        self.__scaler = self.__read_sli(self.__path + '/scaler.txt')
        self.__in_soil = self.__read_sli(self.__path + '/in_soil.txt')
        self.__in_variables = self.__read_sli(self.__path + '/in_variables.txt')
        self.__in_fluxes = self.__read_sli(self.__path + '/in_fluxes.txt')
        self.__in_parameter = self.__read_sli(self.__path + '/in_parameter.txt')
        self.__in_iso = self.__read_sli(self.__path + '/in_iso.txt')
        self.__inout_water = self.__read_sli(self.__path + '/inout_water.txt')
        self.__inout_iso = self.__read_sli(self.__path + '/inout_iso.txt')
        self.__out_iso = self.__read_sli(self.__path + '/out_iso.txt')

    def get_path(self):
        return self.__path

    def __read_sli(self, path):
        # Read the text file
        with open(path, 'r') as file:
            # Read all lines into a list
            lines = file.readlines()

        # Split each line by spaces and flatten the resulting list
        data = [line.strip().split() for line in lines]

        return data

    ###### in Scaler #######
    def get_scaler(self):
        return self.__scaler

    def scaler(self, time):
        return self.get_scaler()[time]

    def isotopologue(self, time):
        return self.scaler(time)[0]

    def testcase(self, time):
        return self.scaler(time)[1]

    def litter(self, time):
        return self.scaler(time)[2]

    def n(self, time):
        return int(self.scaler(time)[3])

    ###### in Soil ######
    def get_in_soil(self):
        return self.__in_soil

    def in_soil(self, time):
        return self.get_in_soil()[time]

    def ns(self, time):
        return int(self.in_soil(time)[0])

    def dx(self, time):
        return [float(soil) for soil in self.in_soil(time)[1:20]]

    def deltaz(self, time):
        return [float(soil) for soil in self.in_soil(time)[20:38]]

    def sig(self, time):
        return float(self.in_soil(time)[38])

    def dt(self, time):
        return float(self.in_soil(time)[39])

    def dxL(self, time):
        return float(self.in_soil(time)[40])

    ###### in variable ######
    def get_in_variables(self):
        return self.__in_variables

    def in_variables(self, time):
        return self.get_in_variables()[time]

    def T_soil0(self, time):
        return [float(temp) for temp in self.in_variables(time)[0:19]]

    def deltaT(self, time):
        return [float(dtm) for dtm in self.in_variables(time)[19:39]]

    def Sliqice(self, time):
        return [float(soil) for soil in self.in_variables(time)[39:58]]

    def deltaSliqice(self, time):
        return [float(soil) for soil in self.in_variables(time)[58:77]]

    def Sliq(self, time):
        return [float(soil) for soil in self.in_variables(time)[77:96]]

    def deltaSliq(self, time):
        return [float(soil) for soil in self.in_variables(time)[96:115]]

    def Sice(self, time):
        return [float(soil) for soil in self.in_variables(time)[115:134]]

    def deltaSice(self, time):
        return [float(soil) for soil in self.in_variables(time)[134:153]]

    def Ts(self, time):
        return float(self.in_variables(time)[153])

    def TL(self, time):
        return float(self.in_variables(time)[154])

    def T0(self, time):
        return float(self.in_variables(time)[155])

    def pond_height(self, time):
        return float(self.in_variables(time)[156])

    def dh0(self, time):
        return float(self.in_variables(time)[157])

    def SLliq(self, time):
        return float(self.in_variables(time)[158])

    def deltaSLliq(self, time):
        return float(self.in_variables(time)[159])

    ###### in fluxes ######
    def get_in_fluxes(self):
        return self.__in_fluxes

    def in_fluxes(self, time):
        return self.get_in_fluxes()[time]

    def q(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[0:20]]

    def ql(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[20:40]]

    def qv(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[40:60]]

    def qprec(self, time):
        return float(self.in_fluxes(time)[60])

    def qevap(self, time):
        return float(self.in_fluxes(time)[61])

    def qrunoff(self, time):
        return float(self.in_fluxes(time)[62])

    def qex(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[63:82]]

    def qd(self, time):
        return float(self.in_fluxes(time)[82])

    def qya(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[83:103]]

    def qyb(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[103:123]]

    def qTa(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[123:143]]

    def qTb(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[143:163]]

    def qvya(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[163:183]]

    def qvyb(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[183:203]]

    def qvTa(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[203:223]]

    def qvTb(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[223:243]]

    def dy(self, time):
        return [float(flux) for flux in self.in_fluxes(time)[243:263]]

    def qsig(self, time):

        q = np.array(self.q(time))

        dy = np.array(self.dy(time))
        qya = np.array(self.qya(time))
        qyb = np.array(self.qyb(time))
        qTa = np.array(self.qTa(time))
        qTb = np.array(self.qTb(time))
        dT = np.array(self.deltaT(time))

        qsig = q[1:-1] + (qya[1:-1] * dy[1:-1] + qyb[1:-1] * dy[2:] + qTa[1:-1] * dT[1:-1] + qTb[1:-1] * dT[2:])
        q0 = q[0] + (qya[0] * dy[0] + qyb[0] * dy[1] + qTa[0] * dT[0] + qTb[0] * dT[1])
        qn = q[-1] + (qya[-1] * dy[-1] + qTa[-1] * dT[-1])

        return np.concatenate(([q0], qsig, [qn]))

    def qvsig(self, time):

        qv = np.array(self.qv(time))

        dy = np.array(self.dy(time))
        qvya = np.array(self.qvya(time))
        qvyb = np.array(self.qvyb(time))
        qvTa = np.array(self.qvTa(time))
        qvTb = np.array(self.qvTb(time))
        dT = np.array(self.deltaT(time))

        qvsig = qv[1:-1] + (qvya[1:-1] * dy[1:-1] + qvyb[1:-1] * dy[2:] + qvTa[1:-1] * dT[1:-1] + qvTb[1:-1] * dT[2:])
        qv0 = qv[0] + (qvyb[0] * dy[1] + qvTb[0] * dT[1])
        qvn = 0

        return np.concatenate(([qv0], qvsig, [qvn]))

    def qlsig(self, time):

        qlsig = self.qsig(time)[1:-1] - self.qvsig(time)[1:-1]
        qln = self.qsig(time)[-1]

        return np.concatenate((qlsig, [qln]))

    def qevapsig(self, time):

        qev = self.qevap(time)
        qyb = self.qyb(time)[0]
        qTb = self.qTb(time)[0]

        dy = self.dy(time)[1]
        dT = self.deltaT(time)[1]

        return qev - (qyb * dy + qTb * dT)



    ###### in parameter ######
    def get_in_parameter(self):
        return self.__in_parameter

    def in_parameter(self, time):
        return self.get_in_parameter()[time]

    def var_cv(self, time):
        return [float(param) for param in self.in_parameter(time)[0:19]]

    def var_dv(self, time):
        return [float(param) for param in self.in_parameter(time)[19:38]]

    def theta(self, time):
        return [float(param) for param in self.in_parameter(time)[38:57]]

    def thetasat(self, time):
        return [float(param) for param in self.in_parameter(time)[57:76]]

    def thetar(self, time):
        return [float(param) for param in self.in_parameter(time)[76:95]]

    def tortuosity(self, time):
        return [float(param) for param in self.in_parameter(time)[95:114]]

    def deltacv(self, time):
        return [float(param) for param in self.in_parameter(time)[114:133]]

    def ram(self, time):
        return float(self.in_parameter(time)[133])

    def rbh(self, time):
        return float(self.in_parameter(time)[134])

    def cva(self, time):
        return float(self.in_parameter(time)[135])

    def civa(self, time):
        return float(self.in_parameter(time)[136])

    def thetasatL(self, time):
        return float(self.in_parameter(time)[137])

    def deltacvL(self, time):
        return float(self.in_parameter(time)[138])

    def R_humidity(self, time):
        return [float(param) for param in self.in_parameter(time)[139:158]]

    def matric_pot(self, time):
        return [float(param) for param in self.in_parameter(time)[158:177]]

    def R_humidity_atm(self, time):
        return float(self.in_parameter(time)[177])

    def Tatm(self, time):
        return float(self.in_parameter(time)[178])

    def wind_speed(self, time):
        return float(self.in_parameter(time)[179])

    def E_pot(self, time):
        return float(self.in_parameter(time)[180])

    def time_step(self, time):
        return float(self.in_parameter(time)[181])


    ###### in iso ######
    def get_in_iso(self):
        return self.__in_iso

    def in_iso(self, time):
        return self.get_in_iso()[time]

    def cprec(self, time):
        return float(self.in_iso(time)[0])

    def cali(self, time):
        return float(self.in_iso(time)[1])

    ###### inout water ######
    def get_inout_water(self):
        return self.__inout_water

    def inout_water(self, time):
        return self.get_inout_water()[time]

    def ql0(self, time):
        return float(self.inout_water(time)[0])

    def qv0(self, time):
        return float(self.inout_water(time)[1])

    ###### inout iso ######
    def get_inout_iso(self):
        return self.__inout_iso

    def inout_iso(self, time):
        return self.get_inout_iso()[time]

    def ciso(self, time):
        return [float(iso) for iso in self.inout_iso(time)[0:20]]

    def cisoice(self, time):
        return [float(iso) for iso in self.inout_iso(time)[20:40]]

    def ciso0(self, time):
        return float(self.inout_iso(time)[40])

    def cisoL(self, time):
        return float(self.inout_iso(time)[41])

    def cisos(self, time):
        return float(self.inout_iso(time)[42])

    ###### out iso ######
    def get_out_iso(self):
        return self.__out_iso

    def out_iso(self, time):
        return self.get_out_iso()[time]

    def qiso_in(self, time):
        return float(self.out_iso(time)[0])

    def qiso_out(self, time):
        return float(self.out_iso(time)[1])

    def qiso_evap(self, time):
        return float(self.out_iso(time)[2])

    def qiso_trans(self, time):
        return float(self.out_iso(time)[3])

    def qiso_liq_adv(self, time):
        return [float(iso) for iso in self.out_iso(time)[4:23]]

    def qiso_vap_adv(self, time):
        return [float(iso) for iso in self.out_iso(time)[23:42]]

    def qiso_liq_diff(self, time):
        return [float(iso) for iso in self.out_iso(time)[42:60]]

    def qiso_vap_diff(self, time):
        return [float(iso) for iso in self.out_iso(time)[60:78]]





