
import numpy as np


class moist(object):

    def __init__(self, path):
        self.__path = path
        self.__moist = self.__read_moist(self.__path + '/output_all.txt')

    def get_path(self):
        return self.__path

    def __read_moist(self, path):
        # Read the text file
        with open(path, 'r') as file:
            # Read all lines into a list
            lines = file.readlines()

        # Split each line by spaces and flatten the resulting list
        data = [line.strip().split() for line in lines]

        return data

    @property
    def all(self):
        return self.__moist

    @property
    def t(self):
        return [float(t[0]) for t in self.all]

    @property
    def dt(self):
        return np.insert(np.diff(self.t), 0, self.t[0])

    def head(self, time):
        return [float(h) for h in self.all[time][1:11]]

    def T(self, time):
        return [float(t) for t in self.all[time][11:21]]

    def theta(self, time):
        return [float(t) for t in self.all[time][31:41]]

    def Ta(self, time):
        return float(self.all[time][41])

    def Rnet(self, time):
        return float( self.all[time][42])

    def hra(self, time):
        return float(self.all[time][43])

    def Ts(self, time):
        return float(self.all[time][44])

    def qev(self, time):
        return float(self.all[time][45])

    def qls(self, time):
        return float(self.all[time][46])

    def qvs(self, time):
        return float(self.all[time][47])

    def ql(self, time):
        return [float(t) for t in self.all[time][48:59]]

    def qv(self, time):
        return [float(t) for t in self.all[time][59:70]]


class moist_v1(object):

    def __init__(self, path):
        self.__path = path
        self.__dt = self.__read_moist(self.__path + '/dt.txt')
        self.__Tatm = self.__read_moist(self.__path + '/Tatm.txt')
        self.__Rhatm = self.__read_moist(self.__path + '/rH_atm.txt')
        self.__Rnet = self.__read_moist(self.__path + '/Rnet.txt')
        self.__Tsur = self.__read_moist(self.__path + '/Tsur.txt')
        self.__qlsur = self.__read_moist(self.__path + '/qlsur.txt')
        self.__qvsur = self.__read_moist(self.__path + '/qvsur.txt')
        self.__qev = self.__read_moist(self.__path + '/qev.txt')
        self.__head = self.__read_moist(self.__path + '/head.txt')
        self.__theta = self.__read_moist(self.__path + '/theta.txt')
        self.__temp = self.__read_moist(self.__path + '/temp.txt')
        self.__ql = self.__read_moist(self.__path + '/ql.txt')
        self.__qv = self.__read_moist(self.__path + '/qv.txt')

    def get_path(self):
        return self.__path

    def __read_moist(self, path):
        # Read the text file
        with open(path, 'r') as file:
            # Read all lines into a list
            lines = file.readlines()

        # Split each line by spaces and flatten the resulting list
        data = [line.strip().split() for line in lines]

        return data

    ######### dt #########
    def get_dt(self):
        return self.__dt[0]

    def dt(self, time):
        return float(self.get_dt()[time])

    ######### atm #########
    def get_Tatm(self):
        return self.__Tatm[0]

    def get_Rhatm(self):
        return self.__Rhatm[0]

    def get_Rnet(self):
        return self.__Rnet[0]

    def T_atm(self, time):
        return float(self.get_Tatm()[time])

    def Rh_atm(self, time):
        return float(self.get_Rhatm()[time])

    def Rnet(self, time):
        return float(self.get_Rnet()[time])

    ######### Surface #########
    def get_Tsur(self):
        return self.__Tsur[0]

    def get_qlsur(self):
        return self.__qlsur[0]

    def get_qvsur(self):
        return self.__qvsur[0]

    def Tsur(self, time):
        return float(self.get_Tsur()[time])

    def qlsur(self, time):
        return float(self.get_qlsur()[time])

    def qvsur(self, time):
        return float(self.get_qvsur()[time])

    ######### Evap #########
    def get_qevap(self):
        return self.__qev[0]

    def q_ev(self, time):
        return float(self.get_qevap()[time])

    ######### States #########

    def get_theta(self):
        return self.__theta

    def get_tmp(self):
        return self.__temp

    def get_head(self):
        return self.__head

    def theta(self, time):
        return [float(th) for th in self.get_theta()[time]]

    def temp(self, time):
        return [float(t) for t in self.get_tmp()[time]]

    def head(self, time):
        return [float(h) for h in self.get_head()[time]]

    ######### Fluxes #########

    def get_ql(self):
        return self.__ql

    def get_qv(self):
        return self.__qv

    def ql(self, time):
        return [float(q) for q in self.get_ql()[time]]

    def qv(self, time):
        return [float(q) for q in self.get_qv()[time]]

path = r'D:\Isotope transport\soil water models\MOIST\8397416\thoritical_test\output'
m = moist(path)
