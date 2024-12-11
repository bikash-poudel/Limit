
class _2D_Soil(object):

    def __init__(self, path):
        self.__path = path
        self.__qevap = self.__read_2dSoil(self.__path + '/q_evap.txt')
        self.__temp = self.__read_2dSoil(self.__path + '/temp.txt')
        self.__theta = self.__read_2dSoil(self.__path + '/theta.txt')
        self.__head = self.__read_2dSoil(self.__path + '/head.txt')
        self.__ql = self.__read_2dSoil(self.__path + '/q_l.txt')
        self.__qv = self.__read_2dSoil(self.__path + '/q_v.txt')

    def get_path(self):
        return self.__path

    def __read_2dSoil(self, path):
        # Read the text file
        with open(path, 'r') as file:
            # Read all lines into a list
            lines = file.readlines()

        # Split each line by spaces and flatten the resulting list
        data = [line.strip().split() for line in lines]

        return data

    ###### Evap #######
    def get_qevap(self):
        return self.__qevap[0]

    def qevap(self, time):
        return float(self.get_qevap()[time]) / 100  # cm/s into m/s

    ###### Head #######

    def get_head(self):
        hd = self.__head
        head = list(map(list, zip(*hd)))

        return head

    def head(self, time):
        return [float(hd) / 100 for hd in self.get_head()[time]][::-1]  # cm to m

    ###### Theta #######

    def get_theta(self):
        th = self.__theta
        theta = list(map(list, zip(*th)))

        return theta

    def theta(self, time):
        return [float(th) for th in self.get_theta()[time]][::-1]

    ###### Temperature #######

    def get_temp(self):
        T = self.__temp
        Temp = list(map(list, zip(*T)))

        return Temp

    def T(self, time):
        return [float(tm) for tm in self.get_temp()[time]][::-1]

    ###### Liquid Fluxes #######

    def get_flux_l(self):
        ql = self.__ql
        q_l = list(map(list, zip(*ql)))

        return q_l

    def q_l(self, time):
        return [float(q)/100 for q in self.get_flux_l()[time]][::-1]  # cm/s into m/s

    ###### Vapor Fluxes #######

    def get_flux_v(self):
        qv = self.__qv
        q_v = list(map(list, zip(*qv)))

        return q_v

    def q_v(self, time):
        return [float(q)/100 for q in self.get_flux_v()[time]][::-1]  # cm/s into m/s

