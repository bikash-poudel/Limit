

import os


class hydrus(object):

    def __init__(self, path):
        self.__path = path
        self.__profile = self.__read_hydrus(self.__path + '/Profile.out')
        self.__Nod_inf = self.__read_hydrus(self.__path + '/Nod_Inf.out')
        self.__Nod_Inf_V = self.__read_hydrus(self.__path + '/Nod_Inf_V.out')
        self.__Atm = self.__read_hydrus(self.__path + '/ATMOSPH.IN')

    def write_selector(self):
        write_selector_in_file(self.__path)

    def get_path(self):
        return self.__path

    def __read_hydrus(self, path):
        # Read the text file
        with open(path, 'r') as file:
            # Read all lines into a list, excluding empty lines
            lines = [line.strip() for line in file if line.strip()]

        # Split each line by spaces and flatten the resulting list
        data = [line.split() for line in lines]

        return data

    def extract_nodes(self, data, time):

        data_blocks = {}
        current_time = None
        current_block = []

        for line in data:
            if "Time:" in line:
                # Store the previous block if it exists
                if current_time is not None:
                    data_blocks[current_time] = current_block

                # Extract the new time
                current_time = float(line[1])
                current_block = []
            elif current_time is not None:
                # Append the line to the current data block
                current_block.append(line)

        # Add the last block if it wasn't added
        if current_time is not None:
            data_blocks[current_time] = current_block

        # Return the data for the specified time
        return data_blocks.get(time, None)

    @property
    def time(self):

        t = []
        for line in self.__Nod_inf:
            if "Time:" in line:
                t.append(float(line[1]))

        return t[1:]  # Ignoring the time stamp


    ###### Profile #######
    @property
    def profile(self):
        return self.__profile[6:-1]

    @property
    def depth(self):
        d = list(zip(*self.profile))[1]
        return [float(x) / 100 for x in d][1:]  # covert [cm] to [m]

    @property
    def thetar(self):
        return list(map(float, list(zip(*self.profile))[2]))[1:]

    @property
    def theta_sat(self):
        return list(map(float, list(zip(*self.profile))[3]))[1:]

    @property
    def Ks(self):

        K = list(zip(*self.profile))[5]

        return [float(x) / 8640000 for x in K][1:]  # covert [cm day-1] to [m s-1]


    ###### Nod_inf #######
    @property
    def nod_inf(self):
        return self.__Nod_inf

    def node_output(self, time):
        return self.extract_nodes(self.nod_inf, time)[2:-1]

    def head(self, time):

        h = list(zip(*self.node_output(time)))[2]
        return [float(x) / 100 for x in h][1:]  # covert [cm] to [m]

    def theta(self, time):
        return list(map(float, list(zip(*self.node_output(time)))[3]))[1:]

    def transpiration(self, time): #TODO: check rtanspiratin before implementing
        return list(map(float, list(zip(*self.node_output(time)))[7]))[1:]

    def T(self, time):
        return list(map(float, list(zip(*self.node_output(time)))[10]))[1:]


    ###### Nod_inf_v #######
    @property
    def node_inf_v(self):
        return self.__Nod_Inf_V

    def node_v_output(self, time):
        return self.extract_nodes(self.node_inf_v, time)[2:-1]

    def cond_iso_thermal_liquid(self, time):
        """"Isothermal hydraulic conductivity of liquid phase"""
        cond = list(zip(*self.node_v_output(time)))[2]

        return [float(x) / 8640000 for x in cond][1:]

    def cond_thermal_liquid(self, time):
        """"Thermal hydraulic conductivity of liquid phase"""
        return list(map(float, list(zip(*self.node_v_output(time)))[3]))[1:]

    def cond_iso_thermal_vapor(self, time):
        """"Isothermal hydraulic conductivity of vapor phase"""

        cond = list(zip(*self.node_v_output(time)))[4]
        return [float(x) / 8640000 for x in cond][1:]

    def cond_thermal_vapor(self, time):
        """"Thermal hydraulic conductivity of vapor phase"""
        return list(map(float, list(zip(*self.node_v_output(time)))[5]))[1:]

    def q_liquid(self, time):
        """"Returns total liquid flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[6]
        return [float(x) / 8640000 for x in q][1:]

    def q_vapor(self, time):
        """"Returns total vapor flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[7]
        return [float(x) / 8640000 for x in q][1:]

    def q_total(self, time):
        """"Returns total water flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[8]
        return [float(x) / 8640000 for x in q][1:]

    def q_vapor_isothermal(self, time):
        """"Returns isothermal vapor flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[9]
        return [float(x) / 8640000 for x in q][1:]

    def q_vapor_thermal(self, time):
        """"Returns thermal vapor flow cm/day"""

        q = list(zip(*self.node_v_output(time)))[10]
        return [float(x) / 8640000 for x in q][1:]


    ###### Atmosphere #######

    @property
    def atm(self):
        return self.__Atm[9]

    def prec(self):
        return float(self.atm[1]) / 86400

    def pot_evaporation(self):
        return float(self.atm[2]) / 86400

    def bot_flux(self):
        return float(self.atm[5])

    def T_s(self):
        return float(self.atm[8])

    def T_bottom(self):
        return float(self.atm[9])


def write_selector_in_file(path):
    """
    """
    # TODO adapt this to parameters of current run

    file_content = f"""Pcp_File_Version=4
*** BLOCK A: BASIC INFORMATION *****************************************
Heading
Welcome to HYDRUS-1D
LUnit  TUnit  MUnit  (indicated units are obligatory for all input data)
cm
days
mmol
lWat   lChem lTemp  lSink lRoot lShort lWDep lScreen lVariabBC lEquil lInverse
 t     f     t      f     f     f      f     t       t         t         f
lSnow  lHP1   lMeteo  lVapor lActiveU lFluxes lIrrig  lDummy  lDummy  lDummy
 f       f       f       t       f       f       f       f       f       f
NMat    NLay  CosAlpha
  1       1       1
*** BLOCK B: WATER FLOW INFORMATION ************************************
MaxIt   TolTh   TolH       (maximum number of iterations and tolerances)
  10    0.001      1
TopInf WLayer KodTop InitCond
 t     f      -1       t
BotInf qGWLF FreeD SeepF KodBot DrainF  hSeep
 f     f     f     f     -1      f      0
         rTop         rBot        rRoot
           0            0            0
    hTab1   hTabN
    1e-006   10000
    Model   Hysteresis
      0          0
   thr     ths    Alfa      n         Ks       l
   0.01 0.35001  0.0518    2.22    1.06272    0.67 
*** BLOCK C: TIME INFORMATION ******************************************
        dt       dtMin       dtMax     DMul    DMul2  ItMin ItMax  MPL
      0.001      1e-005           5     1.3     0.7     3     7   250
      tInit        tMax
          0         250
  lPrintD  nPrintSteps tPrintInterval lEnter
     f           1             1       t
TPrint(1),TPrint(2),...,TPrint(MPL)
          1           2           3           4           5           6 
          7           8           9          10          11          12 
         13          14          15          16          17          18 
         19          20          21          22          23          24 
         25          26          27          28          29          30 
         31          32          33          34          35          36 
         37          38          39          40          41          42 
         43          44          45          46          47          48 
         49          50          51          52          53          54 
         55          56          57          58          59          60 
         61          62          63          64          65          66 
         67          68          69          70          71          72 
         73          74          75          76          77          78 
         79          80          81          82          83          84 
         85          86          87          88          89          90 
         91          92          93          94          95          96 
         97          98          99         100         101         102 
        103         104         105         106         107         108 
        109         110         111         112         113         114 
        115         116         117         118         119         120 
        121         122         123         124         125         126 
        127         128         129         130         131         132 
        133         134         135         136         137         138 
        139         140         141         142         143         144 
        145         146         147         148         149         150 
        151         152         153         154         155         156 
        157         158         159         160         161         162 
        163         164         165         166         167         168 
        169         170         171         172         173         174 
        175         176         177         178         179         180 
        181         182         183         184         185         186 
        187         188         189         190         191         192 
        193         194         195         196         197         198 
        199         200         201         202         203         204 
        205         206         207         208         209         210 
        211         212         213         214         215         216 
        217         218         219         220         221         222 
        223         224         225         226         227         228 
        229         230         231         232         233         234 
        235         236         237         238         239         240 
        241         242         243         244         245         246 
        247         248         249         250 
*** BLOCK E: HEAT TRANSPORT INFORMATION *********************************************************
    Qn      Qo    Disper.    B1          B2          B3          Cn          Co           Cw
    0.6    0.01       0 -1.2706e+016 -6.20464e+016 1.62598e+017 1.43327e+014 1.8737e+014 3.12035e+014 
      tAmpl     tPeriod    Campbell   MeltConst  lDummy  lDummy  lDummy  lDummy  lDummy
          5           1          0        0.43       f       f       f       f       f
      kTopT       TTop      kBotT       TBot
          1         20           1         20
*** END OF INPUT FILE 'SELECTOR.IN' ************************************
    """
    f_path = path / 'SELECTOR.IN'
    with open(f_path, 'w') as writer:
        writer.writelines(file_content)

